#!/usr/bin/env python3
"""serve.py - local-only HTTP server for the nekStab validation gallery.

Replaces `python -m http.server` with a tiny custom server that also
exposes:

  POST /api/resubmit         {case: "<case-path>"}  -> {token}
  GET  /api/jobs/<token>     -> job state + slurm state
  GET  /api/queue            -> live `squeue` snapshot

Security model: binds 127.0.0.1 only. Reachable only via SSH port-forward
or local browser. Cases that can be (re)submitted are restricted to the
allow-list loaded from build_gallery.SECTIONS.

Run:
    python3 validation/serve.py            # default port 8000
    python3 validation/serve.py 8123       # custom port
"""
from __future__ import annotations
import gzip
import json
import subprocess
import sys
import threading
import uuid
from http.server import HTTPServer, SimpleHTTPRequestHandler
from pathlib import Path
from urllib.parse import urlparse

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
EXAMPLE = REPO / 'example'

sys.path.insert(0, str(ROOT))
from build_gallery import SECTIONS  # noqa: E402

ALLOWED: set[str] = {c.path for _, cases in SECTIONS for c in cases}

# In-memory job state: token -> {state, jobid, case, log, slurm_state}
JOBS: dict[str, dict] = {}
JOBS_LOCK = threading.Lock()
JOBS_MAX = 200  # cap memory: prune oldest terminal jobs above this


def _prune_jobs() -> None:
    """Drop oldest terminal jobs when JOBS exceeds JOBS_MAX. Caller holds lock."""
    if len(JOBS) <= JOBS_MAX:
        return
    terminal_tokens = [
        t for t, j in JOBS.items()
        if j.get('state') in ('done', 'failed', 'error')
    ]
    # delete oldest (insertion-order on dict is preserved in Python 3.7+)
    drop = max(0, len(JOBS) - JOBS_MAX)
    for t in terminal_tokens[:drop]:
        del JOBS[t]


def _update_job(token: str, **kwargs) -> None:
    with JOBS_LOCK:
        if token in JOBS:
            JOBS[token].update(kwargs)


def _slurm_state(jobid: str) -> str:
    """Return live state from squeue, falling back to sacct after completion."""
    if not jobid:
        return 'UNKNOWN'
    try:
        r = subprocess.run(
            ['squeue', '-j', jobid, '-h', '-o', '%T'],
            capture_output=True, text=True, timeout=5,
        )
        s = r.stdout.strip()
        if s:
            return s
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return 'UNKNOWN'
    try:
        r = subprocess.run(
            ['sacct', '-j', jobid, '-n', '-o', 'State', '-X', '-P'],
            capture_output=True, text=True, timeout=5,
        )
        toks = r.stdout.strip().split()
        if toks:
            return toks[0]
    except (subprocess.TimeoutExpired, FileNotFoundError):
        pass
    return 'COMPLETED?'


def _read_casename(case_dir: Path) -> str | None:
    sess = case_dir / 'SESSION.NAME'
    if not sess.exists():
        return None
    try:
        return sess.read_text().splitlines()[0].strip()
    except OSError:
        return None


def _run_mks(token: str, case_dir: Path, casename: str) -> bool:
    """Run mks (twice if first attempt fails due to NEKSTAB.inc rewrite)."""
    _update_job(token, state='compiling', log=f'mks {casename} (1/2)')
    for attempt in (1, 2):
        try:
            r = subprocess.run(
                ['mks', casename],
                cwd=str(case_dir),
                capture_output=True, text=True, timeout=600,
            )
        except subprocess.TimeoutExpired:
            _update_job(token, state='error', log='mks timed out (>10 min)')
            return False
        except FileNotFoundError:
            _update_job(token, state='error', log='mks not on PATH')
            return False
        if r.returncode == 0:
            return True
        if attempt == 1:
            _update_job(token, log=f'mks {casename} (2/2 retry)')
            continue
        _update_job(token, state='error', log=f'mks failed: {r.stderr[-400:]}')
        return False
    return False


def _run_sbatch(token: str, case_dir: Path) -> bool:
    """Submit run.local.slurm. Captures the parsable job ID."""
    _update_job(token, state='submitting', log='sbatch run.local.slurm')
    try:
        r = subprocess.run(
            ['sbatch', '--parsable', 'run.local.slurm'],
            cwd=str(case_dir),
            capture_output=True, text=True, timeout=30,
        )
    except subprocess.TimeoutExpired:
        _update_job(token, state='error', log='sbatch timed out')
        return False
    except FileNotFoundError:
        _update_job(token, state='error', log='sbatch not on PATH')
        return False
    if r.returncode != 0:
        _update_job(token, state='error', log=f'sbatch failed: {r.stderr[-400:]}')
        return False
    jobid = r.stdout.strip().split(';')[0]
    _update_job(token, state='queued', jobid=jobid, log='')
    return True


def _action_worker(token: str, case: str, action: str) -> None:
    """Background thread for one of: 'recompile', 'resubmit', 'recompile_submit'."""
    case_dir = EXAMPLE / case
    slurm = case_dir / 'run.local.slurm'
    casename = _read_casename(case_dir)

    if action in ('recompile', 'recompile_submit') and not casename:
        _update_job(token, state='error', log='no SESSION.NAME (need it for mks)')
        return
    if action in ('resubmit', 'recompile_submit') and not slurm.exists():
        _update_job(token, state='error', log='no run.local.slurm')
        return

    if action == 'recompile':
        if _run_mks(token, case_dir, casename):  # type: ignore[arg-type]
            _update_job(token, state='done', log='nek5000 binary built')
        return
    if action == 'resubmit':
        _run_sbatch(token, case_dir)
        return
    # combined
    if _run_mks(token, case_dir, casename):  # type: ignore[arg-type]
        _run_sbatch(token, case_dir)


def _live_queue() -> dict:
    """Snapshot of pending+running jobs in the local queue."""
    try:
        r = subprocess.run(
            ['squeue', '-h', '-o', '%i|%j|%T|%M|%D|%R'],
            capture_output=True, text=True, timeout=5,
        )
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        return {'error': str(e)}
    if r.returncode != 0:
        return {'error': r.stderr.strip()[:200] or f'rc={r.returncode}'}
    jobs = []
    for line in r.stdout.strip().splitlines():
        parts = line.split('|')
        if len(parts) < 6:
            continue
        jobs.append({
            'jobid': parts[0],
            'name': parts[1],
            'state': parts[2],
            'time': parts[3],
            'nodes': parts[4],
            'reason': parts[5],
        })
    return {'jobs': jobs}


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        # Serve the repo root so the gallery (validation/index.html) and the
        # in-place plots it references ("../example/.../plot_*.png") both
        # resolve under one root, matching the Pages artifact layout.
        super().__init__(*args, directory=str(REPO), **kwargs)

    def log_message(self, format, *args):  # noqa: A002 - base API name
        sys.stderr.write(f'[{self.log_date_time_string()}] {format % args}\n')

    def _json(self, status: int, payload: dict) -> None:
        body = json.dumps(payload).encode()
        accepts_gzip = 'gzip' in (self.headers.get('Accept-Encoding') or '')
        # only worth gzipping payloads > ~512 bytes
        if accepts_gzip and len(body) > 512:
            body = gzip.compress(body)
            self.send_response(status)
            self.send_header('Content-Type', 'application/json')
            self.send_header('Content-Encoding', 'gzip')
        else:
            self.send_response(status)
            self.send_header('Content-Type', 'application/json')
        self.send_header('Cache-Control', 'no-store')
        self.send_header('Content-Length', str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    _ACTIONS = {'recompile', 'resubmit', 'recompile_submit'}
    _ALLOWED_ORIGINS = {
        f'http://127.0.0.1:{p}' for p in (8000, 8001, 8080)
    } | {
        f'http://localhost:{p}' for p in (8000, 8001, 8080)
    } | {'null'}  # 'null' allows file:// for local file viewing

    def _is_origin_allowed(self) -> bool:
        """CSRF defense (#1): require Origin header to match our bind."""
        origin = self.headers.get('Origin') or self.headers.get('Referer', '')
        if not origin:
            # no Origin = old browser / same-origin GET. Accept only for GET-shaped requests.
            # For POST we require an Origin.
            return False
        # match prefix to tolerate referer suffixes
        return any(origin == o or origin.startswith(o + '/') for o in self._ALLOWED_ORIGINS)

    def _case_dir_safe(self, case: str) -> Path | None:
        """Path traversal defense (#2): resolve and verify it's under EXAMPLE."""
        try:
            example_real = EXAMPLE.resolve(strict=False)
            candidate = (EXAMPLE / case).resolve(strict=False)
            candidate.relative_to(example_real)
        except (ValueError, OSError):
            return None
        return candidate

    def _has_active_job(self, case: str) -> bool:
        """Dedup (#3): is there already a non-terminal job for this case?"""
        with JOBS_LOCK:
            for job in JOBS.values():
                if job.get('case') == case and job.get('state') not in (
                    'done', 'failed', 'error'
                ):
                    return True
        return False

    def do_POST(self):
        url = urlparse(self.path)

        # CSRF check (#1)
        if not self._is_origin_allowed():
            return self._json(403, {'error': 'origin not allowed'})

        # Map URL -> action. Keep /api/resubmit + /api/recompile + /api/run.
        action: str | None = None
        if url.path == '/api/resubmit':
            action = 'resubmit'
        elif url.path == '/api/recompile':
            action = 'recompile'
        elif url.path == '/api/run':
            action = None  # read from body
        else:
            return self._json(404, {'error': 'not found'})

        length = int(self.headers.get('Content-Length', '0') or '0')
        if length > 4096:  # body should never exceed a few hundred bytes
            return self._json(413, {'error': 'request too large'})
        try:
            data = json.loads(self.rfile.read(length) or b'{}')
        except (json.JSONDecodeError, ValueError):
            return self._json(400, {'error': 'invalid json'})
        case = data.get('case', '')
        if action is None:
            action = data.get('action', 'resubmit')
        if action not in self._ACTIONS:
            return self._json(400, {'error': f'unknown action: {action}'})
        if not isinstance(case, str) or case not in ALLOWED:
            return self._json(403, {'error': f'case not in catalog: {case}'})

        # Path traversal defense (#2)
        case_dir = self._case_dir_safe(case)
        if case_dir is None or not case_dir.is_dir():
            return self._json(404, {'error': f'case dir missing: {case}'})

        # Dedup (#3)
        if self._has_active_job(case):
            return self._json(409, {'error': 'a job for this case is already running'})

        # Stronger token entropy (#4)
        token = uuid.uuid4().hex
        with JOBS_LOCK:
            JOBS[token] = {
                'state': 'starting', 'case': case, 'action': action,
                'jobid': None, 'log': '',
            }
            _prune_jobs()
        t = threading.Thread(
            target=_action_worker, args=(token, case, action), daemon=True,
        )
        t.start()
        return self._json(200, {'token': token, 'action': action})

    def do_GET(self):
        url = urlparse(self.path)
        if url.path == '/':
            self.send_response(302)
            self.send_header('Location', '/validation/index.html')
            self.end_headers()
            return
        if url.path.startswith('/api/jobs/'):
            token = url.path[len('/api/jobs/'):]
            with JOBS_LOCK:
                job = dict(JOBS.get(token, {}))
            if not job:
                return self._json(404, {'error': 'unknown token'})
            if job.get('jobid') and job.get('state') in ('queued', 'running'):
                state = _slurm_state(job['jobid'])
                job['slurm_state'] = state
                if state == 'PENDING':
                    job['state'] = 'queued'
                elif state == 'RUNNING':
                    job['state'] = 'running'
                elif state == 'COMPLETED':
                    job['state'] = 'done'
                elif state in ('FAILED', 'CANCELLED', 'TIMEOUT', 'NODE_FAIL'):
                    job['state'] = 'failed'
                    job['log'] = f'slurm state: {state}'
                _update_job(token, **job)
            return self._json(200, job)
        if url.path == '/api/queue':
            return self._json(200, _live_queue())
        return super().do_GET()


def main():
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
    addr = ('127.0.0.1', port)
    print(f'nekStab gallery on http://{addr[0]}:{addr[1]}/validation/index.html')
    print(f'  serving:        {REPO}')
    print(f'  allowed cases:  {len(ALLOWED)}')
    print(f'  example root:   {EXAMPLE}')
    print(f'  PATH has mks:   {bool(_which("mks"))}')
    print(f'  PATH has sbatch:{bool(_which("sbatch"))}')
    HTTPServer(addr, Handler).serve_forever()


def _which(cmd: str) -> str | None:
    from shutil import which
    return which(cmd)


if __name__ == '__main__':
    main()
