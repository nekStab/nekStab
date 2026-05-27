#!/usr/bin/env python3
"""build_gallery.py - case-catalog-driven HTML gallery for nekStab validation.

Defines the 24 validated v2.0 cases authoritatively and the 15 deferred
ones, then walks validation/figures/ to attach a PNG to each. Result is
a dashboard:

  - validated cases with a matching PNG  -> green check + thumbnail
  - validated cases with stale Re/Ra     -> amber "stale" badge
  - validated cases with no PNG          -> grey "pending" placeholder
  - deferred cases at the bottom (collapsed)

Run:  python validation/build_gallery.py
Open: xdg-open validation/index.html
"""
from __future__ import annotations
from dataclasses import dataclass
from datetime import datetime
from functools import lru_cache
from pathlib import Path
import html
import re
import subprocess

_PNG_NOT_SET = object()  # sentinel for _render_card png_override default

try:
    from validation.catalog import (
        families as catalog_families, cases_for, all_cases as catalog_all_cases,
        family_can, CaseStatus, FlowFamily, MethodCase
    )
except ModuleNotFoundError:
    from catalog import (
        families as catalog_families, cases_for, all_cases as catalog_all_cases,
        family_can, CaseStatus, FlowFamily, MethodCase
    )

ROOT = Path(__file__).resolve().parent
FIG_DIR = ROOT / 'figures'
OUT = ROOT / 'index.html'
EXAMPLE = ROOT.parent / 'example'


@dataclass(frozen=True)
class Case:
    path: str            # 'cylinder/baseflow/newton'
    re_tag: str          # 'Re100' or 'Ra500'
    mode: str            # 'Newton-GMRES (2.0)'
    expected: str        # short physical expectation
    status: str          # 'validated' | 'deferred'
    fig_key: str         # collect_plots.py filename prefix (without .png)


# Sections in display order. Each section has a heading + list[Case].
SECTIONS: list[tuple[str, list[Case]]] = [
    ('Cylinder · Baseflow', [
        Case('cylinder/baseflow/boostconv',     'Re100', 'BoostConv (1.2)',          'converged steady state',                                  'validated', 'cylinder_baseflow_boostconv'),
        Case('cylinder/baseflow/newton',        'Re100', 'Newton-GMRES (2.0)',       'converged steady state',                                  'validated', 'cylinder_baseflow_newton'),
        Case('cylinder/baseflow/newton_dyn',    'Re100', 'Newton + ifdyntol (2.0)',  'converged steady state, fewer GMRES iters',               'validated', 'cylinder_baseflow_newton_dyn'),
        Case('cylinder/baseflow/sfd',           'Re50',  'SFD baseline (1.1)',       'converged steady state',                                  'validated', 'cylinder_baseflow_sfd'),
        Case('cylinder/baseflow/sfd_dyn',       'Re50',  'SFD + dyn-tol (1.1)',      'converged steady state',                                  'validated', 'cylinder_baseflow_sfd_dyn'),
        Case('cylinder/baseflow/sfd_dyn_oifs',  'Re50',  'SFD + OIFS (1.1)',         'converged steady state',                                  'validated', 'cylinder_baseflow_sfd_dyn_oifs'),
    ]),
    ('Cylinder · Stability', [
        Case('cylinder/stability/direct',       'Re100', 'Direct LNSE (3.1)',        'leading eigenvalue near sigma = 0.125 + 0.734i',          'validated', 'cylinder_stability_direct'),
        Case('cylinder/stability/adjoint',      'Re100', 'Adjoint LNSE (3.2)',       'matches direct spectrum',                                 'validated', 'cylinder_stability_adjoint'),
        Case('cylinder/stability/animate_modes','Re100', 'Animate modes (4.50)',     'snapshots per period',                                    'validated', 'cylinder_stability_animate_modes'),
    ]),
    ('Cylinder · Postproc', [
        Case('cylinder/postproc/sensitivity_budget_wavemaker', 'Re100', 'Wavemaker / budget (4.0)', 'wm, sr/si, pr/pi, tr/ti fields', 'validated', 'cylinder_postproc_sensitivity_budget_wavemaker'),
        Case('cylinder/postproc/steady_force_sensitivity',     'Re100', 'Steady force (4.41)',     'sensitivity field fsr',           'validated', 'cylinder_postproc_steady_force_sensitivity'),
    ]),
    ('Cylinder · DNS / smoke', [
        Case('cylinder/dns',             'Re100', 'DNS (0.0)',                'vortex shedding past t=100', 'validated', 'cylinder_dns'),
        Case('cylinder/ci_test',         'Re100', 'DNS 10-step smoke (0.0)',  '10-step build smoke',         'validated', 'cylinder_ci_test'),
        Case('cylinder/moving_cylinder', 'Re100', 'DNS forced osc (0.0)',     'moving-mesh DNS subcase',     'validated', 'cylinder_moving_cylinder'),
    ]),
    ('Other geometries · Baseflow', [
        Case('lid_driven',            'Re3600', 'Newton (2.0)',             'converged steady state',                'validated', 'lid_driven'),
        Case('poiseuille_OTD',        'Re5000', 'OTD smoke',                'builds and runs to endTime=500',        'validated', 'poiseuille'),
        Case('thermosyphon/baseflow',   'Ra500',  'Newton + thermal (2.0)',   'converged thermal steady state',        'validated', 'thermosyphon_baseflow'),
        Case('naca0012',              'Re2000', 'Newton (2.0)',             'converged steady state',                'validated', 'naca0012'),
        Case('flip_flop/baseflow',    'Re62',   'Newton-UPO (2.1)',         'converged UPO',                         'validated', 'flip_flop_baseflow'),
        Case('tpjet/baseflow/newton', 'Re1900', 'Forced UPO (2.2)',         'converged forced UPO',                  'validated', 'tpjet_baseflow_newton'),
    ]),
    ('Other geometries · Stability', [
        Case('thermosyphon/stability/direct',         'Ra500',  'Direct + thermal (3.1)', 'leading real eigenvalue (pitchfork above Ra_c=494)', 'validated', 'thermosyphon_stability_direct'),
        Case('flip_flop/stability/direct_Floquet',  'Re62',   'Floquet direct (3.11)',  'leading multiplier just above unit circle',          'validated', 'flip_flop_stability_direct_Floquet'),
        Case('tpjet/stability/direct_Floquet',      'Re1900', 'Floquet direct (3.11)',  'leading multiplier above unit circle (Re_c~1371)',   'validated', 'tpjet_stability_direct_Floquet'),
    ]),
    ('Deferred for v2.1', [
        Case('cylinder/baseflow/newton_dyn_temp',         'Re100', 'Newton + thermal',         'thermal seed needed',         'deferred', 'cylinder_baseflow_newton_dyn_temp'),
        Case('cylinder/baseflow/newton_dyn_upo',          'Re100', 'Newton-UPO',               'UPO did not converge on workstation budget', 'deferred', 'cylinder_baseflow_newton_upo'),
        Case('cylinder/stability/direct_Floquet',         'Re100', 'Floquet direct',           'needs converged UPO',         'deferred', 'cylinder_stability_direct_Floquet'),
        Case('cylinder/stability/adjoint_Floquet',        'Re100', 'Floquet adjoint',          'needs converged UPO',         'deferred', 'cylinder_stability_adjoint_Floquet'),
        Case('cylinder/stability/animate_modes_with_UPO', 'Re100', 'Animate w/ UPO',           'needs converged UPO',         'deferred', 'cylinder_stability_animate_modes_with_UPO'),
        Case('cylinder/otd',                              'Re100', 'OTD modes',                'NaN-at-init on 2128 mesh',    'deferred', 'cylinder_otd'),
        Case('cylinder/modal',                            'Re100', 'POD/DMD/SPOD',             'needs DNS snapshots',         'deferred', 'cylinder_modal'),
        Case('cylinder/RANS',                             'Re100', 'RANS baseflow',            'turbulence-model integration', 'deferred', 'cylinder_RANS'),
        Case('back_fstep/baseflow',                       'Re500', 'SFD baseline',             'CFL emergency exit during sweep', 'deferred', 'back_fstep_baseflow'),
        Case('back_fstep/transient_growth',               'Re500', 'Transient growth (3.3)',   'depends on baseflow',         'deferred', 'back_fstep_transient_growth'),
        Case('cubic_cavity',                              'Re400', 'Newton',                   'mesh issues',                 'deferred', 'cubic_cavity'),
        Case('cubic_cavity_upo',                          'Re400', 'Newton-UPO',               'UPO budget exceeded',         'deferred', 'cubic_cavity_upo'),
        Case('tpjet/baseflow/tdf',                        'Re1900', 'TDF (1.3)',               'period harmonics tracking',   'deferred', 'tpjet_baseflow_tdf'),
        Case('slot_FST',                                  'Re495', 'DNS with FST',             'IC bootstrap not in v2.0 sweep', 'deferred', 'slot_FST'),
        Case('blasius',                                   'Re1000','OTD on Blasius',           'not in v2.0 sweep',           'deferred', 'blasius'),
        Case('poiseuille_RANS',                           'Re5000','RANS Poiseuille',          'not in v2.0 sweep',           'deferred', 'poiseuille_RANS'),
    ]),
]


@lru_cache(maxsize=None)
def _case_stats(case_path: str) -> dict:
    """Read mesh elements (from .re2), MPI ranks (from slurm), wall time
    (from logfile) for a case. Returns {} when any source missing."""
    case_dir = EXAMPLE / case_path
    out: dict = {}
    if not case_dir.is_dir():
        return out

    # mesh elements: first 132 bytes of .re2 are ASCII header, second token = nelgt
    re2 = next(case_dir.glob('*.re2'), None)
    if re2 is not None:
        try:
            header = re2.read_bytes()[:132].decode('ascii', errors='replace')
            toks = header.split()
            if len(toks) >= 2 and toks[1].isdigit():
                out['nelg'] = int(toks[1])
        except OSError:
            pass

    # MPI ranks from --ntasks=N in slurm script
    slurm = case_dir / 'run.local.slurm'
    if slurm.exists():
        try:
            txt = slurm.read_text(errors='replace')
            m = re.search(r'--ntasks[= ](\d+)', txt)
            if m:
                out['ntasks'] = int(m.group(1))
        except OSError:
            pass
        out['has_slurm'] = True
    else:
        out['has_slurm'] = False

    # wall time: 'total elapsed time : 6.45376E+03 sec' near end of logfile
    log = case_dir / 'logfile'
    if log.exists():
        try:
            # tail-style read of last 8 KB
            with log.open('rb') as f:
                f.seek(0, 2)
                size = f.tell()
                f.seek(max(0, size - 8192))
                tail = f.read().decode('utf-8', errors='replace')
            m = re.search(r'total elapsed time\s*:\s*([0-9.E+\-]+)\s*sec', tail)
            if m:
                out['elapsed_sec'] = float(m.group(1))
        except OSError:
            pass

    # case name from SESSION.NAME (first line)
    sess = case_dir / 'SESSION.NAME'
    if sess.exists():
        try:
            out['casename'] = sess.read_text().splitlines()[0].strip()
        except OSError:
            pass

    return out


def _format_elapsed(sec: float) -> str:
    if sec < 60:
        return f'{sec:.0f} s'
    if sec < 3600:
        return f'{sec / 60:.1f} min'
    return f'{sec / 3600:.2f} h'


def _classify(case: Case) -> tuple[str, str]:
    """Return (geometry, mode) for grouping. Geometry = top dir or flat name.
    Mode = one of: baseflow, stability, postproc, dns, otd, modal, rans, other.
    """
    parts = case.path.split('/')
    geom = parts[0]
    if len(parts) >= 2:
        seg = parts[1].lower()
        if seg in ('baseflow', 'stability', 'postproc', 'modal'):
            return geom, seg
        if seg in ('otd', 'rans'):
            return geom, seg
        if seg in ('dns', 'ci_test', 'moving_cylinder'):
            return geom, 'dns'
        if 'growth' in seg:
            return geom, 'stability'
        return geom, seg
    name = parts[0].lower()
    if 'otd' in name:        return geom, 'otd'
    if 'rans' in name:       return geom, 'rans'
    if 'cavity' in name:     return geom, 'baseflow'
    if 'slot' in name:       return geom, 'dns'
    if name == 'blasius':    return geom, 'otd'
    return geom, 'baseflow'


# Geometry display ordering — richer / canonical first
GEOM_ORDER = [
    'cylinder', 'flip_flop', 'tpjet', 'thermosyphon',
    'naca0012', 'lid_driven', 'poiseuille_OTD',
    'back_fstep', 'cubic_cavity', 'cubic_cavity_upo',
    'blasius', 'slot_FST', 'poiseuille_RANS',
]
# Within a geometry, modes display in this order. Baseflow first (always open).
MODE_ORDER = ['baseflow', 'stability', 'postproc', 'dns', 'otd', 'modal', 'rans', 'other']
MODE_LABEL = {
    'baseflow':  'Baseflow',
    'stability': 'Stability',
    'postproc':  'Postproc',
    'dns':       'DNS',
    'otd':       'OTD',
    'modal':     'Modal (POD/DMD/SPOD)',
    'rans':      'RANS',
    'other':     'Other',
}
GEOM_LABEL = {
    'cylinder':         'Cylinder',
    'flip_flop':        'Flip-flop',
    'tpjet':            'TPJet',
    'thermosyphon':       'Thermosyphon',
    'naca0012':         'NACA0012',
    'lid_driven':       'Lid-driven cavity',
    'poiseuille_OTD':   'Poiseuille (OTD)',
    'back_fstep':       'Backward-facing step',
    'cubic_cavity':     'Cubic cavity',
    'cubic_cavity_upo': 'Cubic cavity (UPO)',
    'blasius':          'Blasius',
    'slot_FST':         'Slot jet + FST',
    'poiseuille_RANS':  'Poiseuille (RANS)',
}

FAMILY_ORDER = [
    'cylinder_re100', 'cylinder_re180', 'flip_flop', 'tpjet', 'thermosyphon',
    'back_fstep', 'lid_driven_cavity', 'naca0012', 'poiseuille',
    'moving_cylinder', 'cubic_cavity', 'slot_FST', 'blasius',
    'cylinder_re1m',
]

LANE_ORDER = ['dns', 'baseflow', 'transient_growth', 'stability_direct', 'direct',
               'stability_adjoint', 'animation', 'floquet_direct', 'floquet_adjoint',
               'wavemaker', 'otd', 'modal', 'rans']

BLOCKED_STATUSES = {
    CaseStatus.BLOCKED,
    CaseStatus.NEEDS_DNS_SEED,
    CaseStatus.NEEDS_SCALAR_CHECKPOINT,
}
PARTIAL_STATUSES = {
    CaseStatus.PARTIAL,
    CaseStatus.LOCAL_SMOKE,
    CaseStatus.AMBIGUOUS,
    CaseStatus.STALE,
    CaseStatus.MISSING_IMAGE,
}


def _case_subpath(current_path: str) -> str:
    return current_path.removeprefix('example/').strip('/')


def _re_tag(mc: MethodCase) -> str:
    params = mc.target_parameters
    if params.Re is not None:
        return f'Re{int(params.Re)}'
    if mc.flow_family == 'thermosyphon' or params.geometry == 'thermosyphon':
        return 'Ra500'
    return 'Re0'


def catalog_to_case(mc: MethodCase) -> Case:
    if mc.status == CaseStatus.VALIDATED:
        status = 'validated'
    elif mc.status in BLOCKED_STATUSES:
        status = 'blocked'
    elif mc.status == CaseStatus.DEFERRED:
        status = 'deferred'
    else:
        status = 'deferred'
    return Case(
        path=_case_subpath(mc.current_path),
        re_tag=_re_tag(mc),
        mode=mc.mode_name or '',
        expected=mc.expected_behavior,
        status=status,
        fig_key=mc.case_id.replace('/', '_'),
    )


def _mode_text(mc: MethodCase) -> str:
    if mc.mode_name and mc.legacy_uparam01:
        return f'{mc.mode_name} ({mc.legacy_uparam01})'
    if mc.mode_name:
        return mc.mode_name
    if mc.legacy_uparam01:
        return f'uparam01 {mc.legacy_uparam01}'
    return mc.method_lane.replace('_', ' ')


def _status_group(status: CaseStatus) -> str:
    if status == CaseStatus.VALIDATED:
        return 'validated'
    if status == CaseStatus.DEFERRED:
        return 'deferred'
    if status in BLOCKED_STATUSES:
        return 'blocked'
    if status in PARTIAL_STATUSES:
        return 'partial'
    return 'partial'


def _status_counts(cases: tuple[MethodCase, ...] | list[MethodCase]) -> dict[str, int]:
    counts = {'validated': 0, 'deferred': 0, 'blocked': 0, 'partial': 0}
    for mc in cases:
        counts[_status_group(mc.status)] += 1
    return counts


def _count_chips(counts: dict[str, int]) -> str:
    bits: list[str] = []
    if counts.get('validated'):
        bits.append(f'<span class="c-ok">{counts["validated"]} validated</span>')
    if counts.get('partial'):
        bits.append(f'<span class="c-pending">{counts["partial"]} partial</span>')
    if counts.get('deferred'):
        bits.append(f'<span class="c-pending">{counts["deferred"]} deferred</span>')
    if counts.get('blocked'):
        bits.append(f'<span class="c-stale">{counts["blocked"]} blocked</span>')
    return ''.join(bits)


def _ordered_families() -> list[FlowFamily]:
    family_by_id = {f.family_id: f for f in catalog_families()}
    return (
        [family_by_id[fid] for fid in FAMILY_ORDER if fid in family_by_id]
        + [family_by_id[fid] for fid in sorted(family_by_id) if fid not in FAMILY_ORDER]
    )


def _git_version() -> str:
    """Best-effort nekStab version: `git describe` w/ tags + dirty marker."""
    repo = ROOT.parent
    try:
        desc = subprocess.check_output(
            ['git', '-C', str(repo), 'describe', '--tags', '--always', '--dirty'],
            stderr=subprocess.DEVNULL,
        ).decode().strip()
        branch = subprocess.check_output(
            ['git', '-C', str(repo), 'rev-parse', '--abbrev-ref', 'HEAD'],
            stderr=subprocess.DEVNULL,
        ).decode().strip()
        return f'{desc} ({branch})' if branch and branch != 'HEAD' else desc
    except (subprocess.CalledProcessError, FileNotFoundError):
        return 'unknown'


@lru_cache(maxsize=None)
def _find_png(case: Case) -> Path | None:
    """Best-fit PNG match in validation/figures/.
    Returns the file that begins with the case's fig_key and (if possible)
    also matches the expected Re/Ra tag. Falls back to any prefix match.
    """
    exact = FIG_DIR / f'{case.fig_key}_{case.re_tag}.png'
    if exact.exists():
        return exact
    prefix_match = sorted(FIG_DIR.glob(f'{case.fig_key}_*.png'))
    if prefix_match:
        return prefix_match[0]
    return None


def _find_artifact_image(artifact) -> "Path | None":
    """Resolve artifact.path (e.g. 'validation/figures/foo.png') to an absolute Path.
    Supports .png and .gif. Returns None if the file does not exist.
    """
    p = ROOT.parent / artifact.path
    if p.exists():
        return p
    return None


def _panel_caption_from_path(art_path: str, case_id: str) -> str:
    """Derive a human-readable panel descriptor from the artifact filename.

    Strips the _Re<N>/_Ra<N> suffix and extension, then strips the longest
    common underscore-token prefix shared with the case_id slug. Title-cases
    the remainder (e.g. 'field', 'direct_mode' -> 'Field', 'Direct Mode').
    Returns empty string for single-panel cases where no suffix remains.
    """
    stem = Path(art_path).stem  # e.g. 'cylinder_postproc_wavemaker_field_Re50'
    # Strip trailing _Re<N> or _Ra<N> (case-insensitive on Ra/Re)
    stem = re.sub(r'_R[aAeE]\d+$', '', stem)
    # The case_id slug (e.g. 'cylinder_postproc_sensitivity_budget_wavemaker')
    slug = case_id.replace('/', '_')
    # Find longest common prefix on underscore-split tokens
    stem_parts = stem.split('_')
    slug_parts = slug.split('_')
    common = 0
    for sp, cp in zip(stem_parts, slug_parts):
        if sp == cp:
            common += 1
        else:
            break
    panel_parts = stem_parts[common:]
    panel = ' '.join(panel_parts).strip()
    return panel.title() if panel else ''


def _is_stale(case: Case, png: Path) -> bool:
    """A PNG is stale if its filename's Re/Ra tag doesn't match the case's."""
    m = re.search(r'_(R[ae]\d+)\.png$', png.name)
    return bool(m) and m.group(1) != case.re_tag


def _evidence_role_css(role_value: str) -> str:
    mapping = {
        'reference_image': 'ev-image',
        'dns_time_history': 'ev-dns',
        'residual_history': 'ev-residual',
        'spectrum': 'ev-spectrum',
        'eigenvalue_table': 'ev-spectrum',
        'mode_shape': 'ev-image',
        'animation': 'ev-animation',
        'checkpoint_stats': 'ev-checkpoint',
        'provenance': 'ev-provenance',
        'log': 'ev-log',
    }
    return mapping.get(role_value, 'ev-other')


def _evidence_role_label(role_value: str) -> str:
    mapping = {
        'reference_image': 'IMG',
        'dns_time_history': 'DNS-HX',
        'residual_history': 'RESID',
        'spectrum': 'SPEC',
        'eigenvalue_table': 'EIG',
        'mode_shape': 'MODE',
        'animation': 'ANIM',
        'checkpoint_stats': 'CKPT',
        'provenance': 'PROV',
        'log': 'LOG',
    }
    return mapping.get(role_value, role_value.upper()[:6])


CSS = """
:root {
    /* surfaces */
    --bg:         #0e1116;
    --bg-card:    #1a1f2a;
    --bg-defer:   #1a1410;
    --bg-subtle:  rgba(255, 255, 255, 0.025);
    /* text */
    --fg:         #d4d8e0;
    --fg-strong:  #ffffff;
    --fg-mute:    #8088a0;
    --fg-dim:     #5a6478;
    /* interaction accent — used for links, hover, focus */
    --accent:     #66ccff;
    --accent-2:   rgba(102, 204, 255, 0.5);
    --accent-bg:  rgba(102, 204, 255, 0.1);
    /* status semantics */
    --success:    #6ee08a;     /* alias --ok */
    --warning:    #e0b86e;     /* alias --stale */
    --danger:     #e07070;
    --info:       #6e7be0;
    --neutral:    #555c6b;     /* alias --pending */
    /* compatibility aliases used throughout */
    --ok:         var(--success);
    --stale:      var(--warning);
    --pending:    var(--neutral);
    /* lines */
    --border:     #2a3140;
    --border-faint: rgba(255, 255, 255, 0.04);
    /* motion curves */
    --ease-snap:  cubic-bezier(0.2, 0, 0, 1);
    --ease-out:   cubic-bezier(0.16, 1, 0.3, 1);
}
* { box-sizing: border-box; }
body {
    background: var(--bg);
    color: var(--fg);
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    margin: 0;
    padding: 2em 3em;
    line-height: 1.45;
    -webkit-font-smoothing: antialiased;
    -moz-osx-font-smoothing: grayscale;
}
h1 {
    font-size: 1.7em;
    font-weight: 700;
    letter-spacing: -0.015em;
    margin: 0 0 0.15em;
    color: var(--fg-strong);
    text-wrap: balance;
}
h2 {
    font-size: 1.05em;
    font-weight: 600;
    letter-spacing: -0.005em;
    color: var(--fg-strong);
    border-bottom: 1px solid var(--border);
    padding-bottom: 0.5em;
    margin: 2.6em 0 1em;
    display: flex;
    align-items: baseline;
    gap: 0.7em;
    text-wrap: balance;
}
h2 .counts {
    font-size: 0.78em;
    color: var(--fg-mute);
    font-weight: 500;
    font-variant-numeric: tabular-nums;
}
/* status chips — work in any container (h2 counts, details summary, etc.) */
.c-ok,
.c-stale,
.c-pending {
    display: inline-block;
    padding: 0.1em 0.55em;
    border-radius: 999px;
    margin-right: 0.25em;
    font-size: 0.95em;
}
.c-ok      { background: rgba(110, 224, 138, 0.14); color: var(--success); }
.c-stale   { background: rgba(224, 184, 110, 0.14); color: var(--warning); }
.c-pending { background: rgba(255, 255, 255, 0.05); color: var(--fg-mute); }
:focus-visible {
    outline: 2px solid var(--accent);
    outline-offset: 2px;
    border-radius: 4px;
}
::selection {
    background: var(--accent-bg);
    color: var(--fg-strong);
}
/* webkit scrollbar — subtle, matches palette */
::-webkit-scrollbar { width: 12px; height: 12px; }
::-webkit-scrollbar-track { background: transparent; }
::-webkit-scrollbar-thumb {
    background: var(--border);
    border-radius: 999px;
    border: 3px solid var(--bg);
}
::-webkit-scrollbar-thumb:hover { background: var(--fg-dim); }
h2 .counts {
    font-size: 0.75em;
    color: var(--fg-mute);
    font-weight: normal;
}
.meta {
    color: var(--fg-mute);
    margin: 0 0 1.2em;
    font-size: 0.85em;
    font-variant-numeric: tabular-nums;
}

.page-head { margin: 0 0 1.2em; }
.meta code {
    font-family: ui-monospace, SF Mono, Monaco, Cascadia Code, monospace;
    color: var(--accent);
    background: var(--accent-bg);
    padding: 0.1em 0.45em;
    border-radius: 4px;
    font-size: 0.95em;
}
.legend {
    display: flex;
    gap: 1.5em;
    font-size: 0.85em;
    color: var(--fg-mute);
    margin: 1em 0 2em;
    flex-wrap: wrap;
}
.legend span::before {
    display: inline-block;
    width: 0.9em;
    height: 0.9em;
    margin-right: 0.4em;
    vertical-align: -0.1em;
    border-radius: 50%;
    content: '';
}
.legend .ok::before      { background: var(--ok); }
.legend .stale::before   { background: var(--stale); }
.legend .pending::before { background: var(--pending); }

/* Section nav chips (I3) */
.section-nav {
    display: flex;
    flex-wrap: wrap;
    gap: 0.5em;
    margin: 0 0 1.2em;
    padding: 0.6em 0.8em;
    background: rgba(255, 255, 255, 0.02);
    border: 1px solid rgba(255, 255, 255, 0.04);
    border-radius: 10px;
    position: sticky;
    top: 0.5em;
    z-index: 50;
    backdrop-filter: blur(8px);
    -webkit-backdrop-filter: blur(8px);
}
.section-nav .nav-chip {
    color: var(--fg-mute);
    text-decoration: none;
    font-size: 0.82em;
    padding: 0.35em 0.7em;
    border-radius: 6px;
    border: 1px solid rgba(255, 255, 255, 0.06);
    transition: background 0.15s ease-out, color 0.15s ease-out, border-color 0.15s ease-out, transform 0.1s ease-out;
    display: inline-flex;
    align-items: center;
    gap: 0.45em;
    white-space: nowrap;
}
.section-nav .nav-chip:hover {
    color: var(--accent);
    border-color: rgba(102, 204, 255, 0.4);
    background: rgba(102, 204, 255, 0.08);
}
.section-nav .nav-chip:active {
    transform: scale(0.97);
}
.section-nav .nav-chip .count {
    color: var(--fg-mute);
    font-size: 0.85em;
    font-variant-numeric: tabular-nums;
    background: rgba(255, 255, 255, 0.04);
    padding: 0 0.45em;
    border-radius: 999px;
}
.section-nav .nav-chip:hover .count {
    background: rgba(102, 204, 255, 0.15);
    color: var(--accent);
}
/* Section heading scroll-offset (so sticky nav doesn't cover) */
h2[id^="sec-"], details[id^="sec-"] {
    scroll-margin-top: 4.5em;
}

/* Geometry sections (case-first restructure) */
section.geom { margin: 2.6em 0 0; }
section.geom > h2 { margin: 0 0 0.8em; scroll-margin-top: 4.5em; }
.mode-label {
    font-size: 0.85em;
    font-weight: 600;
    color: var(--fg-mute);
    margin: 0.4em 0 0.6em;
    letter-spacing: 0.04em;
    text-transform: uppercase;
    display: flex;
    align-items: baseline;
    gap: 0.7em;
}
.mode-label .counts {
    font-size: 0.95em;
    color: var(--fg-dim);
    text-transform: none;
    letter-spacing: normal;
}

/* Collapsible mode sub-sections inside a geometry */
details.mode-section {
    margin: 1em 0 0;
    border: 1px solid var(--border-faint);
    border-radius: 10px;
    background: rgba(255, 255, 255, 0.012);
    overflow: hidden;
    scroll-margin-top: 4.5em;
}
details.mode-section > summary {
    list-style: none;
    cursor: pointer;
    padding: 0.7em 1em;
    user-select: none;
    display: flex;
    align-items: baseline;
    gap: 0.8em;
    transition: background 0.15s var(--ease-snap);
}
details.mode-section > summary::-webkit-details-marker { display: none; }
details.mode-section > summary::before {
    content: '▸';
    color: var(--fg-mute);
    transition: transform 0.18s var(--ease-snap);
    font-size: 0.85em;
}
details.mode-section[open] > summary::before { transform: rotate(90deg); }
details.mode-section > summary:hover { background: rgba(255, 255, 255, 0.025); }
details.mode-section .summary-label {
    color: var(--fg);
    font-weight: 600;
    font-size: 0.98em;
}
details.mode-section .counts {
    color: var(--fg-mute);
    font-size: 0.84em;
    font-weight: normal;
}
details.mode-section > .grid { padding: 0.4em 1em 1em; }

/* Deferred section as collapsible details (I2) */
details.deferred-section {
    margin: 2em 0;
    border: 1px solid var(--border);
    border-radius: 12px;
    background: rgba(255, 255, 255, 0.015);
    overflow: hidden;
}
details.deferred-section > summary {
    list-style: none;
    cursor: pointer;
    padding: 0.9em 1.2em;
    user-select: none;
    display: flex;
    align-items: baseline;
    gap: 0.8em;
    transition: background 0.15s ease-out;
}
details.deferred-section > summary::-webkit-details-marker { display: none; }
details.deferred-section > summary::before {
    content: '▸';
    color: var(--fg-mute);
    transition: transform 0.18s ease-out;
    font-size: 0.9em;
}
details.deferred-section[open] > summary::before {
    transform: rotate(90deg);
}
details.deferred-section > summary:hover {
    background: rgba(255, 255, 255, 0.03);
}
details.deferred-section .summary-label {
    color: var(--fg);
    font-weight: 600;
    font-size: 1.05em;
}
details.deferred-section .counts {
    color: var(--fg-mute);
    font-size: 0.85em;
    font-weight: normal;
}
details.deferred-section > .grid {
    padding: 0.5em 1.2em 1.2em;
}
.grid {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(310px, 1fr));
    gap: 1em;
}
.card {
    background: var(--bg-card);
    border: 1px solid rgba(255, 255, 255, 0.04);
    border-radius: 12px;
    overflow: hidden;
    box-shadow:
        0 1px 2px rgba(0, 0, 0, 0.4),
        0 4px 12px rgba(0, 0, 0, 0.25);
    transition: border-color 0.22s var(--ease-snap), transform 0.22s var(--ease-snap), box-shadow 0.22s var(--ease-snap);
    position: relative;
    animation: cardEnter 320ms cubic-bezier(0.2, 0, 0, 1) backwards;
    /* skip render work for offscreen cards (perf) */
    content-visibility: auto;
    contain-intrinsic-size: 0 380px;
}
.card.deferred { background: var(--bg-defer); }
.card:hover {
    border-color: rgba(102, 204, 255, 0.4);
    transform: translateY(-3px);
    box-shadow:
        0 4px 6px rgba(0, 0, 0, 0.45),
        0 12px 24px rgba(0, 0, 0, 0.35);
}
@keyframes cardEnter {
    from { opacity: 0; transform: translateY(8px); }
    to   { opacity: 1; transform: translateY(0); }
}
.grid .card:nth-child(1) { animation-delay:   0ms; }
.grid .card:nth-child(2) { animation-delay:  40ms; }
.grid .card:nth-child(3) { animation-delay:  80ms; }
.grid .card:nth-child(4) { animation-delay: 120ms; }
.grid .card:nth-child(5) { animation-delay: 160ms; }
.grid .card:nth-child(6) { animation-delay: 200ms; }
.grid .card:nth-child(n+7) { animation-delay: 240ms; }
@media (prefers-reduced-motion: reduce) {
    .card { animation: none; }
    .card:hover { transform: none; }
}
.card a {
    display: block;
    text-decoration: none;
    color: inherit;
}
.card .img-wrap {
    position: relative;
    background: #fff;
    height: 200px;
    display: flex;
    align-items: center;
    justify-content: center;
}
.card img {
    max-width: 100%;
    max-height: 200px;
    object-fit: contain;
    display: block;
    outline: 1px solid rgba(255, 255, 255, 0.1);
    outline-offset: -1px;
}
.card .placeholder {
    color: var(--pending);
    font-size: 0.95em;
    font-style: italic;
}
.card .badge {
    position: absolute;
    top: 0.5em;
    right: 0.5em;
    padding: 0.15em 0.55em;
    font-size: 0.75em;
    font-weight: 600;
    border-radius: 4px;
    background: rgba(0,0,0,0.65);
    color: #fff;
}
.card .badge.ok      { background: var(--ok);      color: #06200d; }
.card .badge.stale   { background: var(--stale);   color: #2b1d04; }
.card .badge.pending { background: var(--pending); }
.card .label {
    padding: 0.7em 0.9em 0.55em;
    font-size: 0.86em;
    color: var(--fg);
    display: flex;
    flex-direction: column;
    gap: 0.2em;
}
.card .label .path {
    font-family: ui-monospace, SF Mono, Monaco, Cascadia Code, monospace;
    color: var(--accent);
    word-break: break-word;
    font-size: 0.95em;
}
.card .label .meta-row {
    color: var(--fg-mute);
    font-size: 0.85em;
}
.card .label .expected {
    margin-top: 0.25em;
    color: var(--fg);
    font-size: 0.85em;
    text-wrap: pretty;
    line-height: 1.4;
}

/* ---- new lightbox (Phase 1) ---- */
.lightbox { position: fixed; inset: 0; z-index: 9999; display: flex; align-items: center; justify-content: center; }
.lightbox[hidden] { display: none; }
.lightbox-backdrop { position: absolute; inset: 0; background: rgba(0,0,0,0.92); cursor: pointer; }
.lightbox-img { position: relative; max-width: 92vw; max-height: 84vh; object-fit: contain; box-shadow: 0 12px 36px rgba(0,0,0,0.6); border-radius: 4px; background: #1a1f2a; }
.lightbox-caption { position: absolute; bottom: 4vh; left: 50%; transform: translateX(-50%); color: #d4d8e0; font-family: monospace; font-size: 0.85em; background: rgba(0,0,0,0.5); padding: 0.4em 0.8em; border-radius: 3px; max-width: 80vw; text-align: center; }
.lightbox-nav-hint { position: absolute; top: 2vh; left: 50%; transform: translateX(-50%); color: #8088a0; font-size: 0.72em; white-space: nowrap; }
.lightbox-close { position: absolute; top: 2vh; right: 2vw; background: rgba(255,255,255,0.08); border: 1px solid rgba(255,255,255,0.2); color: #d4d8e0; font-size: 2em; line-height: 1; width: 1.6em; height: 1.6em; border-radius: 50%; cursor: pointer; }
.lightbox-close:hover { background: rgba(255,255,255,0.18); }
button.img-wrap { all: unset; display: block; cursor: pointer; width: 100%; box-sizing: border-box; }
button.img-wrap:focus-visible { outline: 2px solid var(--accent); outline-offset: 2px; }
/* Multi-thumb strip inside a case card */
.thumbs-strip {
    position: relative;
    display: flex;
    gap: 0.5em;
    flex-wrap: wrap;
    align-items: flex-start;
    padding: 0.6em;
    background: var(--bg-subtle);
    border-top: 1px solid var(--border-faint);
    border-bottom: 1px solid var(--border-faint);
}
.thumbs-strip .strip-badge {
    position: absolute;
    top: 0.4em;
    right: 0.4em;
    z-index: 2;
    font-size: 0.65em;
    padding: 0.1em 0.4em;
}
button.thumb {
    all: unset;
    cursor: pointer;
    display: flex;
    flex-direction: column;
    align-items: center;
    gap: 0.2em;
    min-width: 140px;
    max-width: 220px;
    flex: 1 1 160px;
    box-sizing: border-box;
    border: 1px solid var(--border-faint);
    border-radius: 4px;
    overflow: hidden;
    background: var(--bg-card);
    transition: transform 80ms var(--ease-snap), border-color 80ms var(--ease-snap);
}
button.thumb:hover {
    transform: translateY(-1px);
    border-color: var(--accent-2);
}
button.thumb:focus-visible {
    outline: 2px solid var(--accent);
    outline-offset: 2px;
}
.thumb-img {
    width: 100%;
    height: auto;
    max-height: 140px;
    object-fit: contain;
    background: #0e1116;
    display: block;
}
.thumb-placeholder {
    width: 100%;
    height: 90px;
    display: flex;
    align-items: center;
    justify-content: center;
    color: var(--fg-dim);
    font-size: 0.78em;
    background: #0e1116;
}
.thumb-cap {
    font-size: 0.72em;
    color: var(--fg-mute);
    padding: 0.2em 0.4em 0.4em;
    text-align: center;
    line-height: 1.15;
}
.label.label-top { padding: 0.6em 0.8em 0.4em; }

/* click-to-copy */
.copyable {
    cursor: copy;
    transition: background 0.15s ease-out, color 0.15s ease-out, transform 0.1s ease-out;
    border-radius: 3px;
    padding: 1px 3px;
    margin: -1px -3px;
}
.copyable:hover {
    background: rgba(102, 204, 255, 0.08);
}
.copyable:active {
    transform: scale(0.98);
}
.copyable.copied {
    background: var(--ok);
    color: #06200d !important;
}
.copyable.copied::after {
    content: ' ✓ copied';
    font-size: 0.8em;
    margin-left: 0.3em;
}

/* per-case stats row + action button */
.card .stats {
    display: flex;
    gap: 0.8em;
    font-family: ui-monospace, SF Mono, Monaco, monospace;
    font-size: 0.78em;
    color: var(--fg-mute);
    margin-top: 0.4em;
    flex-wrap: wrap;
    font-variant-numeric: tabular-nums;
}
.card .stats .stat { white-space: nowrap; }
.card .stats .stat b { color: var(--fg); font-weight: 600; }
.card .actions {
    border-top: 1px solid var(--border);
    padding: 0.55em 0.9em;
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 0.6em;
    background: rgba(0, 0, 0, 0.15);
    opacity: 0.55;
    transition: opacity 0.18s ease-out;
}
.card:hover .actions,
.card:focus-within .actions,
.card.has-live-job .actions {
    opacity: 1;
}
@media (hover: none) {
    /* touch devices: always show actions, no hover to reveal */
    .card .actions { opacity: 1; }
}
.card .actions .casename {
    font-family: ui-monospace, SF Mono, Monaco, monospace;
    font-size: 0.78em;
    color: var(--fg-mute);
}
.card .actions .casename.missing {
    color: #666;
    font-style: italic;
}
.card .actions .btn-group {
    display: flex;
    gap: 0.35em;
}
.card .actions button {
    background: transparent;
    color: var(--accent);
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 0.5em 0.95em;
    min-height: 32px;
    font-size: 0.82em;
    font-weight: 500;
    cursor: pointer;
    transition: background 0.15s var(--ease-snap), border-color 0.15s var(--ease-snap), transform 0.08s var(--ease-snap);
}
.card .actions button.primary {
    background: rgba(102, 204, 255, 0.12);
    border-color: rgba(102, 204, 255, 0.5);
}
.card .actions button:hover:not(:disabled) {
    border-color: var(--accent);
    background: rgba(102, 204, 255, 0.18);
}
.card .actions button:active:not(:disabled) {
    transform: scale(0.96);
}
.card .actions button:disabled {
    color: var(--fg-mute);
    cursor: not-allowed;
    background: transparent;
    border-color: var(--border);
    opacity: 0.6;
}
.card .job-status {
    font-family: ui-monospace, SF Mono, Monaco, monospace;
    font-size: 0.76em;
    color: var(--fg-mute);
    padding: 0 0.9em 0.6em;
    display: none;
}
.card .job-status.live { display: block; }
.card .job-status .pill {
    display: inline-block;
    padding: 0.1em 0.5em;
    border-radius: 3px;
    background: var(--border);
    color: var(--fg);
    margin-right: 0.4em;
}
.card .job-status .pill.working    { background: #6e7be0; color: #fff; }
.card .job-status .pill.queued     { background: var(--stale); color: #2b1d04; }
.card .job-status .pill.running    { background: #66ccff; color: #003a52; }
.card .job-status .pill.done       { background: var(--ok); color: #06200d; }
.card .job-status .pill.failed     { background: #e07070; color: #2b0808; }
.card .job-status .job-log { color: var(--fg-mute); }

/* live queue panel */
.queue-panel {
    background: var(--bg-card);
    border: 1px solid rgba(255, 255, 255, 0.04);
    border-radius: 12px;
    padding: 0.8em 1em;
    margin: 1em 0 2em;
    font-size: 0.85em;
    box-shadow:
        0 1px 2px rgba(0, 0, 0, 0.4),
        0 4px 12px rgba(0, 0, 0, 0.25);
}
.queue-panel .qhead {
    display: flex;
    justify-content: space-between;
    align-items: baseline;
    margin-bottom: 0.5em;
    gap: 1em;
}
.queue-panel .qhead .qtitle {
    color: var(--fg);
    font-size: 1.05em;
    font-weight: 600;
    letter-spacing: 0.01em;
    display: inline-flex;
    align-items: center;
    gap: 0.55em;
}
.status-dot {
    width: 8px;
    height: 8px;
    border-radius: 50%;
    background: var(--fg-dim);
    flex-shrink: 0;
    transition: background 0.2s var(--ease-snap), box-shadow 0.2s var(--ease-snap);
}
.status-dot[data-state="active"] {
    background: var(--success);
    box-shadow: 0 0 0 0 rgba(110, 224, 138, 0.6);
    animation: pulse-dot 1.8s cubic-bezier(0.4, 0, 0.6, 1) infinite;
}
.status-dot[data-state="idle"]    { background: var(--success); }
.status-dot[data-state="offline"] { background: var(--danger); }
.status-dot[data-state="connecting"] {
    background: var(--info);
    animation: pulse-dot-soft 1.4s ease-in-out infinite;
}
@keyframes pulse-dot {
    0%   { box-shadow: 0 0 0 0 rgba(110, 224, 138, 0.5); }
    100% { box-shadow: 0 0 0 8px rgba(110, 224, 138, 0); }
}
@keyframes pulse-dot-soft {
    50% { opacity: 0.4; }
}
@media (prefers-reduced-motion: reduce) {
    .status-dot { animation: none; }
}
.queue-panel .qhead .tick {
    color: var(--fg-mute);
    font-size: 0.82em;
    font-weight: normal;
    font-variant-numeric: tabular-nums;
}
.queue-panel table {
    width: 100%;
    border-collapse: collapse;
    font-family: ui-monospace, SF Mono, Monaco, monospace;
    font-size: 0.85em;
    font-variant-numeric: tabular-nums;
}
.queue-panel th, .queue-panel td {
    text-align: left;
    padding: 0.25em 0.6em;
    border-bottom: 1px solid var(--border);
}
.queue-panel th { color: var(--fg-mute); font-weight: 600; }
.queue-panel .empty {
    color: var(--fg-mute);
    font-style: italic;
    padding: 0.4em 0;
}
.queue-panel .offline {
    color: var(--danger);
    font-size: 0.85em;
}

/* Keyboard shortcut help overlay */
.kbd-help {
    position: fixed;
    inset: 0;
    background: rgba(0, 0, 0, 0.7);
    display: flex;
    align-items: center;
    justify-content: center;
    z-index: 2000;
    opacity: 0;
    visibility: hidden;
    pointer-events: none;
    transition: opacity 0.15s var(--ease-snap), visibility 0s linear 0.15s;
}
.kbd-help.open {
    opacity: 1;
    visibility: visible;
    pointer-events: auto;
    transition: opacity 0.2s var(--ease-snap), visibility 0s linear 0s;
}
.kbd-help .kbd-panel {
    background: var(--bg-card);
    border: 1px solid var(--border);
    border-radius: 14px;
    padding: 1.5em 1.8em;
    min-width: 340px;
    box-shadow:
        0 4px 12px rgba(0, 0, 0, 0.4),
        0 16px 40px rgba(0, 0, 0, 0.5);
    transform: translateY(8px) scale(0.985);
    transition: transform 0.22s var(--ease-snap);
}
.kbd-help.open .kbd-panel { transform: translateY(0) scale(1); }
.kbd-help h3 {
    margin: 0 0 1em;
    color: var(--fg-strong);
    font-size: 1.05em;
}
.kbd-help dl {
    margin: 0;
    display: grid;
    grid-template-columns: max-content 1fr;
    column-gap: 1.4em;
    row-gap: 0.45em;
    align-items: baseline;
}
.kbd-help dt { color: var(--fg-mute); }
.kbd-help dd { margin: 0; color: var(--fg); }
.kbd-help kbd {
    background: rgba(255, 255, 255, 0.05);
    border: 1px solid var(--border);
    border-bottom-width: 2px;
    border-radius: 4px;
    padding: 0.1em 0.45em;
    font-family: ui-monospace, SF Mono, Monaco, monospace;
    font-size: 0.85em;
    color: var(--fg-strong);
    margin: 0 0.15em;
}
.kbd-help #kbd-help-close {
    margin-top: 1.3em;
    background: transparent;
    color: var(--accent);
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 0.45em 1em;
    font-size: 0.85em;
    cursor: pointer;
    transition: background 0.15s var(--ease-snap), border-color 0.15s var(--ease-snap);
}
.kbd-help #kbd-help-close:hover {
    background: var(--accent-bg);
    border-color: var(--accent);
}

/* filter bar */
.filter-bar { display:flex; gap:0.5rem; flex-wrap:wrap; padding:0.5rem 0; margin:0.5em 0 1em; align-items:center; }
.filter-bar .filter-group { display:flex; gap:0.4rem; flex-wrap:wrap; align-items:center; }
.filter-bar .filter-label { font-size:0.78em; color:var(--fg-mute); padding-right:0.3rem; white-space:nowrap; }
.filter-chip { padding:0.3rem 0.7rem; border:1px solid #888; border-radius:1rem; cursor:pointer; user-select:none; font-size:0.85rem; background:transparent; color:var(--fg-mute); transition:background 0.15s ease-out, color 0.15s ease-out, border-color 0.15s ease-out; }
.filter-chip:hover { border-color:var(--accent); color:var(--accent); }
.filter-chip.active { background:#3366aa; color:white; border-color:#3366aa; }
.filter-chip.clear-btn { border-color:var(--fg-dim); color:var(--fg-dim); }
.filter-chip.clear-btn:hover { border-color:var(--danger); color:var(--danger); }
.card.filtered-out { display:none; }

/* Evidence badges and status pills on cards */
.card .evidence-row {
    display: flex;
    flex-wrap: wrap;
    gap: 0.3em;
    margin-top: 0.4em;
}
.badge-status {
    display: inline-block;
    padding: 0.1em 0.45em;
    border-radius: 999px;
    font-size: 0.72em;
    font-weight: 600;
    letter-spacing: 0.02em;
    text-transform: uppercase;
    white-space: nowrap;
}
.badge-status.ev-image       { background: rgba(110,224,138,0.15); color: #6ee08a; }
.badge-status.ev-residual    { background: rgba(102,204,255,0.15); color: #66ccff; }
.badge-status.ev-spectrum    { background: rgba(110,138,224,0.15); color: #8e9ee0; }
.badge-status.ev-animation   { background: rgba(224,184,110,0.15); color: #e0b86e; }
.badge-status.ev-dns         { background: rgba(224,110,138,0.15); color: #e06e8a; }
.badge-status.ev-checkpoint  { background: rgba(138,110,224,0.15); color: #8a6ee0; }
.badge-status.ev-provenance  { background: rgba(110,200,200,0.15); color: #6ec8c8; }
.badge-status.ev-log         { background: rgba(160,160,160,0.15); color: #a0a0a0; }
.badge-status.ev-other       { background: rgba(100,100,100,0.15); color: #888; }
.badge-status.ev-missing     { background: rgba(224,112,112,0.10); color: #e07070; border: 1px dashed rgba(224,112,112,0.4); }
.badge-stale   { background: rgba(224,184,110,0.18); color: #e0b86e; border: 1px solid rgba(224,184,110,0.4); }
.badge-blocked { background: rgba(224,112,112,0.18); color: #e07070; border: 1px solid rgba(224,112,112,0.4); }
.card .provenance-row {
    font-size: 0.75em;
    color: var(--fg-dim);
    margin-top: 0.25em;
    font-style: italic;
}

/* Didactic ladder — sequential variant groups */
.ladder {
    display: flex;
    flex-direction: column;
    gap: 0;
    margin: 0.6em 0 1.2em;
}
.ladder-step {
    display: flex;
    align-items: flex-start;
    gap: 0.9em;
    position: relative;
}
.ladder-step + .ladder-step::before {
    content: '';
    position: absolute;
    left: 1.05em;
    top: -0.7em;
    width: 2px;
    height: 0.7em;
    background: var(--border);
}
.ladder-step__num {
    width: 2.1em;
    height: 2.1em;
    flex-shrink: 0;
    border-radius: 50%;
    border: 2px solid var(--border);
    display: flex;
    align-items: center;
    justify-content: center;
    font-size: 0.78em;
    font-weight: 700;
    color: var(--fg-mute);
    background: var(--bg-card);
    margin-top: 0.35em;
    font-variant-numeric: tabular-nums;
}
.ladder-step__card {
    flex: 1;
    min-width: 0;
}

/* Responsive grid breakpoints */
@media (max-width: 599px) {
    .grid { grid-template-columns: 1fr; }
    body { padding: 1em 1.2em; }
}
@media (min-width: 600px) and (max-width: 999px) {
    .grid { grid-template-columns: 1fr; }
}
@media (min-width: 1000px) and (max-width: 1399px) {
    .grid { grid-template-columns: repeat(2, 1fr); }
}
@media (min-width: 1400px) {
    .grid { grid-template-columns: repeat(3, 1fr); }
}
/* Text overflow — prevent long labels from breaking layout */
.card .label .path {
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
}
.card .label .meta-row {
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
}
.card .label .expected {
    overflow: hidden;
    display: -webkit-box;
    -webkit-line-clamp: 3;
    -webkit-box-orient: vertical;
}
/* Focus visible for interactive elements */
button:focus-visible,
.filter-chip:focus-visible,
.nav-chip:focus-visible,
a:focus-visible {
    outline: 2px solid var(--accent);
    outline-offset: 2px;
    border-radius: 4px;
}
"""


def _should_ladder(active_cases) -> bool:
    sfd = [mc for mc in active_cases if (mc.sort_prefix or '').startswith('1')]
    newton = [mc for mc in active_cases if (mc.sort_prefix or '').startswith('2')]
    return len(sfd) >= 2 or len(newton) >= 2


def _render_ladder_group(pairs, label, _render_card_fn, step_offset=1):
    bits = [f'<div class="mode-label" style="margin-top:0.8em">{html.escape(label)}</div>']
    bits.append('<div class="ladder">')
    for i, (c, mc) in enumerate(pairs, start=step_offset):
        bits.append(
            f'<div class="ladder-step">'
            f'<div class="ladder-step__num">{i}</div>'
            f'<div class="ladder-step__card">{_render_card_fn(c, mc)}</div>'
            f'</div>'
        )
    bits.append('</div>')
    return ''.join(bits)


def _render_ladder_group_multi(triples, label, _render_card_fn, step_offset=1):
    """Like _render_ladder_group but triples=(c, mc, art_or_None); expands multi-image cases."""
    bits = [f'<div class="mode-label" style="margin-top:0.8em">{html.escape(label)}</div>']
    bits.append('<div class="ladder">')
    step = step_offset
    for i, (c, mc, art) in enumerate(triples):
        if art is not None:
            png = _find_artifact_image(art)
            panel = _panel_caption_from_path(art.path, mc.case_id)
        else:
            png = None
            panel = None
        card_html = _render_card_fn(
            c,
            mc,
            panel_caption=panel,
            png_override=png,
            panel_index=i,
            panel_of=len(triples),
        )
        bits.append(
            f'<div class="ladder-step">'
            f'<div class="ladder-step__num">{step}</div>'
            f'<div class="ladder-step__card">{card_html}</div>'
            f'</div>'
        )
        step += 1
    bits.append('</div>')
    return ''.join(bits)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    now = datetime.now().strftime('%Y-%m-%d %H:%M')
    version = _git_version()

    catalog_cases = list(catalog_all_cases())
    n_validated = sum(1 for c in catalog_cases if c.status == CaseStatus.VALIDATED)
    n_deferred = sum(1 for c in catalog_cases if c.status == CaseStatus.DEFERRED)
    n_blocked = sum(1 for c in catalog_cases if c.status in (
        CaseStatus.BLOCKED, CaseStatus.NEEDS_DNS_SEED, CaseStatus.NEEDS_SCALAR_CHECKPOINT
    ))
    n_partial = sum(1 for c in catalog_cases if c.status in (
        CaseStatus.PARTIAL, CaseStatus.LOCAL_SMOKE, CaseStatus.AMBIGUOUS,
        CaseStatus.STALE, CaseStatus.MISSING_IMAGE, CaseStatus.COMPUTE_HEAVY
    ))
    n_total = len(catalog_cases)

    parts: list[str] = []
    parts.append('<!doctype html>')
    parts.append('<html lang="en"><head>')
    parts.append('<meta charset="utf-8">')
    parts.append('<meta name="viewport" content="width=device-width, initial-scale=1">')
    parts.append('<meta name="color-scheme" content="dark">')
    parts.append('<meta name="robots" content="noindex, nofollow">')
    parts.append(f'<title>nekStab validation gallery — {n_total} cases</title>')
    parts.append(f'<style>{CSS}</style>')
    parts.append('</head><body>')
    parts.append('<header class="page-head">')
    parts.append('<h1>nekStab validation gallery</h1>')
    parts.append(
        f'<div class="meta">'
        f'nekStab <code>{html.escape(version)}</code> &middot; '
        f'{n_validated}/{n_total} validated catalog cases &middot; '
        f'generated {now}'
        f'</div>'
    )
    parts.append('</header>')
    parts.append(
        '<div class="legend">'
        '<span class="ok">match (case Re/Ra)</span>'
        '<span class="stale">stale (Re/Ra mismatch)</span>'
        '<span class="pending">pending re-run</span>'
        '</div>'
    )
    parts.append(
        '<div class="status-summary" style="display:flex;gap:1em;flex-wrap:wrap;margin:0.5em 0 1em;">'
        f'<span class="c-ok">{n_validated} validated</span>'
        f'<span class="c-pending">{n_deferred} deferred</span>'
        f'<span class="c-stale">{n_blocked} blocked / needs seed</span>'
        f'<span class="c-pending">{n_partial} partial</span>'
        f'<span style="color:var(--fg-mute);font-size:0.85em;">&middot; {n_total} total &middot; {now}</span>'
        '</div>'
    )

    parts.append(
        '<div class="filter-bar" id="filter-bar" aria-label="Filter cards">'
        '<div class="filter-group">'
        '<span class="filter-label">Status:</span>'
        '<button class="filter-chip" aria-label="Filter by validated status" data-filter-group="status" data-filter-value="validated">validated</button>'
        '<button class="filter-chip" aria-label="Filter by deferred status" data-filter-group="status" data-filter-value="deferred">deferred</button>'
        '<button class="filter-chip" aria-label="Filter by blocked status" data-filter-group="status" data-filter-value="blocked">blocked</button>'
        '</div>'
        '<div class="filter-group">'
        '<span class="filter-label">Lane:</span>'
        '<button class="filter-chip" aria-label="Filter by baseflow lane" data-filter-group="lane" data-filter-value="baseflow">baseflow</button>'
        '<button class="filter-chip" aria-label="Filter by stability lane" data-filter-group="lane" data-filter-value="stability">stability</button>'
        '<button class="filter-chip" aria-label="Filter by postproc lane" data-filter-group="lane" data-filter-value="postproc">postproc</button>'
        '<button class="filter-chip" aria-label="Filter by dns lane" data-filter-group="lane" data-filter-value="dns">dns</button>'
        '<button class="filter-chip" aria-label="Filter by modal lane" data-filter-group="lane" data-filter-value="modal">modal</button>'
        '<button class="filter-chip" aria-label="Filter by otd lane" data-filter-group="lane" data-filter-value="otd">otd</button>'
        '</div>'
        '<button class="filter-chip clear-btn" aria-label="Clear all filters" id="filter-clear">Clear filters</button>'
        '</div>'
    )

    parts.append(
        '<div class="queue-panel" id="queue-panel">'
        '<div class="qhead">'
        '<span class="qtitle">'
        '<span class="status-dot" id="queue-status-dot" data-state="connecting" aria-hidden="true"></span>'
        'Running jobs'
        '</span>'
        '<span class="tick" id="queue-tick">connecting&hellip;</span>'
        '</div>'
        '<div id="queue-body"><div class="empty">waiting for first poll&hellip;</div></div>'
        '</div>'
    )

    def _render_card(c: Case, mc=None, panel_caption: str | None = None, png_override=_PNG_NOT_SET, panel_index: int = 0, panel_of: int = 1) -> str:
        if png_override is _PNG_NOT_SET:
            png = _find_png(c)
        else:
            png = png_override
        if c.status == 'validated' and png is not None and not _is_stale(c, png):
            state = 'ok'
        elif c.status == 'blocked':
            state = 'stale'
        elif png is None:
            state = 'pending'
        elif _is_stale(c, png):
            state = 'stale'
        else:
            state = 'pending'
        card_classes = 'card'
        if c.status != 'validated':
            card_classes += ' deferred'
        # Derive data attributes for filter chips
        family_val = html.escape(mc.flow_family if mc is not None else c.path.split('/')[0])
        lane_raw = mc.method_lane if mc is not None else _classify(c)[1]
        # Normalize method_lane values to display groups
        _lane_map = {
            'baseflow': 'baseflow', 'dns': 'dns', 'modal': 'modal',
            'otd': 'otd', 'rans': 'rans', 'wavemaker': 'postproc',
            'stability': 'stability', 'stability_direct': 'stability',
            'stability_adjoint': 'stability',
        }
        if lane_raw in _lane_map:
            lane_val = _lane_map[lane_raw]
        elif 'floquet' in lane_raw or 'direct' in lane_raw or 'adjoint' in lane_raw or 'transient' in lane_raw or 'animation' in lane_raw:
            lane_val = 'stability'
        elif 'wavemaker' in lane_raw:
            lane_val = 'postproc'
        elif 'post' in lane_raw:
            lane_val = 'postproc'
        else:
            lane_val = lane_raw
        lane_val = html.escape(lane_val)
        re_val = html.escape(c.re_tag)
        case_id_val = html.escape(mc.case_id if mc is not None else c.path)
        cap_text = f"{mc.case_id if mc is not None else c.path} · {panel_caption}" if panel_caption else (mc.case_id if mc is not None else c.path)
        img_src_val = f"figures/{html.escape(png.name)}" if png is not None else ""
        chunks: list[str] = [
            f'<div class="{card_classes}"'
            f' data-case-id="{case_id_val}"'
            f' data-panel-index="{panel_index}"'
            f' data-panel-of="{panel_of}"'
            f' data-img-src="{html.escape(img_src_val)}"'
            f' data-caption="{html.escape(cap_text)}"'
            f' data-family="{family_val}"'
            f' data-lane="{lane_val}"'
            f' data-status="{html.escape(c.status)}"'
            f' data-re="{re_val}">'
        ]
        caption = mc.current_path if mc is not None else c.path
        if png is not None:
            chunks.append(
                f'<button class="img-wrap" aria-label="Open figure" '
                f'data-img-src="figures/{html.escape(png.name)}" '
                f'data-caption="{html.escape(cap_text)}">'
            )
        else:
            chunks.append(
                f'<button class="img-wrap" aria-label="Open figure" '
                f'data-img-src="" data-caption="{html.escape(cap_text)}">'
            )
        badge_label = {'ok': 'OK', 'stale': 'STALE', 'pending': 'PENDING'}[state]
        chunks.append(f'<span class="badge {state}">{badge_label}</span>')
        if png is not None:
            chunks.append(
                f'<img loading="lazy" decoding="async" '
                f'src="figures/{html.escape(png.name)}" alt="{html.escape(caption)}">'
            )
        else:
            msg = 'pending re-run' if c.status == 'validated' else c.status
            chunks.append(f'<div class="placeholder">{msg}</div>')
        chunks.append('</button><div class="label">')
        chunks.append(
            f'<div class="copyable path" data-action="copy" data-copy="{html.escape(c.path)}">'
            f'{html.escape(caption)}</div>'
        )
        label = mc.label if mc is not None else c.mode
        mode_name = mc.mode_name if mc is not None else c.mode
        legacy_uparam01 = mc.legacy_uparam01 if mc is not None else None
        meta_bits = [c.re_tag, label]
        if mode_name:
            meta_bits.append(mode_name)
        if legacy_uparam01:
            meta_bits.append(f'uparam01 {legacy_uparam01}')
        meta_text = ' · '.join(meta_bits)
        chunks.append(
            f'<div class="copyable meta-row" data-action="copy" data-copy="{html.escape(meta_text)}">'
            f'{html.escape(c.re_tag)} &middot; {html.escape(label)}'
            f'{(" &middot; " + html.escape(mode_name)) if mode_name else ""}'
            f'{(" &middot; uparam01 " + html.escape(legacy_uparam01)) if legacy_uparam01 else ""}</div>'
        )
        if panel_caption:
            chunks.append(
                f'<div class="meta-row" style="color:var(--accent);font-size:0.82em;">'
                f'{html.escape(panel_caption)}</div>'
            )
        chunks.append(
            f'<div class="copyable expected" data-action="copy" data-copy="{html.escape(c.expected)}">'
            f'{html.escape(c.expected)}</div>'
        )
        if mc is not None and mc.why_note:
            chunks.append(
                f'<div class="expected" style="color:var(--fg-dim);font-size:0.8em;'
                f'font-style:italic;margin-top:0.15em">{html.escape(mc.why_note[:120])}</div>'
            )
        if mc is not None and mc.evidence_requirements:
            badge_html_parts = []
            for req in mc.evidence_requirements:
                role_val = req.role.value
                css_mod = _evidence_role_css(role_val)
                label = _evidence_role_label(role_val)
                req_cls = 'badge-status ' + css_mod
                badge_html_parts.append(
                    f'<span class="{req_cls}" title="{html.escape(role_val)}">{label}</span>'
                )
            if badge_html_parts:
                chunks.append(
                    '<div class="evidence-row">' + ''.join(badge_html_parts) + '</div>'
                )
        if mc is not None and any(getattr(a, 'freshness', None) == 'stale' for a in mc.artifacts):
            chunks.append(
                '<div class="evidence-row">'
                '<span class="badge-status badge-stale" title="One or more artifacts are stale">STALE</span>'
                '</div>'
            )
        if mc is not None and mc.blockers:
            tooltip = html.escape(' | '.join(str(b)[:80] for b in mc.blockers))
            chunks.append(
                f'<div class="evidence-row">'
                f'<span class="badge-status badge-blocked" title="{tooltip}">BLOCKED</span>'
                f'</div>'
            )
        prov_parts = []
        if mc is not None and mc.seed_parameters and mc.seed_parameters.seed_status:
            prov_parts.append(f'seed: {html.escape(mc.seed_parameters.seed_status)}')
        if mc is not None:
            for a in mc.artifacts:
                if a.role.value == 'provenance' and a.provenance:
                    prov_parts.append(html.escape(a.provenance[:60]))
                    break
        if prov_parts:
            chunks.append(
                f'<div class="provenance-row">{" &middot; ".join(prov_parts)}</div>'
            )
        stats = _case_stats(c.path)
        stat_bits: list[str] = []
        if 'nelg' in stats:
            stat_bits.append(f'<span class="stat"><b>{stats["nelg"]}</b> el</span>')
        if 'ntasks' in stats:
            stat_bits.append(f'<span class="stat"><b>{stats["ntasks"]}</b> mpi</span>')
        if 'elapsed_sec' in stats:
            stat_bits.append(
                f'<span class="stat"><b>'
                f'{html.escape(_format_elapsed(stats["elapsed_sec"]))}</b></span>'
            )
        if stat_bits:
            chunks.append('<div class="stats">' + ''.join(stat_bits) + '</div>')
        chunks.append('</div>')
        casename = stats.get('casename', '')
        has_slurm = stats.get('has_slurm', False)
        chunks.append('<div class="actions">')
        if casename:
            chunks.append(f'<span class="casename">{html.escape(casename)}</span>')
        else:
            chunks.append('<span class="casename missing">no session</span>')
        recompile_attrs = '' if casename else (
            ' disabled data-static-disabled="1"'
            ' title="case has no SESSION.NAME — cannot compile"'
        )
        resubmit_attrs = '' if has_slurm else (
            ' disabled data-static-disabled="1"'
            ' title="case has no Slurm script — cannot submit"'
        )
        elapsed_attr = ''
        if 'elapsed_sec' in stats:
            elapsed_attr = f' data-elapsed="{int(stats["elapsed_sec"])}"'
        chunks.append('<div class="btn-group">')
        chunks.append(
            f'<button class="action-btn" data-action="recompile" '
            f'aria-label="Recompile {html.escape(c.path)}" '
            f'data-case="{html.escape(c.path)}"{elapsed_attr}{recompile_attrs}>'
            f'recompile</button>'
        )
        chunks.append(
            f'<button class="action-btn primary" data-action="resubmit" '
            f'aria-label="Resubmit {html.escape(c.path)}" '
            f'data-case="{html.escape(c.path)}"{elapsed_attr}{resubmit_attrs}>'
            f'resubmit</button>'
        )
        chunks.append('</div></div>')
        chunks.append(
            f'<div class="job-status" data-case="{html.escape(c.path)}"'
            f' aria-live="polite" aria-atomic="true"></div></div>'
        )
        return ''.join(chunks)

    def _render_case_card(c: Case, mc, panels: list) -> str:
        """One card per simulation with N thumbnails inside.

        panels: list of (artifact, png_path, panel_caption) tuples.
        Empty list -> single pending card (delegates to _render_card).
        """
        if not panels:
            return _render_card(c, mc, panel_caption=None, png_override=None,
                                panel_index=0, panel_of=1)
        first_png = next((p[1] for p in panels if p[1] is not None), None)
        if c.status == 'validated' and first_png is not None and not _is_stale(c, first_png):
            state = 'ok'
        elif c.status == 'blocked':
            state = 'stale'
        elif first_png is None:
            state = 'pending'
        elif _is_stale(c, first_png):
            state = 'stale'
        else:
            state = 'pending'
        card_classes = 'card'
        if c.status != 'validated':
            card_classes += ' deferred'
        family_val = html.escape(mc.flow_family if mc is not None else c.path.split('/')[0])
        lane_raw = mc.method_lane if mc is not None else _classify(c)[1]
        _lane_map_local = {
            'baseflow': 'baseflow', 'dns': 'dns', 'modal': 'modal',
            'otd': 'otd', 'rans': 'rans', 'wavemaker': 'postproc',
            'stability': 'stability', 'stability_direct': 'stability',
            'stability_adjoint': 'stability',
        }
        if lane_raw in _lane_map_local:
            lane_val = _lane_map_local[lane_raw]
        elif ('floquet' in lane_raw or 'direct' in lane_raw or 'adjoint' in lane_raw
              or 'transient' in lane_raw or 'animation' in lane_raw):
            lane_val = 'stability'
        elif 'wavemaker' in lane_raw or 'post' in lane_raw:
            lane_val = 'postproc'
        else:
            lane_val = lane_raw
        lane_val = html.escape(lane_val)
        re_val = html.escape(c.re_tag)
        case_id_val = html.escape(mc.case_id if mc is not None else c.path)
        caption = mc.current_path if mc is not None else c.path

        chunks: list[str] = [
            f'<div class="{card_classes}"'
            f' data-case-id="{case_id_val}"'
            f' data-family="{family_val}"'
            f' data-lane="{lane_val}"'
            f' data-status="{html.escape(c.status)}"'
            f' data-re="{re_val}">'
        ]
        chunks.append('<div class="label label-top">')
        chunks.append(
            f'<div class="copyable path" data-action="copy" data-copy="{html.escape(c.path)}">'
            f'{html.escape(caption)}</div>'
        )
        label = mc.label if mc is not None else c.mode
        mode_name = mc.mode_name if mc is not None else c.mode
        legacy_uparam01 = mc.legacy_uparam01 if mc is not None else None
        meta_bits = [c.re_tag, label]
        if mode_name:
            meta_bits.append(mode_name)
        if legacy_uparam01:
            meta_bits.append(f'uparam01 {legacy_uparam01}')
        meta_text = ' · '.join(meta_bits)
        chunks.append(
            f'<div class="copyable meta-row" data-action="copy" data-copy="{html.escape(meta_text)}">'
            f'{html.escape(c.re_tag)} &middot; {html.escape(label)}'
            f'{(" &middot; " + html.escape(mode_name)) if mode_name else ""}'
            f'{(" &middot; uparam01 " + html.escape(legacy_uparam01)) if legacy_uparam01 else ""}</div>'
        )
        chunks.append(
            f'<div class="copyable expected" data-action="copy" data-copy="{html.escape(c.expected)}">'
            f'{html.escape(c.expected)}</div>'
        )
        if mc is not None and mc.why_note:
            chunks.append(
                f'<div class="expected" style="color:var(--fg-dim);font-size:0.8em;'
                f'font-style:italic;margin-top:0.15em">{html.escape(mc.why_note[:120])}</div>'
            )
        if mc is not None and mc.evidence_requirements:
            badge_html_parts = []
            for req in mc.evidence_requirements:
                role_val = req.role.value
                css_mod = _evidence_role_css(role_val)
                lbl = _evidence_role_label(role_val)
                req_cls = 'badge-status ' + css_mod
                badge_html_parts.append(
                    f'<span class="{req_cls}" title="{html.escape(role_val)}">{lbl}</span>'
                )
            if badge_html_parts:
                chunks.append('<div class="evidence-row">' + ''.join(badge_html_parts) + '</div>')
        if mc is not None and mc.blockers:
            tooltip = html.escape(' | '.join(str(b)[:80] for b in mc.blockers))
            chunks.append(
                f'<div class="evidence-row">'
                f'<span class="badge-status badge-blocked" title="{tooltip}">BLOCKED</span>'
                f'</div>'
            )
        stats_local = _case_stats(c.path)
        stat_bits: list[str] = []
        if 'nelg' in stats_local:
            stat_bits.append(f'<span class="stat"><b>{stats_local["nelg"]}</b> el</span>')
        if 'ntasks' in stats_local:
            stat_bits.append(f'<span class="stat"><b>{stats_local["ntasks"]}</b> mpi</span>')
        if 'elapsed_sec' in stats_local:
            stat_bits.append(
                f'<span class="stat"><b>'
                f'{html.escape(_format_elapsed(stats_local["elapsed_sec"]))}</b></span>'
            )
        if stat_bits:
            chunks.append('<div class="stats">' + ''.join(stat_bits) + '</div>')
        chunks.append('</div>')

        # ----- Multi-panel thumbnail strip -----
        badge_label = {'ok': 'OK', 'stale': 'STALE', 'pending': 'PENDING'}[state]
        chunks.append(f'<div class="thumbs-strip"><span class="badge {state} strip-badge">{badge_label}</span>')
        for i, (art, png, panel_caption) in enumerate(panels):
            thumb_src = f'figures/{html.escape(png.name)}' if png is not None else ''
            base_label = mc.case_id if mc is not None else c.path
            thumb_cap = f'{base_label} · {panel_caption}' if panel_caption else base_label
            chunks.append(
                f'<button class="thumb" '
                f'data-case-id="{case_id_val}" '
                f'data-panel-index="{i}" '
                f'data-panel-of="{len(panels)}" '
                f'data-img-src="{html.escape(thumb_src)}" '
                f'data-caption="{html.escape(thumb_cap)}" '
                f'aria-label="Open panel {i+1} of {len(panels)}">'
            )
            if png is not None:
                chunks.append(
                    f'<img class="thumb-img" loading="lazy" decoding="async" '
                    f'src="{thumb_src}" alt="{html.escape(thumb_cap)}">'
                )
            else:
                chunks.append('<div class="thumb-placeholder">pending</div>')
            if panel_caption:
                chunks.append(
                    f'<div class="thumb-cap">{html.escape(panel_caption)}</div>'
                )
            chunks.append('</button>')
        chunks.append('</div>')

        # ----- ACTIONS -----
        casename = stats_local.get('casename', '')
        has_slurm = stats_local.get('has_slurm', False)
        chunks.append('<div class="actions">')
        if casename:
            chunks.append(f'<span class="casename">{html.escape(casename)}</span>')
        else:
            chunks.append('<span class="casename missing">no session</span>')
        recompile_attrs = '' if casename else (
            ' disabled data-static-disabled="1"'
            ' title="case has no SESSION.NAME — cannot compile"'
        )
        resubmit_attrs = '' if has_slurm else (
            ' disabled data-static-disabled="1"'
            ' title="case has no Slurm script — cannot submit"'
        )
        elapsed_attr = ''
        if 'elapsed_sec' in stats_local:
            elapsed_attr = f' data-elapsed="{int(stats_local["elapsed_sec"])}"'
        chunks.append('<div class="btn-group">')
        chunks.append(
            f'<button class="action-btn" data-action="recompile" '
            f'aria-label="Recompile {html.escape(c.path)}" '
            f'data-case="{html.escape(c.path)}"{elapsed_attr}{recompile_attrs}>'
            f'recompile</button>'
        )
        chunks.append(
            f'<button class="action-btn primary" data-action="resubmit" '
            f'aria-label="Resubmit {html.escape(c.path)}" '
            f'data-case="{html.escape(c.path)}"{elapsed_attr}{resubmit_attrs}>'
            f'resubmit</button>'
        )
        chunks.append('</div></div>')
        chunks.append(
            f'<div class="job-status" data-case="{html.escape(c.path)}"'
            f' aria-live="polite" aria-atomic="true"></div></div>'
        )
        return ''.join(chunks)

    def _mc_panels(mc) -> list:
        """Build [(art, png, panel_caption)] for all reference_image artifacts of a MethodCase."""
        result = []
        for art in mc.artifacts:
            if art.role.value != 'reference_image':
                continue
            png = _find_artifact_image(art)
            panel = _panel_caption_from_path(art.path, mc.case_id)
            result.append((art, png, panel))
        return result

    all_families = {f.family_id: f for f in catalog_families()}
    ordered_fids = [fid for fid in FAMILY_ORDER if fid in all_families] + [
        fid for fid in all_families if fid not in FAMILY_ORDER
    ]
    nav_chips: list[str] = []
    for fid in ordered_fids:
        fam = all_families[fid]
        fam_cases = cases_for(fid)
        n_total_fam = len(fam_cases)
        n_val_fam = sum(1 for c in fam_cases if c.status == CaseStatus.VALIDATED)
        nav_chips.append(
            f'<a class="nav-chip" href="#family-{fid}" aria-label="Jump to {html.escape(fam.label)} family">'
            f'{html.escape(fam.label)}'
            f'<span class="count">{n_val_fam}/{n_total_fam}</span>'
            f'</a>'
        )
    parts.append(
        '<nav class="section-nav" id="family-nav" aria-label="Flow family navigation">'
        + ''.join(nav_chips)
        + '</nav>'
    )

    def _family_counter_chips(fam_cases) -> str:
        n_val = sum(1 for c in fam_cases if c.status == CaseStatus.VALIDATED)
        n_def = sum(1 for c in fam_cases if c.status == CaseStatus.DEFERRED)
        n_blk = sum(1 for c in fam_cases if c.status in (
            CaseStatus.BLOCKED, CaseStatus.NEEDS_DNS_SEED, CaseStatus.NEEDS_SCALAR_CHECKPOINT
        ))
        n_prt = sum(1 for c in fam_cases if c.status in (
            CaseStatus.PARTIAL, CaseStatus.LOCAL_SMOKE, CaseStatus.AMBIGUOUS,
            CaseStatus.STALE, CaseStatus.MISSING_IMAGE, CaseStatus.COMPUTE_HEAVY
        ))
        bits = []
        if n_val:
            bits.append(f'<span class="c-ok">{n_val} validated</span>')
        if n_prt:
            bits.append(f'<span class="c-pending">{n_prt} partial</span>')
        if n_blk:
            bits.append(f'<span class="c-stale">{n_blk} blocked</span>')
        if n_def:
            bits.append(f'<span class="c-pending">{n_def} deferred</span>')
        return ''.join(bits)

    for fid in ordered_fids:
        fam = all_families[fid]
        fam_cases = list(cases_for(fid))
        active_cases = [c for c in fam_cases if c.status != CaseStatus.DEFERRED]
        deferred_cases = [c for c in fam_cases if c.status == CaseStatus.DEFERRED]
        active_cases = [mc for mc in active_cases if getattr(mc, 'gallery_visible', True)]
        deferred_cases = [mc for mc in deferred_cases if getattr(mc, 'gallery_visible', True)]
        active_cases.sort(key=lambda mc: LANE_ORDER.index(mc.method_lane) if mc.method_lane in LANE_ORDER else 999)
        parts.append(
            f'<section class="geom" id="family-{fid}" style="scroll-margin-top:4.5em">'
            f'<h2>{html.escape(fam.label)} '
            f'<span style="font-size:0.7em;color:var(--fg-dim);font-weight:normal">'
            f'{html.escape(fam.geometry_kind)}</span> '
            f'<span class="counts">{_family_counter_chips(fam_cases)}</span></h2>'
        )
        if fam.why_note:
            parts.append(
                f'<p style="color:var(--fg-mute);font-size:0.85em;font-style:italic;'
                f'margin:0 0 1em">{html.escape(fam.why_note)}</p>'
            )
        if active_cases:
            # DNS cases (lane='dns') render first, before SFD/Newton ladders.
            # Per gallery design: DNS at top of each family.
            dns_mcs = sorted(
                [mc for mc in active_cases if mc.method_lane == 'dns'],
                key=lambda m: (m.sort_prefix or '')
            )
            non_dns_active = [mc for mc in active_cases if mc.method_lane != 'dns']
            if dns_mcs:
                parts.append('<div class="grid">')
                for mc in dns_mcs:
                    c = catalog_to_case(mc)
                    parts.append(_render_case_card(c, mc, _mc_panels(mc)))
                parts.append('</div>')
            if _should_ladder(non_dns_active):
                sfd_mcs = sorted(
                    [mc for mc in non_dns_active if (mc.sort_prefix or '').startswith('1')],
                    key=lambda m: (m.sort_prefix or '')
                )
                newton_mcs = sorted(
                    [mc for mc in non_dns_active if (mc.sort_prefix or '').startswith('2')],
                    key=lambda m: (m.sort_prefix or '')
                )
                other_mcs = [
                    mc for mc in non_dns_active
                    if not (mc.sort_prefix or '').startswith('1')
                    and not (mc.sort_prefix or '').startswith('2')
                ]
                if len(sfd_mcs) >= 2:
                    parts.append(
                        '<div class="mode-label" style="margin-top:0.8em">'
                        'SFD variants — baseline → dynamic tolerance → OIFS'
                        '</div><div class="ladder">'
                    )
                    for step_i, mc in enumerate(sfd_mcs, 1):
                        c = catalog_to_case(mc)
                        parts.append(
                            f'<div class="ladder-step">'
                            f'<div class="ladder-step__num">{step_i}</div>'
                            f'<div class="ladder-step__card">{_render_case_card(c, mc, _mc_panels(mc))}</div>'
                            f'</div>'
                        )
                    parts.append('</div>')
                elif sfd_mcs:
                    other_mcs = sfd_mcs + other_mcs
                if len(newton_mcs) >= 2:
                    parts.append(
                        '<div class="mode-label" style="margin-top:0.8em">'
                        'Newton variants — standard → UPO → forced'
                        '</div><div class="ladder">'
                    )
                    for step_i, mc in enumerate(newton_mcs, 1):
                        c = catalog_to_case(mc)
                        parts.append(
                            f'<div class="ladder-step">'
                            f'<div class="ladder-step__num">{step_i}</div>'
                            f'<div class="ladder-step__card">{_render_case_card(c, mc, _mc_panels(mc))}</div>'
                            f'</div>'
                        )
                    parts.append('</div>')
                elif newton_mcs:
                    other_mcs = newton_mcs + other_mcs
                if other_mcs:
                    parts.append('<div class="grid">')
                    for mc in other_mcs:
                        c = catalog_to_case(mc)
                        parts.append(_render_case_card(c, mc, _mc_panels(mc)))
                    parts.append('</div>')
            elif non_dns_active:
                parts.append('<div class="grid">')
                for mc in non_dns_active:
                    c = catalog_to_case(mc)
                    parts.append(_render_case_card(c, mc, _mc_panels(mc)))
                parts.append('</div>')
        if deferred_cases:
            n_def = len(deferred_cases)
            parts.append(
                f'<details class="deferred-section">'
                f'<summary><span class="summary-label">Deferred</span> '
                f'<span class="counts">{n_def} case{"s" if n_def != 1 else ""}</span></summary>'
                f'<div class="grid">'
            )
            for mc in deferred_cases:
                img_arts = [a for a in mc.artifacts if a.role.value == 'reference_image']
                if img_arts:
                    for i, art in enumerate(img_arts):
                        c = catalog_to_case(mc)
                        png = _find_artifact_image(art)
                        panel = _panel_caption_from_path(art.path, mc.case_id)
                        parts.append(_render_card(c, mc, panel_caption=panel, png_override=png, panel_index=i, panel_of=len(img_arts)))
                else:
                    c = catalog_to_case(mc)
                    parts.append(_render_card(c, mc, panel_index=0, panel_of=1))
            parts.append('</div></details>')
        parts.append('</section>')

    parts.append(
        '<div class="lightbox" id="lightbox" hidden role="dialog" aria-modal="true" aria-label="Figure viewer">'
        '<div class="lightbox-backdrop"></div>'
        '<button class="lightbox-close" aria-label="Close">&times;</button>'
        '<img class="lightbox-img" alt="">'
        '<div class="lightbox-caption"></div>'
        '<div class="lightbox-nav-hint">&#8592; prev panel &middot; &#8594; next panel &middot; &#8593; prev case &middot; &#8595; next case &middot; esc close</div>'
        '</div>'
        '<div class="kbd-help" id="kbd-help" role="dialog" aria-modal="true"'
        ' aria-labelledby="kbd-help-title" aria-hidden="true">'
        '<div class="kbd-panel">'
        '<h3 id="kbd-help-title">Keyboard shortcuts</h3>'
        '<dl>'
        '<dt><kbd>?</kbd></dt><dd>Toggle this help</dd>'
        '<dt><kbd>Esc</kbd></dt><dd>Close lightbox / help</dd>'
        '<dt><kbd>g</kbd> then <kbd>1</kbd>&hellip;<kbd>9</kbd></dt>'
        '<dd>Jump to N-th geometry</dd>'
        '<dt><kbd>t</kbd></dt><dd>Scroll to top</dd>'
        '<dt><kbd>r</kbd></dt><dd>Refresh queue panel</dd>'
        '</dl>'
        '<button id="kbd-help-close" type="button">close</button>'
        '</div>'
        '</div>'
    )
    parts.append(r"""<script>
(function () {
    // ---- new lightbox (Phase 1) ----
    (function() {
        const lb = document.getElementById('lightbox');
        if (!lb) return;
        const lbImg = lb.querySelector('.lightbox-img');
        const lbCap = lb.querySelector('.lightbox-caption');
        const lbBackdrop = lb.querySelector('.lightbox-backdrop');
        const lbClose = lb.querySelector('.lightbox-close');
        function rebuild() {
            // Cases in DOM order (one .card per case in the new layout)
            const cards = Array.from(document.querySelectorAll('.card[data-case-id]'));
            const caseOrder = [];
            const byCase = new Map();
            for (const c of cards) {
                const cid = c.dataset.caseId;
                if (byCase.has(cid)) continue;  // skip dupes; should be 1:1 now
                caseOrder.push(cid);
                // Each case's panels are the thumb buttons inside its card.
                const thumbs = Array.from(c.querySelectorAll('button.thumb[data-case-id]'));
                thumbs.sort((a, b) => parseInt(a.dataset.panelIndex||'0') - parseInt(b.dataset.panelIndex||'0'));
                // Fallback: a card with no thumbs (pending) — represent itself as a single "panel"
                byCase.set(cid, thumbs.length ? thumbs : [c]);
            }
            return { caseOrder, byCase };
        }
        let state = { caseOrder: [], byCase: new Map(), currentCase: null, currentPanel: 0 };
        function show(caseId, panelIdx) {
            const panels = state.byCase.get(caseId) || [];
            if (!panels.length) return;
            panelIdx = Math.max(0, Math.min(panelIdx, panels.length - 1));
            const card = panels[panelIdx];
            const src = card.dataset.imgSrc;
            if (!src) {
                lbImg.removeAttribute('src');
                lbCap.textContent = (card.dataset.caption || caseId) + ' (pending)';
            } else {
                lbImg.src = src;
                lbImg.alt = card.dataset.caption || '';
                lbCap.textContent = card.dataset.caption || caseId;
            }
            state.currentCase = caseId;
            state.currentPanel = panelIdx;
            lb.hidden = false;
            document.body.style.overflow = 'hidden';
        }
        function closeLb() {
            lb.hidden = true;
            document.body.style.overflow = '';
        }
        function navigate(dCase, dPanel) {
            if (dCase) {
                const idx = state.caseOrder.indexOf(state.currentCase);
                const newIdx = (idx + dCase + state.caseOrder.length) % state.caseOrder.length;
                show(state.caseOrder[newIdx], 0);
            } else if (dPanel) {
                show(state.currentCase, state.currentPanel + dPanel);
            }
        }
        document.addEventListener('click', function(e) {
            // Clicking a multi-thumb opens lightbox at that panel.
            const thumbBtn = e.target.closest('button.thumb');
            if (thumbBtn) {
                e.preventDefault();
                const fresh = rebuild();
                state.caseOrder = fresh.caseOrder;
                state.byCase = fresh.byCase;
                show(thumbBtn.dataset.caseId, parseInt(thumbBtn.dataset.panelIndex||'0'));
                return;
            }
            // Legacy single-panel button (img-wrap inside a pending card).
            const imgBtn = e.target.closest('button.img-wrap');
            if (imgBtn) {
                e.preventDefault();
                const fresh = rebuild();
                state.caseOrder = fresh.caseOrder;
                state.byCase = fresh.byCase;
                const card = imgBtn.closest('.card');
                show(card.dataset.caseId, 0);
            }
        });
        lbBackdrop.addEventListener('click', closeLb);
        lbClose.addEventListener('click', closeLb);
        document.addEventListener('keydown', function(e) {
            if (lb.hidden) return;
            if (e.key === 'Escape')     { e.preventDefault(); closeLb(); }
            else if (e.key === 'ArrowLeft')  { e.preventDefault(); navigate(0, -1); }
            else if (e.key === 'ArrowRight') { e.preventDefault(); navigate(0, +1); }
            else if (e.key === 'ArrowUp')    { e.preventDefault(); navigate(-1, 0); }
            else if (e.key === 'ArrowDown')  { e.preventDefault(); navigate(+1, 0); }
        });
    })();
    // ---- end new lightbox ----

    // ----- click-to-copy on .copyable spans -----
    function flashCopied(el) {
        el.classList.add('copied');
        clearTimeout(el._flashTimer);
        el._flashTimer = setTimeout(function () { el.classList.remove('copied'); }, 900);
    }
    document.body.addEventListener('click', function (e) {
        const el = e.target.closest('.copyable');
        if (!el) return;
        // don't hijack clicks inside image buttons
        if (e.target.closest('button.img-wrap')) return;
        e.preventDefault();
        e.stopPropagation();
        const text = el.dataset.copy || el.textContent;
        if (navigator.clipboard) {
            navigator.clipboard.writeText(text).then(function () { flashCopied(el); });
        } else {
            // fallback for non-HTTPS contexts
            const ta = document.createElement('textarea');
            ta.value = text;
            document.body.appendChild(ta);
            ta.select();
            try { document.execCommand('copy'); flashCopied(el); } catch (_) {}
            document.body.removeChild(ta);
        }
    });

    // ----- resubmit + per-job polling -----
    const jobsByCase = new Map(); // case path -> {token, pollId, sibs}

    // I4: collapse server vocabulary -> user vocabulary
    // server: starting|compiling|submitting|queued|running|done|failed|error
    // user:   working | queued | running | done | failed
    function uiState(serverState) {
        const s = (serverState || '').toLowerCase();
        if (s === 'starting' || s === 'compiling' || s === 'submitting') return 'working';
        if (s === 'queued')   return 'queued';
        if (s === 'running')  return 'running';
        if (s === 'done')     return 'done';
        if (s === 'failed' || s === 'error') return 'failed';
        return s || 'unknown';
    }
    function uiLabel(uiS) {
        return ({
            working: 'Working',
            queued:  'Queued',
            running: 'Running',
            done:    'Done',
            failed:  'Failed',
        })[uiS] || uiS;
    }
    // C2: human-friendly error mapping
    function friendlyError(raw) {
        const msg = String(raw || '');
        if (/Unexpected token|JSON\.parse|is not valid JSON/i.test(msg))
            return "Can't reach server — is validation/serve.py running and tunnel up?";
        if (/Failed to fetch|NetworkError|TypeError/i.test(msg))
            return "Server unreachable — check the SSH tunnel.";
        if (/403/.test(msg) && /catalog/i.test(msg))
            return "Case isn't in the catalog — can't resubmit.";
        if (/403/.test(msg) && /origin/i.test(msg))
            return "Origin blocked — open the gallery via http://localhost:8000/, not file://.";
        if (/409/.test(msg) || /already running/i.test(msg))
            return "A job for this case is already running.";
        if (/404/.test(msg) && /case dir/i.test(msg))
            return "Case directory missing on disk.";
        if (/413/.test(msg))
            return "Request too large.";
        if (/no SESSION\.NAME/i.test(msg))
            return "Case has no SESSION.NAME file.";
        if (/no run\.local\.slurm/i.test(msg))
            return "Case has no Slurm script.";
        if (/mks failed/i.test(msg))
            return "Build failed — see logfile in case dir.";
        if (/sbatch failed/i.test(msg))
            return "sbatch rejected the job — check Slurm config.";
        return msg.length > 80 ? msg.slice(0, 80) + '…' : msg;
    }

    function renderJob(caseDiv, info) {
        const ui  = uiState(info.state);
        const lab = uiLabel(ui);
        caseDiv.replaceChildren();
        const pill = document.createElement('span');
        pill.className = 'pill ' + ui;
        pill.textContent = lab;
        caseDiv.appendChild(pill);
        if (info.jobid) {
            caseDiv.appendChild(document.createTextNode(' job '));
            const b = document.createElement('b');
            b.textContent = info.jobid;
            caseDiv.appendChild(b);
        }
        if (info.slurm_state && ui !== 'done' && ui !== 'failed') {
            caseDiv.appendChild(document.createTextNode(' · ' + info.slurm_state.toLowerCase()));
        }
        if (info.log) {
            const friendly = ui === 'failed' ? friendlyError(info.log) : info.log;
            caseDiv.appendChild(document.createTextNode(' · '));
            const s = document.createElement('span');
            s.className = 'job-log';
            s.textContent = friendly;
            caseDiv.appendChild(s);
        }
        caseDiv.classList.add('live');
        // toggle parent card class so .actions stays revealed during a live job
        const card = caseDiv.closest('.card');
        if (card) {
            const terminal = (ui === 'done' || ui === 'failed');
            card.classList.toggle('has-live-job', !terminal);
        }
    }

    async function pollJob(casePath, token, caseDiv) {
        // pause polling when tab hidden — server still tracks state, we'll catch up
        if (document.hidden) return;
        try {
            const r = await fetch('/api/jobs/' + encodeURIComponent(token));
            if (!r.ok) throw new Error('http ' + r.status);
            const info = await r.json();
            renderJob(caseDiv, info);
            const terminal = ['done', 'failed', 'error'].includes(info.state);
            if (terminal) {
                const existing = jobsByCase.get(casePath);
                if (existing && existing.pollId) clearInterval(existing.pollId);
                if (existing && existing.sibs) {
                    existing.sibs.forEach(function (b) {
                        b.disabled = b.dataset.staticDisabled === '1';
                    });
                }
                jobsByCase.delete(casePath);
            }
        } catch (e) {
            renderJob(caseDiv, { state: 'failed', log: 'poll: ' + e.message });
        }
    }

    // C1: confirm-before-fire on jobs estimated > LONG_JOB_THRESHOLD_SEC
    const LONG_JOB_THRESHOLD_SEC = 300; // 5 min
    function formatMinutes(sec) {
        if (sec < 60) return sec + ' s';
        if (sec < 3600) return Math.round(sec / 60) + ' min';
        return (sec / 3600).toFixed(1) + ' h';
    }
    function confirmIfLong(btn) {
        const elapsed = parseInt(btn.dataset.elapsed || '0', 10);
        const action  = btn.dataset.action;
        // recompile is bounded (~60s); only confirm resubmit on long-running cases
        if (action !== 'resubmit' || elapsed < LONG_JOB_THRESHOLD_SEC) return true;
        const casePath = btn.dataset.case;
        return window.confirm(
            'Resubmit ' + casePath + '?\n\n' +
            'Last run took ~' + formatMinutes(elapsed) + '. ' +
            'This will queue another job of similar length.'
        );
    }

    async function submitAction(btn) {
        const casePath = btn.dataset.case;
        const action   = btn.dataset.action || 'resubmit';
        const caseDiv  = document.querySelector(
            '.job-status[data-case="' + CSS.escape(casePath) + '"]'
        );
        const sibs = btn.parentElement.querySelectorAll('button.action-btn');
        sibs.forEach(function (b) { b.disabled = true; });
        renderJob(caseDiv, { state: 'starting', log: action + ' requested' });
        try {
            const r = await fetch('/api/run', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ case: casePath, action: action }),
            });
            let data;
            try { data = await r.json(); }
            catch (parseErr) { throw new Error('Unexpected token: server returned non-JSON'); }
            if (!r.ok) throw new Error(data.error || ('http ' + r.status));
            const initialLog = action === 'recompile' ? 'mks running...' : 'sbatch...';
            renderJob(caseDiv, { state: 'working', log: initialLog });
            const pollId = setInterval(function () {
                pollJob(casePath, data.token, caseDiv);
            }, 4000);
            const existing = jobsByCase.get(casePath);
            if (existing && existing.pollId) clearInterval(existing.pollId);
            jobsByCase.set(casePath, { token: data.token, pollId: pollId, sibs: sibs });
            pollJob(casePath, data.token, caseDiv);
        } catch (e) {
            renderJob(caseDiv, { state: 'failed', log: e.message });
            sibs.forEach(function (b) { b.disabled = b.dataset.staticDisabled === '1'; });
        }
    }

    document.querySelectorAll('button.action-btn').forEach(function (btn) {
        btn.addEventListener('click', function () {
            if (!confirmIfLong(btn)) return;  // user said no
            submitAction(btn);
        });
    });

    // ----- live queue panel (safe DOM, no innerHTML) -----
    const queueBody = document.getElementById('queue-body');
    const queueTick = document.getElementById('queue-tick');

    function makeEmpty(text, cls) {
        const d = document.createElement('div');
        d.className = cls || 'empty';
        d.textContent = text;
        return d;
    }

    function renderQueue(data) {
        queueBody.replaceChildren();
        if (data.error) {
            queueBody.appendChild(makeEmpty('slurm unreachable: ' + data.error, 'offline'));
            return;
        }
        const rows = data.jobs || [];
        if (rows.length === 0) {
            queueBody.appendChild(makeEmpty('queue idle — nothing running or pending'));
            return;
        }
        const table = document.createElement('table');
        const thead = document.createElement('thead');
        const htr = document.createElement('tr');
        ['JobID', 'Name', 'State', 'Time', 'Nodes', 'Reason'].forEach(function (h) {
            const th = document.createElement('th');
            th.textContent = h;
            htr.appendChild(th);
        });
        thead.appendChild(htr);
        table.appendChild(thead);
        const tbody = document.createElement('tbody');
        rows.forEach(function (j) {
            const tr = document.createElement('tr');
            ['jobid', 'name', 'state', 'time', 'nodes', 'reason'].forEach(function (k) {
                const td = document.createElement('td');
                td.textContent = j[k] != null ? j[k] : '';
                tr.appendChild(td);
            });
            tbody.appendChild(tr);
        });
        table.appendChild(tbody);
        queueBody.appendChild(table);
    }

    // Adaptive polling: faster when there's queue activity, slow when idle,
    // paused when tab is hidden. Saves CPU/network when nobody's looking.
    let queueIdleStreak = 0;
    let queuePollId = null;

    function nextQueueInterval() {
        if (document.hidden) return 30000;       // tab hidden
        if (queueIdleStreak >= 6) return 15000;  // 30+ s idle -> slow poll
        if (queueIdleStreak >= 2) return 8000;   // 10+ s idle -> medium
        return 4000;                              // active -> 4 s
    }
    const queueDot = document.getElementById('queue-status-dot');
    function setDot(state) {
        if (queueDot) queueDot.dataset.state = state;
    }
    let queueStaticMode = false;  // becomes true on first failure -> no backend
    function enterStaticMode() {
        if (queueStaticMode) return;
        queueStaticMode = true;
        const qp = document.getElementById('queue-panel');
        if (qp) qp.style.display = 'none';
        document.querySelectorAll('button.action-btn').forEach(function (b) {
            b.disabled = true;
            b.dataset.staticDisabled = '1';
            b.title = 'live actions need `python validation/serve.py`';
        });
    }
    async function pollQueueAdaptive() {
        try {
            const r = await fetch('/api/queue');
            const ct = r.headers.get('content-type') || '';
            if (!r.ok || !ct.includes('application/json')) throw new Error('no backend');
            const data = await r.json();
            renderQueue(data);
            queueTick.textContent = 'updated ' + new Date().toLocaleTimeString();
            const hasJobs = data.jobs && data.jobs.length > 0;
            queueIdleStreak = hasJobs ? 0 : (queueIdleStreak + 1);
            setDot(hasJobs ? 'active' : 'idle');
        } catch (e) {
            // First-time failure on page load -> assume static (Live Preview, file://)
            // and hide live-job UI cleanly instead of spamming console errors.
            enterStaticMode();
            return;  // stop polling
        }
        queuePollId = setTimeout(pollQueueAdaptive, nextQueueInterval());
    }
    pollQueueAdaptive();
    // Re-poll immediately when tab comes back into focus
    document.addEventListener('visibilitychange', function () {
        if (!document.hidden) {
            if (queuePollId) clearTimeout(queuePollId);
            queueIdleStreak = 0;
            pollQueueAdaptive();
        }
    });

    // ----- keyboard shortcuts -----
    const kbdHelp      = document.getElementById('kbd-help');
    const kbdHelpClose = document.getElementById('kbd-help-close');
    const sectionIds = Array.from(document.querySelectorAll('[id^="family-"]')).map(function (el) { return el.id; });

    function openHelp() {
        kbdHelp.classList.add('open');
        kbdHelp.setAttribute('aria-hidden', 'false');
        requestAnimationFrame(function () { kbdHelpClose.focus(); });
    }
    function closeHelp() {
        kbdHelp.classList.remove('open');
        kbdHelp.setAttribute('aria-hidden', 'true');
    }
    kbdHelpClose.addEventListener('click', closeHelp);
    kbdHelp.addEventListener('click', function (e) {
        if (e.target === kbdHelp) closeHelp();
    });

    let pendingGoto = false;
    let pendingGotoTimer = null;

    function jumpToSection(idx) {
        const id = sectionIds[idx];
        if (!id) return;
        const el = document.getElementById(id);
        if (!el) return;
        if (el.tagName === 'DETAILS' && !el.open) el.open = true;
        el.scrollIntoView({ behavior: 'smooth', block: 'start' });
        history.replaceState(null, '', '#' + id);
    }

    document.addEventListener('keydown', function (e) {
        // ignore when typing in inputs (none currently, but safe)
        const tag = (e.target && e.target.tagName) || '';
        if (tag === 'INPUT' || tag === 'TEXTAREA' || (e.target && e.target.isContentEditable)) return;

        if (e.key === '?' && !e.ctrlKey && !e.metaKey) {
            e.preventDefault();
            if (kbdHelp.classList.contains('open')) closeHelp(); else openHelp();
            return;
        }
        if (e.key === 'Escape' && kbdHelp.classList.contains('open')) {
            closeHelp();
            return;
        }

        // 'g' starts a two-key chord: g + digit -> section
        if (pendingGoto) {
            const n = parseInt(e.key, 10);
            if (!isNaN(n) && n >= 1 && n <= sectionIds.length) {
                jumpToSection(n - 1);
            }
            pendingGoto = false;
            clearTimeout(pendingGotoTimer);
            return;
        }
        if (e.key === 'g') {
            pendingGoto = true;
            pendingGotoTimer = setTimeout(function () { pendingGoto = false; }, 1200);
            return;
        }
        if (e.key === 't') {
            window.scrollTo({ top: 0, behavior: 'smooth' });
            return;
        }
        if (e.key === 'r') {
            if (queuePollId) clearTimeout(queuePollId);
            queueIdleStreak = 0;
            pollQueueAdaptive();
            return;
        }
    });

    // ----- filter chips -----
    (function () {
        var bar = document.getElementById('filter-bar');
        if (!bar) return;
        var clearBtn = document.getElementById('filter-clear');
        // active filters: object mapping group -> Set<value>
        var activeFilters = {};

        function applyFilters() {
            var cards = document.querySelectorAll('.card');
            cards.forEach(function (card) {
                var show = true;
                for (var group in activeFilters) {
                    var vals = activeFilters[group];
                    if (vals.size === 0) continue;
                    var cardVal = card.dataset[group] || '';
                    if (!vals.has(cardVal)) { show = false; break; }
                }
                card.classList.toggle('filtered-out', !show);
            });
        }

        bar.querySelectorAll('.filter-chip[data-filter-group]').forEach(function (chip) {
            chip.addEventListener('click', function () {
                var group = chip.dataset.filterGroup;
                var value = chip.dataset.filterValue;
                if (!activeFilters[group]) activeFilters[group] = new Set();
                if (activeFilters[group].has(value)) {
                    activeFilters[group].delete(value);
                    chip.classList.remove('active');
                } else {
                    activeFilters[group].add(value);
                    chip.classList.add('active');
                }
                applyFilters();
            });
        });

        if (clearBtn) {
            clearBtn.addEventListener('click', function () {
                activeFilters = {};
                bar.querySelectorAll('.filter-chip.active').forEach(function (fc) {
                    fc.classList.remove('active');
                });
                document.querySelectorAll('.card.filtered-out').forEach(function (fc) {
                    fc.classList.remove('filtered-out');
                });
            });
        }
    })();
})();
</script>""")
    parts.append('</body></html>')
    OUT.write_text('\n'.join(parts), encoding='utf-8')
    print(f'Wrote {OUT}')
    print(f'  total cases: {n_total} ({n_validated} validated, {n_deferred} deferred, {n_blocked} blocked)')
    print(f'  families: {len(ordered_fids)}')


if __name__ == '__main__':
    main()
