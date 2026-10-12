#!/usr/bin/env python3
"""build_gallery.py - HTML gallery of the nekStab example stages.

The gallery is built from the stage folders of example/. A stage is a folder
with a .par file and a three-digit stage number in its path (000_dns,
310_stability_direct, ...). Each stage becomes one card:

  - the figures are the image files of the ref/ folder of the stage;
  - the reference table is ref/reference.json, the numbers that
    scripts/check_against_ref.py compares a new run with;
  - the one-line text is the "summary" of that file, and its "description" says what
    the reference values check.

Nothing is listed by hand. A stage without a reference, a summary, a description or
a figure is not drawn, and validation/tests/test_gallery_complete.py fails for it. The
page depends on tracked files only, not on the date or on the machine, so
validation/index.html in git is the output of this script.

Image src is emitted relative to validation/index.html (for example
"../example/cylinder_re100/000_dns/ref/snapshot.png"). The Pages artifact of
.github/workflows/pages-validation.yml ships the png and gif files of example/
next to validation/index.html, so the same paths work there.

The same run rewrites the table of verified stages in EXAMPLES.md between the markers
VERIFIED_START and VERIFIED_END: one row per drawn stage, with the summary of its reference.

Run:  python validation/build_gallery.py [output.html]
Open: xdg-open validation/index.html
"""
from __future__ import annotations

import html
import json
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'index.html'
EXAMPLE = ROOT.parent / 'example'
EXAMPLES_MD = ROOT.parent / 'EXAMPLES.md'

FIGURE_SUFFIXES = ('.png', '.gif')

# A folder that holds a .par file is a stage when its path has a number NNN_.
STAGE_NUMBER = re.compile(r'(?:^|/)(\d{3})_')
# Path parts that never belong to a stage: templates, reference folders, build
# folders, archives and the smoke test of the build (not a validation stage).
NOT_STAGE_PARTS = frozenset(['_templates', 'ref', 'obj', 'ci_test', '__pycache__'])
# Cases that stay in example/ as unvalidated work: no card, and the completeness test asks for nothing.
EXPERIMENTAL_FAMILIES = frozenset(['cylinder_re1m'])

# First digit of the stage number -> lane. The lanes are the filter chips and the
# sections inside a family, in page order.
STAGE_LANE = {
    '0': 'dns', '1': 'baseflow', '2': 'baseflow', '3': 'stability',
    '4': 'postproc', '5': 'otd', '6': 'modal',
}
LANES = [
    ('dns', 'DNS'),
    ('baseflow', 'Base flow'),
    ('stability', 'Stability'),
    ('postproc', 'Post-processing'),
    ('otd', 'OTD'),
    ('modal', 'Modal decomposition'),
]
# Base-flow stages split by the first digit of the stage number.
BASEFLOW_GROUPS = {
    '1': 'time-stepping methods',
    '2': 'Newton methods',
}

# Example folder -> (title, one sentence). The dictionary order is the page order.
# validation/tests/test_gallery_complete.py fails for a folder that is missing here.
FAMILIES: dict[str, tuple[str, str]] = {
    'cylinder_re100': (
        'Cylinder, Re = 100',
        'The reference family: the two-dimensional wake of a circular cylinder on one mesh of '
        '2128 elements. DNS, the steady-state solvers (SFD, BoostConv, DMT, Newton), direct and '
        'adjoint stability, transient growth, wavemaker, OTD and modal decomposition. BoostConv '
        'runs at Re = 50 and transient growth at Re = 40.',
    ),
    'cylinder_re180': (
        'Cylinder, Re = 180',
        'The periodic cylinder wake below the three-dimensional secondary instability: Newton for '
        'the periodic orbit, direct and adjoint Floquet stability, animation of the Floquet mode '
        'and OTD along the orbit. The Floquet analysis is two-dimensional and finds the phase '
        'mode only; modes A and B of Barkley and Henderson (1996) are three-dimensional.',
    ),
    'cylinder_re30_thermal': (
        'Heated cylinder, Re = 30',
        'The cylinder wake with an active temperature field (Richardson number 0.2): a DNS seed '
        'and the Newton steady state.',
    ),
    'moving_cylinder_re100': (
        'Moving cylinder, Re = 100',
        'A cylinder in forced transverse oscillation on a moving mesh: the lock-in of the wake, its '
        'periodic orbit by Newton, the Floquet multiplier of that orbit and the DMT filter bank on the orbit.',
    ),
    'flip_flop_re62': (
        'Flip-flop, Re = 62',
        'Two side-by-side cylinders just above the critical Reynolds number: the unstable periodic '
        'orbit, its direct and adjoint Floquet stability and the wavemaker.',
    ),
    'tpjet_re2000': (
        'Harmonically forced jet, Re = 2000',
        'The axisymmetric jet with harmonic forcing at St = 0.6: the forced periodic orbit by '
        'time-delayed feedback and by Newton, and its direct and adjoint Floquet stability.',
    ),
    'thermosyphon_ra500': (
        'Thermosyphon, Ra = 500',
        'The buoyancy-driven annular thermosyphon just above the critical Rayleigh number: DNS, the '
        'Newton steady state and direct stability.',
    ),
    'back_fstep_re500': (
        'Backward-facing step, Re = 500',
        'A stable but non-normal open shear flow: the Newton base flow and the optimal transient '
        'growth at the horizons of Blackburn, Barkley and Sherwin (2008).',
    ),
    'lid_driven_re3600': (
        'Lid-driven cavity, Re = 3600',
        'The two-dimensional cavity of aspect ratio 1.5 with a regularized lid: DNS, the Newton '
        'steady state, direct and adjoint stability.',
    ),
    'cubic_cavity_re1914': (
        'Cubic cavity, Re = 1914',
        'The three-dimensional lid-driven cavity just below the critical Reynolds number: the '
        'steady state, direct and adjoint stability, transient growth, the wavemaker and OTD.',
    ),
    'cubic_cavity_re1920': (
        'Cubic cavity, Re = 1920',
        'The cubic cavity just above the critical Reynolds number: the unstable steady state and '
        'its direct stability.',
    ),
    'cubic_cavity_re1950': (
        'Cubic cavity, Re = 1950',
        'The saturated limit cycle of the cubic cavity: its period, direct and adjoint Floquet '
        'stability, a DNS check of the unstable pair, POD and DMD of the cycle, and the DMT filter bank on the cycle.',
    ),
    'naca0012_re2000': (
        'NACA 0012, Re = 2000',
        'The wake of a NACA 0012 airfoil: DNS, the Newton steady state with sponge, direct and '
        'adjoint stability and the wavemaker.',
    ),
    'poiseuille_re5k': (
        'Poiseuille channel, Re = 5000',
        'OTD modes about the laminar channel profile, compared with theory.',
    ),
    'slot_fst_re495': (
        'Slot jet with free-stream turbulence, Re = 495',
        'DNS of a slot jet with a synthetic free-stream turbulence inflow.',
    ),
}

# userParam01 -> mode of nekStab. The keys are the codes of src/mode_codes.f90 divided by 100;
# validation/tests/test_gallery_complete.py compares the two lists.
MODE_NAMES = {
    '0': 'DNS',
    '0.1': 'linearized DNS',
    '1.1': 'selective frequency damping',
    '1.2': 'BoostConv',
    '1.3': 'dynamic mode tracking',
    '1.4': 'time-delayed feedback',
    '2': 'Newton, fixed point',
    '2.1': 'Newton, periodic orbit',
    '2.2': 'Newton, forced periodic orbit',
    '3.1': 'direct eigenmodes',
    '3.11': 'direct Floquet',
    '3.2': 'adjoint eigenmodes',
    '3.21': 'adjoint Floquet',
    '3.3': 'transient growth',
    '3.31': 'Floquet transient growth',
    '4': 'all post-processing',
    '4.1': 'energy budget',
    '4.11': 'Floquet energy budget',
    '4.2': 'wavemaker',
    '4.3': 'base-flow sensitivity',
    '4.41': 'force sensitivity, real part',
    '4.42': 'force sensitivity, imaginary part',
    '4.43': 'delta forcing',
    '4.5': 'mode animation',
    '4.51': 'animation with base-flow deformation',
    '4.52': 'Floquet mode animation',
    '5': 'OTD, frozen base flow',
    '5.1': 'OTD, base flow advanced by the DNS',
    '6': 'POD, DMD and SPOD',
    '6.1': 'POD',
    '6.2': 'DMD',
    '6.3': 'SPOD',
    '6.4': 'DMT filter bank',
}

# Names of the figure files of ref/ that a reader should not have to decode.
CAPTIONS = {
    'bf': 'Base flow',
    'baseflow': 'Base flow',
    'ic': 'Initial condition',
    'residual': 'Residual',
    'residual_vs_walltime': 'Residual against wall time',
    'spectrum': 'Spectrum',
    'mode': 'Leading mode',
    'direct_mode': 'Direct mode',
    'adjoint_mode': 'Adjoint mode',
    'adjoint': 'Adjoint mode',
    'wavemaker': 'Wavemaker',
    'snapshot': 'Snapshot',
    'signal': 'Probe signal',
    'fft': 'Spectrum of the signal',
    'phase': 'Phase portrait',
    'upo_snapshot': 'Orbit snapshot',
    'plot': 'Animation',
}


# Words of a figure file name that are written in capitals.
ACRONYMS = {'dns': 'DNS', 'pod': 'POD', 'dmd': 'DMD', 'spod': 'SPOD', 'otd': 'OTD', 'svd': 'SVD',
            'fft': 'FFT', 'ns': 'NS'}


@dataclass
class Stage:
    """One stage folder of example/ with the files that the gallery reads."""

    path: str                  # relative to example/, for example 'cylinder_re100/000_dns'
    family: str                # first part of the path
    number: str                # last stage number of the path, for example '310'
    lane: str                  # key of LANES
    re_tag: str                # 'Re 100', 'Ra 500' or ''
    mode: str                  # userParam01 as a key of MODE_NAMES, or ''
    summary: str               # one line: what the stage computes (ref/reference.json)
    description: str           # what the reference values check (ref/reference.json)
    quantities: dict           # quantities of ref/reference.json
    produced: dict             # produced of ref/reference.json
    figures: list[Path] = field(default_factory=list)
    drawable: bool = False     # a summary, a description, a quantity and a figure: enough for a card
    nelg: int | None = None    # elements of the mesh
    ranks: str = ''            # MPI ranks of the run behind the reference
    wall_s: float | None = None
    casename: str = ''
    has_slurm: bool = False
    problems: list[str] = field(default_factory=list)


def _stage_number(rel: str) -> str:
    found = STAGE_NUMBER.findall(rel)
    return found[-1] if found else ''


def _read_par(folder: Path) -> tuple[str, str, str]:
    """Return (case name, Reynolds number text, userParam01 key) of the .par files of a stage.

    A stage with several .par files (a DNS phase and a Newton phase) shows a value only
    when all files agree on it.
    """
    pars = sorted(folder.glob('*.par'))
    if not pars:
        return '', '', ''
    texts = [p.read_text(encoding='utf-8', errors='replace') for p in pars]
    # A negative viscosity is minus the Reynolds number in the .par convention.
    reynolds = {
        f'Re {-float(m.group(1)):g}'
        for m in (re.search(r'^\s*viscosity\s*=\s*(-[0-9.eE+]+)', t, re.M) for t in texts) if m
    }
    modes = {
        f'{float(m.group(1)):g}'
        for m in (re.search(r'^\s*userParam01\s*=\s*([0-9.]+)', t, re.M | re.I) for t in texts) if m
    }
    return (
        pars[0].stem,
        reynolds.pop() if len(reynolds) == 1 else '',
        modes.pop() if len(modes) == 1 else '',
    )


def _mesh_elements(folder: Path) -> int | None:
    """Return the number of elements from the ASCII header of the .re2 file."""
    re2 = next(iter(sorted(folder.glob('*.re2'))), None)
    if re2 is None:
        return None
    with re2.open('rb') as handle:
        tokens = handle.read(132).decode('ascii', errors='replace').split()
    return int(tokens[1]) if len(tokens) >= 2 and tokens[1].isdigit() else None


def _slurm_ranks(folder: Path) -> str:
    slurm = folder / 'run.local.slurm'
    if not slurm.is_file():
        return ''
    found = re.search(r'--ntasks[= ](\d+)', slurm.read_text(encoding='utf-8', errors='replace'))
    return found.group(1) if found else ''


def read_stage(folder: Path) -> Stage:
    """Read one stage folder. Missing parts are listed in Stage.problems."""
    rel = folder.relative_to(EXAMPLE).as_posix()
    problems: list[str] = []
    reference: dict = {}
    ref_dir = folder / 'ref'
    ref_file = ref_dir / 'reference.json'
    if ref_file.is_file():
        try:
            reference = json.loads(ref_file.read_text(encoding='utf-8'))
        except json.JSONDecodeError as exc:
            problems.append(f'ref/reference.json is not valid JSON ({exc})')
    else:
        problems.append('no ref/reference.json')
    summary = str(reference.get('summary', '')).strip()
    if ref_file.is_file() and not summary:
        problems.append('ref/reference.json has no summary')
    description = str(reference.get('description', '')).strip()
    if ref_file.is_file() and not description:
        problems.append('ref/reference.json has no description')
    quantities = reference.get('quantities', {})
    if ref_file.is_file() and not quantities:
        problems.append('ref/reference.json lists no quantity')
    figures = sorted(
        p for p in ref_dir.iterdir() if p.suffix.lower() in FIGURE_SUFFIXES
    ) if ref_dir.is_dir() else []
    if not figures:
        problems.append('no figure in ref/')
    if not (folder / 'README.md').is_file():
        problems.append('no README.md')
    casename, re_tag, mode = _read_par(folder)
    if not re_tag:
        # The viscosity of a buoyancy case is a Prandtl number: take Re or Ra from the folder name.
        named = re.search(r'_(re|ra)(\d+)(?:_|$)', rel.split('/')[0])
        re_tag = f'{named.group(1).capitalize()} {named.group(2)}' if named else ''
    produced = reference.get('produced', {})
    number = _stage_number(rel)
    wall = produced.get('wall_s')
    return Stage(
        path=rel,
        family=rel.split('/')[0],
        number=number,
        lane=STAGE_LANE.get(number[:1], 'other'),
        re_tag=re_tag,
        mode=mode,
        summary=summary,
        description=description,
        quantities=quantities,
        produced=produced,
        figures=figures,
        drawable=bool(summary and description and quantities and figures),
        nelg=_mesh_elements(folder),
        ranks=str(produced.get('ranks', produced.get('nproc', ''))) or _slurm_ranks(folder),
        wall_s=float(wall) if isinstance(wall, int | float) else None,
        casename=casename,
        has_slurm=(folder / 'run.local.slurm').is_file(),
        problems=problems,
    )


def _sort_key(stage: Stage) -> tuple[str, int, str]:
    # Inside one stage number the fixed-tolerance variant ("fp") comes first, then the rest by name.
    return (stage.number, 0 if stage.path.endswith('/fp') else 1, stage.path)


def discover_stages() -> list[Stage]:
    """Return every stage folder of example/, in page order, with or without figures."""
    stages: list[Stage] = []
    seen: set[str] = set()
    for par in sorted(EXAMPLE.rglob('*.par')):
        folder = par.parent
        rel = folder.relative_to(EXAMPLE).as_posix()
        parts = rel.split('/')
        if rel in seen or not STAGE_NUMBER.search(rel):
            continue
        if any(p in NOT_STAGE_PARTS or p.startswith('.') for p in parts):
            continue
        if parts[0] in EXPERIMENTAL_FAMILIES:
            continue
        seen.add(rel)
        stages.append(read_stage(folder))
    families = list(FAMILIES)
    stages.sort(key=lambda s: (
        families.index(s.family) if s.family in families else len(families), s.family,
        _sort_key(s),
    ))
    return stages


def drawn_stages(stages: list[Stage]) -> list[Stage]:
    """Return the stages that the page draws: the ones with a description, a quantity and a figure."""
    return [s for s in stages if s.drawable]


def _format_elapsed(sec: float) -> str:
    if sec < 60:
        return f'{sec:.0f} s'
    if sec < 3600:
        return f'{sec / 60:.1f} min'
    return f'{sec / 3600:.2f} h'


def _format_number(value: object) -> str:
    if isinstance(value, bool) or not isinstance(value, int | float):
        return str(value)
    return str(value) if isinstance(value, int) else f'{value:.8g}'


def _tolerance_text(spec: dict) -> str:
    if 'tol_rel' in spec:
        return f'±{spec["tol_rel"]:g} rel'
    if 'tol_abs' in spec:
        return f'±{spec["tol_abs"]:g}'
    return ''


def _produced_text(produced: dict) -> str:
    bits = []
    if produced.get('nekstab'):
        bits.append(f'nekStab {produced["nekstab"]}')
    ranks = produced.get('ranks', produced.get('nproc'))
    if ranks:
        bits.append(f'{ranks} ranks')
    for key in ('compiler', 'date'):
        if produced.get(key):
            bits.append(str(produced[key]))
    return 'produced with ' + ' · '.join(bits) if bits else ''


def _figure_caption(figure: Path) -> str:
    """Return the caption of a figure from its file name."""
    stem = figure.stem
    if stem in CAPTIONS:
        return CAPTIONS[stem]
    words = re.sub(r'(?<=[a-z])(?=\d)', ' ', stem.removeprefix('plot_')).replace('_', ' ').split()
    text = ' '.join(ACRONYMS.get(w, w) for w in words)
    return text[:1].upper() + text[1:]


def _img_url(png: Path) -> str:
    """URL of a figure, relative to validation/index.html."""
    return Path(os.path.relpath(png, ROOT)).as_posix()


def _versions(stages: list[Stage]) -> str:
    found = sorted({str(s.produced['nekstab']) for s in stages if s.produced.get('nekstab')})
    return ', '.join(found)


CSS = """
:root {
    /* surfaces */
    --bg:         #0e1116;
    --bg-card:    #1a1f2a;
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
.c-pending {
    display: inline-block;
    padding: 0.1em 0.55em;
    border-radius: 999px;
    margin-right: 0.25em;
    font-size: 0.95em;
}
.c-ok      { background: rgba(110, 224, 138, 0.14); color: var(--success); }
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
/* The case chips wrap to many rows on a phone, in portrait and in landscape. A sticky panel
   of that height covers most of the screen, so it scrolls with the page there. */
@media (max-width: 599px), (max-height: 500px) {
    .section-nav { position: static; }
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

/* Intro line, family blurb and the reference table of a card */
p.intro {
    color: var(--fg-mute);
    font-size: 0.88em;
    max-width: 62em;
    margin: 0.6em 0 0;
    text-wrap: pretty;
}
p.intro code {
    font-family: ui-monospace, SF Mono, Monaco, Cascadia Code, monospace;
    color: var(--accent);
    font-size: 0.95em;
}
.card details.ref {
    margin-top: 0.45em;
    font-size: 0.82em;
    color: var(--fg-mute);
}
.card details.ref > summary {
    cursor: pointer;
    color: var(--accent);
    user-select: none;
}
.card details.ref .ref-desc {
    margin: 0.5em 0 0.4em;
    color: var(--fg);
    line-height: 1.4;
}
.card table.ref-table {
    width: 100%;
    border-collapse: collapse;
    font-variant-numeric: tabular-nums;
}
.card table.ref-table td {
    padding: 0.15em 0.5em 0.15em 0;
    vertical-align: top;
    border-top: 1px solid var(--border-faint);
}
.card table.ref-table td.name {
    font-family: ui-monospace, SF Mono, Monaco, monospace;
    font-size: 0.92em;
    word-break: break-word;
}
.card table.ref-table td.num {
    text-align: right;
    white-space: nowrap;
    color: var(--fg);
}
"""

LIGHTBOX_AND_HELP = (
    '<div class="lightbox" id="lightbox" hidden role="dialog" aria-modal="true" aria-label="Figure viewer">'
    '<div class="lightbox-backdrop"></div>'
    '<button class="lightbox-close" aria-label="Close">&times;</button>'
    '<img class="lightbox-img" alt="">'
    '<div class="lightbox-caption"></div>'
    '<div class="lightbox-nav-hint">&#8592; prev figure &middot; &#8594; next figure &middot; &#8593; prev stage '
    '&middot; &#8595; next stage &middot; esc close</div>'
    '</div>'
    '<div class="kbd-help" id="kbd-help" role="dialog" aria-modal="true"'
    ' aria-labelledby="kbd-help-title" aria-hidden="true">'
    '<div class="kbd-panel">'
    '<h3 id="kbd-help-title">Keyboard shortcuts</h3>'
    '<dl>'
    '<dt><kbd>?</kbd></dt><dd>Toggle this help</dd>'
    '<dt><kbd>Esc</kbd></dt><dd>Close lightbox / help</dd>'
    '<dt><kbd>g</kbd> then <kbd>1</kbd>&hellip;<kbd>9</kbd></dt>'
    '<dd>Jump to N-th case</dd>'
    '<dt><kbd>t</kbd></dt><dd>Scroll to top</dd>'
    '<dt><kbd>r</kbd></dt><dd>Refresh queue panel</dd>'
    '</dl>'
    '<button id="kbd-help-close" type="button">close</button>'
    '</div>'
    '</div>'
)

SCRIPT = r"""<script>
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
        if (/403/.test(msg) && /stage folder/i.test(msg))
            return "Not a stage folder of example/ — can't resubmit.";
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
        // Without the server the buttons cannot work, so a reader of the static page does not see them.
        document.querySelectorAll('.card .actions, .card .job-status').forEach(function (el) {
            el.style.display = 'none';
        });
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
</script>"""


def _reference_table(stage: Stage) -> str:
    """Return the collapsed table of the reference values of a stage."""
    rows = []
    for name, spec in stage.quantities.items():
        note = html.escape(str(spec.get('note', '')), quote=True)
        source = html.escape(str(spec.get('source', '')), quote=True)
        rows.append(
            f'<tr><td class="name" title="{source}">{html.escape(name)}</td>'
            f'<td class="num" title="{note}">{html.escape(_format_number(spec.get("value")))}</td>'
            f'<td class="num">{html.escape(_tolerance_text(spec))}</td></tr>'
        )
    return (
        '<details class="ref"><summary>Reference values</summary>'
        f'<p class="ref-desc">{html.escape(stage.description)}</p>'
        f'<table class="ref-table"><tbody>{"".join(rows)}</tbody></table>'
        '</details>'
    )


def render_card(stage: Stage) -> str:
    """Return the HTML of one stage card."""
    path = html.escape(stage.path, quote=True)
    lane = html.escape(stage.lane, quote=True)
    re_tag = html.escape(stage.re_tag, quote=True)
    n_checks = len(stage.quantities)
    meta_bits = [b for b in (stage.re_tag, f'mode {stage.mode} · {MODE_NAMES[stage.mode]}'
                             if stage.mode in MODE_NAMES else '') if b]
    meta_text = ' · '.join(meta_bits)

    chunks: list[str] = [
        f'<div class="card" data-case-id="{path}" data-family="{html.escape(stage.family, quote=True)}"'
        f' data-lane="{lane}" data-re="{re_tag}">',
        '<div class="label label-top">',
        f'<div class="copyable path" data-action="copy" data-copy="{path}">{html.escape(stage.path)}</div>',
    ]
    if meta_text:
        chunks.append(
            f'<div class="copyable meta-row" data-action="copy" data-copy="{html.escape(meta_text, quote=True)}">'
            f'{html.escape(meta_text)}</div>'
        )
    chunks.append(
        f'<div class="copyable expected" data-action="copy" data-copy="{html.escape(stage.summary, quote=True)}">'
        f'{html.escape(stage.summary)}</div>'
    )
    chunks.append(
        '<div class="evidence-row">'
        f'<span class="badge-status ev-image" title="numbers that scripts/check_against_ref.py compares a run with">'
        f'{n_checks} reference check{"" if n_checks == 1 else "s"}</span></div>'
    )
    chunks.append(_reference_table(stage))
    produced = _produced_text(stage.produced)
    if produced:
        chunks.append(f'<div class="provenance-row">{html.escape(produced)}</div>')
    stat_bits: list[str] = []
    if stage.nelg is not None:
        stat_bits.append(f'<span class="stat"><b>{stage.nelg}</b> el</span>')
    if stage.ranks:
        stat_bits.append(f'<span class="stat"><b>{html.escape(stage.ranks)}</b> mpi</span>')
    if stage.wall_s is not None:
        stat_bits.append(f'<span class="stat"><b>{html.escape(_format_elapsed(stage.wall_s))}</b></span>')
    if stat_bits:
        chunks.append('<div class="stats">' + ''.join(stat_bits) + '</div>')
    chunks.append('</div>')

    chunks.append('<div class="thumbs-strip">')
    for index, figure in enumerate(stage.figures):
        src = html.escape(_img_url(figure), quote=True)
        caption = html.escape(f'{stage.path} · {_figure_caption(figure)}', quote=True)
        chunks.append(
            f'<button class="thumb" data-case-id="{path}" data-panel-index="{index}"'
            f' data-panel-of="{len(stage.figures)}" data-img-src="{src}" data-caption="{caption}"'
            f' aria-label="Open figure {index + 1} of {len(stage.figures)}">'
            f'<img class="thumb-img" loading="lazy" decoding="async" src="{src}" alt="{caption}">'
            f'<div class="thumb-cap">{html.escape(_figure_caption(figure))}</div></button>'
        )
    chunks.append('</div>')

    elapsed = f' data-elapsed="{int(stage.wall_s)}"' if stage.wall_s is not None else ''
    recompile = '' if stage.casename else ' disabled data-static-disabled="1" title="stage has no .par file"'
    resubmit = '' if stage.has_slurm else (
        ' disabled data-static-disabled="1" title="stage has no Slurm script, so it cannot be submitted"'
    )
    chunks.append('<div class="actions">')
    chunks.append(f'<span class="casename">{html.escape(stage.casename)}</span>')
    chunks.append('<div class="btn-group">')
    chunks.append(
        f'<button class="action-btn" data-action="recompile" aria-label="Recompile {path}"'
        f' data-case="{path}"{elapsed}{recompile}>recompile</button>'
    )
    chunks.append(
        f'<button class="action-btn primary" data-action="resubmit" aria-label="Resubmit {path}"'
        f' data-case="{path}"{elapsed}{resubmit}>resubmit</button>'
    )
    chunks.append('</div></div>')
    chunks.append(f'<div class="job-status" data-case="{path}" aria-live="polite" aria-atomic="true"></div>')
    chunks.append('</div>')
    return ''.join(chunks)


def _grid(stages: list[Stage]) -> str:
    return '<div class="grid">' + ''.join(render_card(s) for s in stages) + '</div>'


def _ladder(stages: list[Stage]) -> str:
    """Numbered steps, one full-width card per stage, for variants that follow each other."""
    steps = [
        f'<div class="ladder-step"><div class="ladder-step__num">{i}</div>'
        f'<div class="ladder-step__card">{render_card(s)}</div></div>'
        for i, s in enumerate(stages, start=1)
    ]
    return '<div class="ladder">' + ''.join(steps) + '</div>'


def _mode_label(text: str, count: int) -> str:
    return (
        f'<div class="mode-label" style="margin-top:0.8em">{html.escape(text)}'
        f'<span class="counts">{count}</span></div>'
    )


def render_family(family: str, stages: list[Stage]) -> str:
    """Return the section of one example folder: its stages grouped by lane."""
    title, blurb = FAMILIES.get(family, (family.replace('_', ' '), ''))
    parts = [
        f'<section class="geom" id="family-{html.escape(family, quote=True)}" style="scroll-margin-top:4.5em">',
        f'<h2>{html.escape(title)} <span class="counts"><span class="c-ok">{len(stages)} '
        f'stage{"" if len(stages) == 1 else "s"}</span></span></h2>',
    ]
    if blurb:
        parts.append(
            '<p style="color:var(--fg-mute);font-size:0.85em;font-style:italic;margin:0 0 1em">'
            f'{html.escape(blurb)}</p>'
        )
    for lane, lane_label in LANES:
        group = [s for s in stages if s.lane == lane]
        if not group:
            continue
        if lane != 'baseflow':
            parts.append(_mode_label(lane_label, len(group)))
            parts.append(_grid(group))
            continue
        subgroups = [(key, [s for s in group if s.number[:1] == key]) for key in BASEFLOW_GROUPS]
        for key, members in subgroups:
            if not members:
                continue
            label = lane_label if len([m for _, m in subgroups if m]) == 1 else (
                f'{lane_label}: {BASEFLOW_GROUPS[key]}'
            )
            parts.append(_mode_label(label, len(members)))
            parts.append(_ladder(members) if len(members) >= 2 else _grid(members))
    parts.append('</section>')
    return ''.join(parts)


VERIFIED_START = '<!-- verified:start (written by validation/build_gallery.py from the summary of each ref/reference.json) -->'
VERIFIED_END = '<!-- verified:end -->'


def verified_table(stages: list[Stage]) -> str:
    """Return the Markdown table of the drawn stages for EXAMPLES.md."""
    rows = ['| Stage | Re / Ra | Mode | Result |', '|---|---|---|---|']
    for stage in stages:
        mode = f'{MODE_NAMES[stage.mode]} ({stage.mode})' if stage.mode in MODE_NAMES else ''
        summary = ' '.join(stage.summary.split()).replace('|', '\\|')
        rows.append(f'| `example/{stage.path}/` | {stage.re_tag} | {mode} | {summary} |')
    return '\n'.join(rows)


def examples_md_with_table(text: str, stages: list[Stage]) -> str:
    """Return the text of EXAMPLES.md with the table between the markers replaced."""
    start = text.index(VERIFIED_START) + len(VERIFIED_START)
    end = text.index(VERIFIED_END)
    return text[:start] + '\n' + verified_table(stages) + '\n' + text[end:]


def build_page() -> str:
    """Return the HTML of the gallery for the stage folders in example/."""
    stages = drawn_stages(discover_stages())
    families = [f for f in FAMILIES if any(s.family == f for s in stages)]
    families += sorted({s.family for s in stages} - set(FAMILIES))
    n_checks = sum(len(s.quantities) for s in stages)
    versions = _versions(stages)

    parts: list[str] = [
        '<!doctype html>',
        '<html lang="en"><head>',
        '<meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1">',
        '<meta name="color-scheme" content="dark">',
        '<meta name="robots" content="noindex, nofollow">',
        f'<title>nekStab validation gallery — {len(stages)} stages</title>',
        f'<style>{CSS}</style>',
        '</head><body>',
        '<header class="page-head">',
        '<h1>nekStab validation gallery</h1>',
        '<div class="meta">'
        + (f'references produced with nekStab <code>{html.escape(versions)}</code> &middot; ' if versions else '')
        + f'{len(stages)} stages in {len(families)} cases &middot; {n_checks} reference checks'
        '</div>',
        '<p class="intro">Each card is one stage folder of <code>example/</code>. The figures are the files of '
        'its <code>ref/</code> folder. The reference values are the numbers that '
        '<code>scripts/check_against_ref.py</code> compares a new run with. Open a figure to enlarge it; '
        'the arrow keys move between figures.</p>',
        '</header>',
        '<div class="filter-bar" id="filter-bar" aria-label="Filter cards">'
        '<div class="filter-group"><span class="filter-label">Lane:</span>'
        + ''.join(
            f'<button class="filter-chip" aria-label="Filter by {label} lane" data-filter-group="lane"'
            f' data-filter-value="{lane}">{html.escape(label.lower())}</button>'
            for lane, label in LANES
        )
        + '</div>'
        '<button class="filter-chip clear-btn" aria-label="Clear all filters" id="filter-clear">Clear filters</button>'
        '</div>',
        '<div class="queue-panel" id="queue-panel"><div class="qhead"><span class="qtitle">'
        '<span class="status-dot" id="queue-status-dot" data-state="connecting" aria-hidden="true"></span>'
        'Running jobs</span><span class="tick" id="queue-tick">connecting&hellip;</span></div>'
        '<div id="queue-body"><div class="empty">waiting for first poll&hellip;</div></div></div>',
        '<nav class="section-nav" id="family-nav" aria-label="Flow family navigation">'
        + ''.join(
            f'<a class="nav-chip" href="#family-{html.escape(f, quote=True)}"'
            f' aria-label="Jump to {html.escape(FAMILIES.get(f, (f,))[0], quote=True)}">'
            f'{html.escape(FAMILIES.get(f, (f,))[0])}'
            f'<span class="count">{sum(1 for s in stages if s.family == f)}</span></a>'
            for f in families
        )
        + '</nav>',
    ]
    for family in families:
        parts.append(render_family(family, [s for s in stages if s.family == family]))
    parts.append(LIGHTBOX_AND_HELP)
    parts.append(SCRIPT)
    parts.append('</body></html>')
    return '\n'.join(parts) + '\n'


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    out = Path(argv[0]) if argv else OUT
    stages = discover_stages()
    out.write_text(build_page(), encoding='utf-8')
    drawn = drawn_stages(stages)
    EXAMPLES_MD.write_text(examples_md_with_table(EXAMPLES_MD.read_text(encoding='utf-8'), drawn), encoding='utf-8')
    print(f'Wrote {out}')
    print(f'  stages drawn: {len(drawn)} of {len(stages)} in {len({s.family for s in drawn})} cases')
    for stage in stages:
        for problem in stage.problems:
            print(f'  not drawn or incomplete: {stage.path}: {problem}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
