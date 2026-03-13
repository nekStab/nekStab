#!/usr/bin/env python3
"""collect_plots.py: Gather all plot.png into validation/figures/ with descriptive names.

Walks example directories, finds plot.png files, extracts Re/Ra from .par files,
and copies them to a single validation/figures/ folder with names like:
    cylinder_baseflow_sfd_Re50.png
    thersyphon_stability_direct_Ra500.png

Can also (re-)generate plots before collecting by running each plot.py.

OUTPUTS: validation/figures/*.png
USAGE:
    python example/collect_plots.py               # collect existing plot.png
    python example/collect_plots.py --generate     # run plot.py then collect
    python example/collect_plots.py --list         # show what would be collected
"""
from pathlib import Path
import os
import re
import shutil
import subprocess
import sys

EXAMPLE_DIR = Path(__file__).resolve().parent
NEKSTAB_ROOT = EXAMPLE_DIR.parent
OUTPUT_DIR = NEKSTAB_ROOT / 'validation' / 'figures'


def _find_python():
    """Find a Python interpreter that has pymech installed."""
    # Try current interpreter first
    candidates = [sys.executable]
    # Try common venv locations
    for venv in [Path.home() / '.venv', NEKSTAB_ROOT / '.venv',
                 EXAMPLE_DIR / '.venv']:
        p = venv / 'bin' / 'python'
        if p.exists():
            candidates.insert(0, str(p))
    for py in candidates:
        try:
            r = subprocess.run([py, '-c', 'import pymech'],
                               capture_output=True, timeout=5)
            if r.returncode == 0:
                return py
        except Exception:
            continue
    return sys.executable  # fallback


PYTHON = _find_python()

# Directories to skip (3D cases, venvs, non-case dirs)
SKIP_DIRS = {'.venv', '__pycache__', 'modal'}


def find_par_file(case_dir):
    """Find the .par file in case_dir or nearest parent under example/."""
    d = case_dir
    while d != EXAMPLE_DIR.parent:
        pars = list(d.glob('*.par'))
        if pars:
            return pars[0]
        d = d.parent
    return None


def extract_params(par_file):
    """Extract Re and/or Ra from a .par file."""
    if par_file is None:
        return {}
    text = par_file.read_text()
    params = {}

    # Viscosity: negative = 1/Re
    m = re.search(r'(?i)viscosity\s*=\s*([-.\deE]+)', text)
    if m:
        v = float(m.group(1))
        if v < 0:
            params['Re'] = int(abs(v))

    # userParam06: Rayleigh number (thermal cases)
    m = re.search(r'(?i)userParam06\s*=\s*([-.\deE]+)', text)
    if m:
        val = float(m.group(1))
        if val > 0:
            params['Ra'] = int(val)

    # Drop Ra=0 (non-thermal cases sometimes have it)
    if params.get('Ra') == 0:
        del params['Ra']

    return params


def make_name(case_dir, params):
    """Build descriptive filename from directory path and parameters.

    example/cylinder/baseflow/sfd  +  {Re: 50}  →  cylinder_baseflow_sfd_Re50.png
    """
    rel = case_dir.relative_to(EXAMPLE_DIR)
    parts = [p for p in rel.parts if p not in SKIP_DIRS]
    base = '_'.join(parts)

    # Append parameters (skip if already in dir name)
    for key in ('Re', 'Ra'):
        if key in params:
            tag = f'{key}{params[key]}'
            if tag.lower() not in base.lower():
                base += f'_{tag}'

    return base + '.png'


def collect_plots(generate=False, list_only=False):
    """Main collection logic."""
    # Find all plot.py files (our scripts, not third-party)
    plot_scripts = sorted(EXAMPLE_DIR.rglob('plot.py'))
    plot_scripts = [p for p in plot_scripts
                    if not any(skip in p.parts for skip in SKIP_DIRS)]

    if not list_only:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    collected = []
    skipped = []

    for script in plot_scripts:
        case_dir = script.parent
        plot_png = case_dir / 'plot.png'

        # Optionally run plot.py to generate/refresh plot.png
        if generate:
            print(f'  Running {script.relative_to(NEKSTAB_ROOT)} ...', end=' ',
                  flush=True)
            try:
                result = subprocess.run(
                    [PYTHON, str(script)],
                    capture_output=True, text=True, timeout=120,
                    cwd=str(case_dir),
                    env={**os.environ, 'MPLBACKEND': 'Agg'}
                )
                if result.returncode == 0:
                    print('OK')
                else:
                    err = result.stderr.strip().split('\n')[-1] if result.stderr else 'unknown'
                    print(f'FAIL ({err})')
            except subprocess.TimeoutExpired:
                print('TIMEOUT')
            except Exception as e:
                print(f'ERROR ({e})')

        if not plot_png.exists():
            skipped.append(case_dir.relative_to(EXAMPLE_DIR))
            continue

        # Extract parameters and build name
        par = find_par_file(case_dir)
        params = extract_params(par)
        name = make_name(case_dir, params)

        if list_only:
            print(f'  {name:55s} ← {case_dir.relative_to(EXAMPLE_DIR)}/')
        else:
            dest = OUTPUT_DIR / name
            shutil.copy2(plot_png, dest)
            collected.append(name)

    if not list_only:
        print(f'\nCollected {len(collected)} plots → {OUTPUT_DIR.relative_to(NEKSTAB_ROOT)}/')
        for name in collected:
            print(f'  {name}')
        if skipped:
            print(f'\nSkipped {len(skipped)} (no plot.png):')
            for s in skipped:
                print(f'  {s}/')
    else:
        if skipped:
            print(f'\nNo plot.png yet ({len(skipped)}):')
            for s in skipped:
                print(f'  {s}/')


if __name__ == '__main__':
    args = sys.argv[1:]
    generate = '--generate' in args or '-g' in args
    list_only = '--list' in args or '-l' in args

    if list_only:
        print('Would collect:')
    elif generate:
        print('Generating and collecting plots...')
    else:
        print('Collecting existing plots...')

    collect_plots(generate=generate, list_only=list_only)
