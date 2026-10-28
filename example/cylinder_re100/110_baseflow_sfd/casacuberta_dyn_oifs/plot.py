#!/usr/bin/env python3
"""Visualize the SFD dynamic + OIFS testcase.

This figure must compare against both non-OIFS baselines.  OIFS reduces
the number of timesteps substantially, but it is a case-specific acceleration;
the residual panel makes the speed/accuracy tradeoff explicit.

OUTPUTS: plot_ic.png, plot_residual.png, plot_residual_walltime.png, plot_bf.png
"""
from pathlib import Path
import sys, re

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import nekplot as nk

CASE_DIR = Path(__file__).resolve().parent
BASEFLOW_DIR = CASE_DIR.parent

CASES = [
    {
        'rel_path': Path('casacuberta') / 'residu.dat',
        'label': 'Casacuberta fixed',
        'color': '0.35',
        'linestyle': '--',
        'marker': 'x',
        'lw': 2.2,
        'zorder': 1,
    },
    {
        'rel_path': Path('casacuberta_dyn') / 'residu.dat',
        'label': 'Casacuberta dyn',
        'color': 'C1',
        'linestyle': '-',
        'marker': 'o',
        'lw': 1.4,
        'zorder': 2,
    },
    {
        'rel_path': Path('casacuberta_dyn_oifs') / 'residu.dat',
        'label': 'Casacuberta dyn OIFS',
        'color': 'C3',
        'linestyle': '-',
        'marker': '^',
        'lw': 0.9,
        'zorder': 3,
    },
]


def load_residual(path):
    data = np.genfromtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data[:, 0], data[:, 1]


def load_case(case):
    path = BASEFLOW_DIR / case['rel_path']
    if not path.exists():
        return None
    t, r = load_residual(path)
    return t, r


def plot_residual(ax, case, x, label_suffix=''):
    style = {
        'color': case['color'],
        'linestyle': case['linestyle'],
        'lw': case['lw'],
        'marker': case['marker'],
        'ms': 3.0,
        # about 12 markers per curve; the thickest curve (fixed) lies at the back
        'markevery': max(1, len(x) // 12),
        'zorder': case['zorder'],
    }
    label = f"{case['label']}{label_suffix}"
    ax.semilogy(x, case['residual'], label=label, **style)
    imin = int(np.argmin(case['residual']))
    ax.plot(x[imin], case['residual'][imin], marker='o', ms=2.4,
            color=case['color'], zorder=case['zorder'])


STEP = re.compile(r"^Step\s+\d+, t=\s*\S+, DT=\s*\S+, C=\s*\S+\s+(\S+)\s+\S+")


def load_walltime(case):
    """Wall time (s) of each time step, from the Step lines of the logfile next to residu.dat; None without a logfile."""
    log = (BASEFLOW_DIR / case['rel_path']).parent / 'logfile'
    if not log.exists():
        return None
    wall = [float(m.group(1)) for m in map(STEP.match, log.read_text(errors='replace').splitlines()) if m]
    return np.array(wall) if len(wall) >= case['steps'] else None


def load_available_cases():
    cases = []
    for case in CASES:
        data = load_case(case)
        if data is None:
            continue
        t, r = data
        entry = dict(case)
        entry['time'] = t
        entry['residual'] = r
        entry['steps'] = len(t)
        entry['wall'] = load_walltime(entry)
        cases.append(entry)
    return cases


def find_result_field():
    # Prefer the actual end-time checkpoint.  BFRe40_1cyl0.f00001 is the
    # shared initial condition and must not be plotted as the testcase result.
    for pattern in ('1cyl0.f*', 'BF_1cyl0.f*', 'BF*1cyl0.f*'):
        matches = nk.find_fields(pattern, CASE_DIR)
        if matches:
            return matches[0]
    return None


def _read_startfrom(par_path):
    """Return the startFrom filename from a .par file, or None."""
    try:
        with open(par_path) as f:
            for line in f:
                stripped = line.split('#')[0].strip()
                if stripped.lower().startswith('startfrom'):
                    val = stripped.split('=', 1)[1].strip()
                    return val
    except Exception:
        pass
    return None


def residual_gate(case_dir):
    """Stop gate of the run: the larger residualTol of the .par file (SFD and BoostConv stop below it)."""
    tols = [float(v) for f in case_dir.glob('*.par')
            for v in re.findall(r'^\s*residualTol\s*=\s*([0-9.eE+-]+)', f.read_text(), re.M)]
    return max(tols) if tols else None


def main():
    nk.configure_style()

    # PANEL: IC — initial condition from startFrom
    ic_fname = _read_startfrom(CASE_DIR / '1cyl.par')
    if ic_fname:
        ic_path = CASE_DIR / ic_fname
        if ic_path.exists():
            x, y, fields, _t = nk.read_field(ic_path)
            triang = nk.make_triangulation(x, y)
            umag = np.sqrt(fields['vx']**2 + fields['vy']**2)
            fig, ax_ic = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2, nk.COL_WIDTH * 0.4))
            cf = nk.tricontourf(ax_ic, triang, umag, levels=257,
                                cmap='Blues', vmin=0, vmax=1.5, extend='max')
            nk.inset_colorbar(ax_ic, cf, orientation='horizontal',
                              width='50%', height='5%', loc=9,
                              ticks=[0, 1.5], tick_labels=['0', '1.5'])
            ax_ic.set_aspect('equal')
            nk.add_cylinder_patches(ax_ic)
            ax_ic.set_xlim(-2, 20)
            ax_ic.set_ylim(-4, 4)
            ax_ic.set_xlabel(r'$x$', labelpad=-1)
            ax_ic.set_ylabel(r'$y$', labelpad=1)
            ax_ic.spines['right'].set_visible(False)
            ax_ic.spines['top'].set_visible(False)
            output = CASE_DIR / 'plot_ic.png'
            fig.savefig(output, dpi=600, bbox_inches='tight')
            print(f'Saved {output}')
            plt.close(fig)

    cases = load_available_cases()

    # PANEL: RESIDUAL — residual vs simulation time
    if cases:
        fig, ax_time = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.8,
                                                    nk.COL_WIDTH * 0.55))
        for case in cases:
            plot_residual(ax_time, case, case['time'],
                          label_suffix=f" ({case['steps']} steps)")
        gate = residual_gate(CASE_DIR)
        if gate:
            ax_time.axhline(gate, lw=0.5, c='k', ls=':', label=f'gate {gate:.0e}')
        ax_time.grid(True, which='both', lw=0.25, alpha=0.35)
        ax_time.set_xlabel(r'$t$')
        ax_time.set_ylabel(r'$\|r\|$')
        ax_time.legend(fontsize=4.8, loc='lower left',
                       handlelength=1.2, borderpad=0.25)
        output = CASE_DIR / 'plot_residual.png'
        fig.savefig(output, dpi=600, bbox_inches='tight')
        print(f'Saved {output}')
        plt.close(fig)

    # PANEL: RESIDUAL AGAINST WALL TIME — the cost of each variant; runs without a logfile are left out
    timed = [case for case in cases if case['wall'] is not None]
    if timed:
        fig, ax_wall = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.8,
                                                    nk.COL_WIDTH * 0.55))
        for case in timed:
            wall = case['wall'][:case['steps']]
            plot_residual(ax_wall, case, wall,
                          label_suffix=f" ({wall[-1]:.0f} s)")
        gate = residual_gate(CASE_DIR)
        if gate:
            ax_wall.axhline(gate, lw=0.5, c='k', ls=':', label=f'gate {gate:.0e}')
        ax_wall.grid(True, which='both', lw=0.25, alpha=0.35)
        ax_wall.set_xlabel('wall time (s)')
        ax_wall.set_ylabel(r'$\|r\|$')
        ax_wall.legend(fontsize=4.8, loc='upper right',
                       handlelength=1.2, borderpad=0.25)
        output = CASE_DIR / 'plot_residual_walltime.png'
        fig.savefig(output, dpi=400, bbox_inches='tight')
        print(f'Saved {output}')
        plt.close(fig)

    # PANEL: BF — converged base flow
    field_file = find_result_field()
    if field_file is not None:
        fig, ax_bf = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 2,
                                                  nk.COL_WIDTH * 0.4))
        x, y, fields, _time = nk.read_field(field_file)
        triang = nk.make_triangulation(x, y)
        umag = np.sqrt(fields['vx']**2 + fields['vy']**2)
        cf = nk.tricontourf(ax_bf, triang, umag, levels=257,
                            cmap='Blues', vmin=0, vmax=1.5, extend='max')
        nk.inset_colorbar(ax_bf, cf, orientation='horizontal',
                          width='50%', height='5%', loc=9,
                          ticks=[0, 1.5], tick_labels=['0', '1.5'])
        ax_bf.set_aspect('equal')
        nk.add_cylinder_patches(ax_bf)
        ax_bf.set_xlim(-2, 20)
        ax_bf.set_ylim(-4, 4)
        ax_bf.set_xlabel(r'$x$', labelpad=-1)
        ax_bf.set_ylabel(r'$y$', labelpad=1)
        ax_bf.spines['right'].set_visible(False)
        ax_bf.spines['top'].set_visible(False)
        output = CASE_DIR / 'plot_bf.png'
        fig.savefig(output, dpi=600, bbox_inches='tight')
        print(f'Saved {output}')
        plt.close(fig)


if __name__ == '__main__':
    main()
