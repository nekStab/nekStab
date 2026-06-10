#!/usr/bin/env python3
"""Visualize the SFD dynamic + OIFS testcase.

This figure must compare against both non-OIFS baselines.  OIFS reduces
the number of timesteps substantially, but it is a case-specific acceleration;
the residual and measured-cost panels make the speed/accuracy tradeoff explicit.

OUTPUTS: plot_ic.png, plot_residual.png, plot_residual_vs_walltime.png, plot_bf.png
"""
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import nekplot as nk

CASE_DIR = Path(__file__).resolve().parent
BASEFLOW_DIR = CASE_DIR.parent
N_RANKS = 8

# Elapsed times are parsed dynamically from each case's logfile.
# The OIFS run takes many fewer timesteps, but each OIFS step is more expensive
# than a standard BDF step, so a timestep count alone exaggerates the practical gain.
CASES = [
    {
        'rel_path': Path('casacuberta') / 'residu.dat',
        'label': 'Casacuberta fixed',
        'color': '0.35',
        'linestyle': '--',
        'zorder': 1,
    },
    {
        'rel_path': Path('casacuberta_dyn') / 'residu.dat',
        'label': 'Casacuberta dyn',
        'color': 'C1',
        'linestyle': '-',
        'zorder': 2,
    },
    {
        'rel_path': Path('casacuberta_dyn_oifs') / 'residu.dat',
        'label': 'Casacuberta dyn OIFS',
        'color': 'C3',
        'linestyle': '-',
        'zorder': 3,
    },
]


def parse_elapsed(logfile_path):
    """Parse 'total elapsed time : X sec' from a Nek5000 logfile."""
    try:
        with open(logfile_path) as f:
            for line in f:
                if 'total elapsed time' in line.lower():
                    parts = line.split(':')
                    return float(parts[-1].split()[0])
    except Exception:
        pass
    return None


# Populate elapsed from logfiles at import time so load_available_cases can use it
for _case in CASES:
    _logfile = BASEFLOW_DIR / _case['rel_path'].parent / 'logfile'
    _case['elapsed'] = parse_elapsed(_logfile)


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


def wall_cost_axis(n_samples, elapsed_seconds):
    # One residual sample is written per timestep.  We distribute the measured
    # total wall time over the samples so the plot reports actual run cost
    # instead of only physical time or timestep count.
    return np.linspace(elapsed_seconds / n_samples, elapsed_seconds, n_samples)


def plot_residual(ax, case, x, label_suffix=''):
    style = {
        'color': case['color'],
        'linestyle': case['linestyle'],
        'lw': 0.65 if case['label'] == 'SFD fixed' else 0.85,
        'zorder': case['zorder'],
    }
    label = f"{case['label']}{label_suffix}"
    ax.semilogy(x, case['residual'], label=label, **style)
    imin = int(np.argmin(case['residual']))
    ax.plot(x[imin], case['residual'][imin], marker='o', ms=2.4,
            color=case['color'], zorder=case['zorder'])


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
        if case['elapsed'] is not None:
            entry['cost'] = wall_cost_axis(len(t), case['elapsed'])
        else:
            entry['cost'] = None
        cases.append(entry)
    return cases


def add_core_min_axis(ax):
    sec_to_core_min = lambda seconds: seconds * N_RANKS / 60.0
    core_min_to_sec = lambda core_min: core_min * 60.0 / N_RANKS
    top = ax.secondary_xaxis('top', functions=(sec_to_core_min,
                                               core_min_to_sec))
    top.set_xlabel('core-min')
    top.tick_params(labelsize=5.5, pad=1)


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
        ax_time.axhline(1e-9, lw=0.5, c='k', ls=':', label='1e-9 target')
        ax_time.grid(True, which='both', lw=0.25, alpha=0.35)
        ax_time.set_xlabel(r'$t$')
        ax_time.set_ylabel(r'$\|r\|$')
        ax_time.legend(fontsize=4.8, loc='lower left',
                       handlelength=1.2, borderpad=0.25)
        output = CASE_DIR / 'plot_residual.png'
        fig.savefig(output, dpi=600, bbox_inches='tight')
        print(f'Saved {output}')
        plt.close(fig)

    # PANEL: RESIDUAL_VS_WALLTIME — only if elapsed times are available
    cases_with_cost = [c for c in cases if c.get('cost') is not None]
    if cases_with_cost:
        fig, ax_cost = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.85,
                                                    nk.COL_WIDTH * 0.55))
        for case in cases_with_cost:
            plot_residual(ax_cost, case, case['cost'],
                          label_suffix=f" ({case['elapsed']:.0f} s)")
        ax_cost.axhline(1e-9, lw=0.5, c='k', ls=':', label='1e-9 target')
        ax_cost.grid(True, which='both', lw=0.25, alpha=0.35)
        ax_cost.set_xlabel('wall time on 8 ranks (s)')
        ax_cost.set_ylabel(r'$\|r\|$')
        add_core_min_axis(ax_cost)
        output = CASE_DIR / 'plot_residual_vs_walltime.png'
        fig.savefig(output, dpi=600, bbox_inches='tight')
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
