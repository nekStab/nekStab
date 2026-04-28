#!/usr/bin/env python3
"""Plot residual and scheduler histories for the SFD/OIFS comparison case."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    'text.usetex': False,
    'font.size': 8,
    'legend.fontsize': 6,
    'legend.handlelength': 1.4,
})
plt.style.use('seaborn-v0_8-white')

CASE_DIR = Path(__file__).resolve().parent
BASEFLOW_DIR = CASE_DIR.parent
N_RANKS = 8

# Measured elapsed wall times from the latest matching 8-rank Re=50 runs.  They
# are intentionally stored next to the plotting code so the cost comparison is
# reproducible from the committed histories.  If the cases are rerun, update
# these numbers together with the README.
CASES = [
    {
        'rel_path': Path('sfd') / 'residu.dat',
        'label': 'SFD fixed',
        'color': '0.35',
        'linestyle': '--',
        'elapsed': 231.736,
        'zorder': 1,
    },
    {
        'rel_path': Path('sfd_dyn') / 'residu.dat',
        'label': 'SFD dyn',
        'color': 'C1',
        'linestyle': '-',
        'elapsed': 199.891,
        'zorder': 2,
    },
    {
        'rel_path': Path('sfd_dyn_oifs') / 'residu.dat',
        'label': 'SFD dyn OIFS',
        'color': 'C3',
        'linestyle': '-',
        'elapsed': 110.790,
        'zorder': 3,
    },
]


def load_table(path):
    data = np.genfromtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def wall_cost_axis(n_samples, elapsed_seconds):
    # One residual sample is written per timestep.  Mapping sample index to the
    # measured total wall time shows the actual 8-rank cost of each method.
    return np.linspace(elapsed_seconds / n_samples, elapsed_seconds, n_samples)


def load_case(case):
    path = BASEFLOW_DIR / case['rel_path']
    if not path.exists():
        return None
    data = load_table(path)
    out = dict(case)
    out['time'] = data[:, 0]
    out['residual'] = data[:, 1]
    out['steps'] = len(data)
    out['cost'] = wall_cost_axis(len(data), case['elapsed'])
    return out


def plot_residual(ax, case, x, label_suffix):
    lw = 0.65 if case['label'] == 'SFD fixed' else 0.85
    ax.semilogy(x, case['residual'], lw=lw, color=case['color'],
                ls=case['linestyle'], zorder=case['zorder'],
                label=f"{case['label']} {label_suffix}")
    imin = int(np.argmin(case['residual']))
    ax.plot(x[imin], case['residual'][imin], marker='o', ms=2.5,
            color=case['color'], zorder=case['zorder'])


def add_core_min_axis(ax):
    sec_to_core_min = lambda seconds: seconds * N_RANKS / 60.0
    core_min_to_sec = lambda core_min: core_min * 60.0 / N_RANKS
    top = ax.secondary_xaxis('top', functions=(sec_to_core_min,
                                               core_min_to_sec))
    top.set_xlabel('core-min')
    top.tick_params(labelsize=6, pad=1)


def panel_label(ax, label):
    ax.text(0.03, 0.96, f'({label})', transform=ax.transAxes,
            ha='left', va='top', fontsize=8, fontweight='bold',
            bbox=dict(facecolor='white', edgecolor='none',
                      alpha=0.78, pad=1.1))


def maybe_plot_scheduler(ax, rel_path, label, color):
    path = BASEFLOW_DIR / rel_path
    if not path.exists():
        return
    data = load_table(path)
    t = data[:, 0]
    requested = data[:, 3]
    used = data[:, 4]
    ax.semilogy(t, requested, lw=0.65, color=color, ls='--',
                label=f'{label} requested')
    ax.semilogy(t, used, lw=0.9, color=color,
                label=f'{label} used')


def main():
    fig, (ax_res, ax_cost, ax_tol) = plt.subplots(
        1, 3, figsize=(9.6, 2.65), sharey=False)

    cases = [case for case in (load_case(case) for case in CASES)
             if case is not None]
    for case in cases:
        plot_residual(ax_res, case, case['time'],
                      f"({case['steps']} steps)")
        plot_residual(ax_cost, case, case['cost'],
                      f"({case['elapsed']:.0f} s)")

    ax_res.axhline(1e-9, lw=0.5, c='k', ls=':',
                   label='target 1e-9')
    ax_res.set_xlabel(r'$t$')
    ax_res.set_ylabel(r'$\|r\|$')
    ax_res.legend(loc='lower left')
    ax_res.grid(True, which='both', lw=0.25, alpha=0.35)
    panel_label(ax_res, 'a')

    ax_cost.axhline(1e-9, lw=0.5, c='k', ls=':',
                    label='target 1e-9')
    ax_cost.set_xlabel('wall time on 8 ranks (s)')
    ax_cost.set_ylabel(r'$\|r\|$')
    ax_cost.legend(loc='lower left')
    ax_cost.grid(True, which='both', lw=0.25, alpha=0.35)
    add_core_min_axis(ax_cost)
    ax_cost.text(0.98, 0.94, r'OIFS wall cost: $0.55\times$ dyn',
                 transform=ax_cost.transAxes, ha='right', va='top',
                 fontsize=6)
    panel_label(ax_cost, 'b')

    maybe_plot_scheduler(ax_tol, Path('sfd_dyn') / 'dyn_tol.dat',
                         'SFD dyn', 'C1')
    maybe_plot_scheduler(ax_tol, Path('sfd_dyn_oifs') / 'dyn_tol.dat',
                         'SFD dyn OIFS', 'C3')
    ax_tol.axhline(1e-9, lw=0.5, c='k', ls=':', label='target 1e-9')
    ax_tol.set_xlabel(r'$t$')
    ax_tol.set_ylabel('solver tolerance')
    ax_tol.legend(loc='best')
    ax_tol.grid(True, which='both', lw=0.25, alpha=0.35)
    panel_label(ax_tol, 'c')

    fname = CASE_DIR / 'residu.png'
    fig.savefig(fname, dpi=500, bbox_inches='tight')
    print(f'Saving {fname.name}')
    plt.close(fig)


if __name__ == '__main__':
    main()
