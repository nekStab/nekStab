#!/usr/bin/env python3
"""Visualize the SFD dynamic + OIFS testcase.

This figure must compare against both non-OIFS baselines.  OIFS reduces
the number of timesteps substantially, but it is a case-specific acceleration;
the residual and measured-cost panels make the speed/accuracy tradeoff explicit.
"""
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import nekplot as nk

CASE_DIR = Path(__file__).resolve().parent
BASEFLOW_DIR = CASE_DIR.parent
OUTPUT = CASE_DIR / 'plot.png'
N_RANKS = 8

# Latest verified 8-rank timings from the matching Re=50 runs.  The OIFS run
# takes many fewer timesteps, but each OIFS step is more expensive than a
# standard BDF step, so a timestep count alone exaggerates the practical gain.
CASES = [
    {
        'rel_path': Path('sfd') / 'residu.dat',
        'label': 'SFD fixed',
        'color': '0.35',
        'linestyle': '--',
        'elapsed': 982.240,
        'zorder': 1,
    },
    {
        'rel_path': Path('sfd_dyn') / 'residu.dat',
        'label': 'SFD dyn',
        'color': 'C1',
        'linestyle': '-',
        'elapsed': 810.623,
        'zorder': 2,
    },
    {
        'rel_path': Path('sfd_dyn_oifs') / 'residu.dat',
        'label': 'SFD dyn OIFS',
        'color': 'C3',
        'linestyle': '-',
        'elapsed': 617.122,
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
        entry['cost'] = wall_cost_axis(len(t), case['elapsed'])
        entry['steps'] = len(t)
        cases.append(entry)
    return cases


def add_core_min_axis(ax):
    sec_to_core_min = lambda seconds: seconds * N_RANKS / 60.0
    core_min_to_sec = lambda core_min: core_min * 60.0 / N_RANKS
    top = ax.secondary_xaxis('top', functions=(sec_to_core_min,
                                               core_min_to_sec))
    top.set_xlabel('core-min')
    top.tick_params(labelsize=5.5, pad=1)


def panel_label(ax, label):
    ax.text(0.03, 0.96, f'({label})', transform=ax.transAxes,
            ha='left', va='top', fontsize=8, fontweight='bold',
            bbox=dict(facecolor='white', edgecolor='none',
                      alpha=0.78, pad=1.1))


def find_result_field():
    # Prefer the actual end-time checkpoint.  BFRe40_1cyl0.f00001 is the
    # shared initial condition and must not be plotted as the testcase result.
    for pattern in ('1cyl0.f*', 'BF_1cyl0.f*', 'BF*1cyl0.f*'):
        matches = nk.find_fields(pattern, CASE_DIR)
        if matches:
            return matches[0]
    return None


def main():
    nk.configure_style()
    fig = plt.figure(figsize=(nk.COL_WIDTH * 2.35, nk.COL_WIDTH * 0.55))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.8, 0.85, 1.0], wspace=0.38)
    ax_time = fig.add_subplot(gs[0])
    ax_cost = fig.add_subplot(gs[1], sharey=ax_time)
    ax_bf = fig.add_subplot(gs[2])

    cases = load_available_cases()

    for case in cases:
        plot_residual(ax_time, case, case['time'],
                      label_suffix=f" ({case['steps']} steps)")
        plot_residual(ax_cost, case, case['cost'],
                      label_suffix=f" ({case['elapsed']:.0f} s)")

    for ax in (ax_time, ax_cost):
        ax.axhline(1e-9, lw=0.5, c='k', ls=':', label='1e-9 target')
        ax.grid(True, which='both', lw=0.25, alpha=0.35)

    ax_time.set_xlabel(r'$t$')
    ax_time.set_ylabel(r'$\|r\|$')
    ax_time.legend(fontsize=4.8, loc='lower left',
                   handlelength=1.2, borderpad=0.25)
    panel_label(ax_time, 'a')

    ax_cost.set_xlabel('wall time on 8 ranks (s)')
    ax_cost.tick_params(labelleft=False)
    add_core_min_axis(ax_cost)
    ax_cost.text(0.98, 0.94, r'OIFS wall cost: $0.76\times$ dyn',
                 transform=ax_cost.transAxes, ha='right', va='top',
                 fontsize=5.2)
    panel_label(ax_cost, 'b')

    field_file = find_result_field()
    if field_file is not None:
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
        panel_label(ax_bf, 'c')
    else:
        ax_bf.axis('off')
        ax_bf.text(0.5, 0.5, 'No field file found',
                   ha='center', va='center', transform=ax_bf.transAxes)

    fig.savefig(OUTPUT, dpi=600, bbox_inches='tight')
    print(f'Saved {OUTPUT}')
    plt.close(fig)


if __name__ == '__main__':
    main()
