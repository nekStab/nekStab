#!/usr/bin/env python3
"""plot.py: Visualize cubic-cavity direct Floquet results (Re=1950, 3D).

Floquet multiplier spectrum on the unit circle (geometry-independent).
The leading-mode field is 3D; a volume/slice rendering is a TODO (the shared
nekplot helpers are 2D-tricontour, so the mode panel is intentionally omitted
here rather than rendered with a wrong 2D projection).

OUTPUTS: plot_spectrum.png
USAGE:   PYTHONPATH=<repo>/example python plot.py
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import nekplot as nk
import matplotlib.pyplot as plt

CASE_DIR = Path(__file__).resolve().parent


def main():
    nk.configure_style()

    spec_h = CASE_DIR / 'Spectre_Hd.dat'
    if spec_h.exists():
        fig, ax = plt.subplots(1, 1, figsize=(nk.COL_WIDTH * 0.5,
                                              nk.COL_WIDTH * 0.5))
        nk.setup_unit_circle_axes(ax, lim=1.5)
        nk.plot_spectrum_H_paper(ax, spec_h, label=r'$Re = 1950$')
        out = CASE_DIR / 'plot_spectrum.png'
        fig.savefig(out, dpi=600, bbox_inches='tight')
        print(f'Saved {out}')
        plt.close(fig)
    else:
        print('Spectre_Hd.dat not found; nothing to plot')


if __name__ == '__main__':
    main()
