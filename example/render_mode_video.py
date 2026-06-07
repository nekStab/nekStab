#!/usr/bin/env python
"""Render a sequence of Nek5000 field frames into an MP4 (spanwise vorticity).

Pedagogical helper for the Floquet mode-animation demo
(cylinder_re180/410_postproc_animate_modes/upo): turns the per-frame field
files written by the `animate modes` mode (userParam01=4.51) into a movie.

Pipeline: Nek SEM frames -> Cartesian grid (reused Delaunay, via
nek_to_pymodal) -> spanwise vorticity w_z = dv/dx - du/dy -> matplotlib
frames -> MP4 (ffmpeg bundled by imageio-ffmpeg; no system ffmpeg needed).

Usage:
    python render_mode_video.py 'CASE/dQm1cyl0.f*' out.mp4 [--fps 20]
"""
import argparse
import glob
import sys
import tempfile

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])  # find nek_to_pymodal
from nek_to_pymodal import nek_snapshots_to_npz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("frames", help="glob of Nek frame files (quote it)")
    ap.add_argument("out", help="output .mp4")
    ap.add_argument("--fps", type=int, default=20)
    ap.add_argument("--bbox", type=float, nargs=4,
                    default=(-2.0, 12.0, -3.0, 3.0), help="xmin xmax ymin ymax")
    ap.add_argument("--nx", type=int, default=400)
    ap.add_argument("--ny", type=int, default=180)
    ap.add_argument("--title", default="Floquet mode ω_z",
                    help="frame title prefix")
    ap.add_argument("--circle", type=float, nargs=3, default=None,
                    metavar=("X", "Y", "R"),
                    help="overlay a solid body (e.g. cylinder: 0 0 0.5)")
    args = ap.parse_args()

    paths = sorted(glob.glob(args.frames))
    if not paths:
        sys.exit(f"no frames match {args.frames!r}")
    print(f"{len(paths)} frames")

    # SEM -> Cartesian (u, v on a regular grid); reuses one triangulation
    with tempfile.NamedTemporaryFile(suffix=".npz") as tmp:
        nek_snapshots_to_npz(paths, tmp.name, bbox=tuple(args.bbox),
                             nx=args.nx, ny=args.ny, fields=("vx", "vy"))
        d = np.load(tmp.name)
        x, y, u, v = d["x"], d["y"], d["u"], d["v"]

    # spanwise vorticity per frame: w_z = dv/dx - du/dy
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    wz = np.empty_like(u)
    for k in range(u.shape[0]):
        dvdx = np.gradient(v[k], dx, axis=0)
        dudy = np.gradient(u[k], dy, axis=1)
        wz[k] = dvdx - dudy

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import imageio.v2 as imageio
    import imageio_ffmpeg  # noqa: F401  (registers the bundled ffmpeg)

    lim = np.percentile(np.abs(wz), 99)
    levels = np.linspace(-lim, lim, 31)
    Xg, Yg = np.meshgrid(x, y, indexing="ij")

    writer = imageio.get_writer(args.out, fps=args.fps, codec="libx264",
                                quality=8, macro_block_size=None)
    for k in range(wz.shape[0]):
        fig, ax = plt.subplots(figsize=(7, 3.2), dpi=130)
        ax.contourf(Xg, Yg, wz[k], levels=levels, cmap="RdBu_r", extend="both")
        if args.circle is not None:
            cx, cy, cr = args.circle
            ax.add_patch(plt.Circle((cx, cy), cr, color="k", zorder=5))
        ax.set_aspect("equal"); ax.set_xlim(args.bbox[0], args.bbox[1])
        ax.set_ylim(args.bbox[2], args.bbox[3]); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{args.title}   frame {k+1}/{wz.shape[0]}", fontsize=9)
        fig.tight_layout(pad=0.2)
        fig.canvas.draw()
        frame = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
        frame = frame.reshape(fig.canvas.get_width_height()[::-1] + (4,))[..., :3]
        writer.append_data(frame)
        plt.close(fig)
    writer.close()
    print(f"wrote {args.out}  ({wz.shape[0]} frames @ {args.fps} fps)")


if __name__ == "__main__":
    main()
