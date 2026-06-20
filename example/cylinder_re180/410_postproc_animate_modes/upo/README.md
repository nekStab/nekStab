# Cylinder Animate Modes with UPO — Floquet mode animation on the periodic orbit

## ⚠ 2D PROXY of a 3D analysis (read first)
The Floquet mode animated here is from a **2D** run (β = 0). Barkley &
Henderson (1996) Mode A/B are **3D** (spanwise-periodic) and are NOT what
this 2D pipeline produces — this is a cost-reduced **proxy/demonstrator**
of the famous case, not a physical secondary-instability animation. Full
explanation: `../../210_baseflow_newton/README.md`.

## Two kinds of mode animation (and why this case only does one)
For a periodic base flow there are, in general, **two** videos one might want:

1. **Animated Floquet mode** — the leading mode evolving over one period.
2. **Mode distorting the base flow** — the base flow with the mode forced on
   top (`U_base + ε·û`), so you see how the instability *deforms* the wake.

**For this 2D cylinder we only produce video (1).** The leading Floquet
multiplier here is **μ = +1 (real, σ ≈ −1.4×10⁻⁵ ≈ 0)** — the trivial
**synchronous / phase mode** of the orbit (≈ ∂U_base/∂t). In the 2D plane
this mode is just the phase shift of the existing shedding: **it cannot
deform the base flow into anything new** — a genuine secondary instability
that distorts the wake (Mode A/B) only exists in **3D** (spanwise
modulation). So video (2) is physically meaningless for this 2D proxy and
is intentionally omitted. In a real 3D Floquet study you would also render
video (2).

## What "animated μ=+1 mode" means here
> **The delivered video shows the BASE FLOW, not the eigenvector field.**
> It is the evolving von Kármán shedding (`1cyl0.f*`), used as the physical
> stand-in for the μ=+1 *synchronous phase mode*. Read on for why that is
> the correct visualization — and what the alternative (the genuine
> eigenfunction) would be.

A Floquet mode is `u'(x,t) = e^{σt} û(x,t)` with `û` T-periodic. For a μ=+1
mode, σ ≈ 0 and **ω = 0**. There are three distinct objects, easy to
confuse:

| Object | File(s) | What it is | Animates? |
|---|---|---|---|
| **Base flow** (what this video uses) | `1cyl0.f*` | the periodic orbit over one period T | **yes** — it's the shedding |
| **Eigenvector reconstruction** | `dQm…` = `Re·cos(ωt) − Im·sin(ωt)` | the stored mode pair rotated in the complex plane | **no** — *frozen* for ω=0 |
| **Genuine eigenfunction** `û(x,t)` | (not produced here) | the leading mode evolved along the orbit by the *linearized* solver | yes — and it's the actual mode |

Because ω = 0, the cheap reconstruction `Re·cos − Im·sin` collapses to a
**frozen field** — it cannot animate. The μ=+1 mode is the trivial phase
mode (≈ ∂U_base/∂t): its entire time-dependence lives in the **base-flow
orbit it rides on**, not in any envelope. So the honest "animation of the
μ=+1 mode" is the **base-flow DNS advanced one period T** — the von Kármán
street shedding and advecting, at the base-flow frequency `ω_base = 2π/T`.
That is exactly what this case renders.

**If you want the genuine eigenfunction `û(x,t)` instead** (the mode's *own*
spatial structure changing through the period, distinct from the base flow),
that is a different run: seed the **linearized perturbation** solver with the
leading eigenvector `dRe`/`dIm` and integrate one period along the orbit,
outputting the perturbation field. For a μ=+1 mode this would re-trace the
phase mode; for an oscillatory mode (μ complex, ω ≠ 0) it is the standard way
to animate the eigenfunction. This case does **not** do that — it uses the
base-flow proxy, which is sufficient and physically faithful for μ=+1.

## nekStab Mode
`userParam01 = 4.52` — `animate_mode_Floquet` (src/sensitivity.f90):
advances the base-flow DNS one period from the on-orbit UPO state and writes
a frame sequence of the evolving flow.

- `userParam07 = 100` — number of output frames over one period.
- `targetCFL = 0.5` — auto `dt` from CFL (variable-dt selection inside the
  routine; here dt ≈ 7.4×10⁻³, 700 steps, 100 frames).
- The routine reads `(σ, ω)` from row 1 of `Spectre_NSd_conv.dat` and sets
  the animation duration to **one period = 2π/ω**. Here ω is set to the
  **base-flow frequency** `ω_base = 2π/T = 1.21049` (T = 5.1906, the
  converged UPO period from `../../210_baseflow_newton/`) so the DNS
  runs exactly one shedding period. (The *physical* direct-Floquet spectrum
  in `../../311_stability_direct_floquet/` keeps the true μ=+1 row with
  ω = 0; this case's copy overrides ω only to set the animation period.)

## Prerequisites (copy from upstream stages)
- `BF_1cyl0.f00001` — on-orbit UPO state (from `../../210_baseflow_newton/`).
- `dRe1cyl0.f00001`, `dIm1cyl0.f00001` — leading Floquet eigenvector pair
  (from `../../311_stability_direct_floquet/`). Loaded by the routine even
  though the μ=+1 reconstruction is static; only the evolving base flow is
  animated.
- `Spectre_NSd_conv.dat` with row 1 ω = `ω_base` (see above).
- `SESSION.NAME` line 2 must be this directory's absolute path (trailing `/`).

## Run
```bash
makeneks 1cyl                       # compile (links this case's .usr)
mpirun -np 8 ./nek5000 > logfile    # advances base flow one period T
```
Output: `1cyl0.f00001 … f00100` — 100 snapshots of the evolving base flow
over one period (vorticity-rate Ω_R carried in the temperature slot).

## Render the MP4
```bash
# from example/  (ffmpeg bundled by imageio-ffmpeg; SEM→Cartesian via pymech)
uv run --with pymech --with imageio --with imageio-ffmpeg --with scipy --with matplotlib \
  python render_mode_video.py \
    'cylinder_re180/410_postproc_animate_modes/upo/1cyl0.f*' \
    cylinder_re180/410_postproc_animate_modes/upo/floquet_mode_animation.mp4 \
    --fps 20 --title "Re=180 UPO  μ=+1 synchronous mode (base-flow shedding ω_z)" \
    --circle 0 0 0.5
```
`render_mode_video.py` is **case-generic**: it interpolates each SEM frame
onto a Cartesian grid, computes spanwise vorticity `ω_z = ∂v/∂x − ∂u/∂y`,
and writes an H.264 MP4. Tune `--bbox xmin xmax ymin ymax`, `--nx/--ny`,
`--title`, and the optional solid-body overlay `--circle X Y R` (omit it for
non-cylinder geometries) per case.

## Expected Output
`floquet_mode_animation.mp4` — one full period of the Re=180 von Kármán
street: vortices born in the near wake and advecting downstream, the
pattern convecting ~one wavelength per period. This is the μ=+1 synchronous
Floquet mode at the base-flow frequency.

## Reference
Barkley & Henderson, *J. Fluid Mech.* **322** (1996) 215 — Floquet stability
of the periodic cylinder wake (Mode A/B are the 3D instabilities this 2D
case proxies).
