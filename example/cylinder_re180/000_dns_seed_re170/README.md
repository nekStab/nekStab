# Re=170 DNS — on-orbit seed + period source for the Re=180 UPO

Purpose: produce (a) a saturated **Re=170** limit-cycle field to seed the
Re=180 Newton-UPO (`../210_baseflow_newton/upo/`), and (b) the **shedding
period T** by FFT of a wake probe, used as the UPO `endTime`.

Why Re=170 (not 150, not 180): close enough to Re=180 that the UPO's
period-Newton converges (small initial residual → controlled ΔT), yet
safely below the 3D Mode-A onset (Re_A≈188) so the flow stays 2D. The
older Re=150 seed is too far off-orbit and makes the period-Newton
diverge (T explodes 5.5→538→…). See the UPO README for the full failure.

## How it was run

- Mesh/SIZE/usr copied from `../000_dns_seed/`; `viscosity = -170.0`.
- Seeded from a saturated Re=150 frame (`ic_re150.f00001`); re-saturates
  to the Re=170 orbit in ~50 t_c.
- `endTime = 1100` (continued from t=600) → ~100 periods of clean signal.
- Wake probe via `hpts`: `1cyl.his` header is
  ```
  1
  5 0 0
  ```
  (one point at x=5, y=0 downstream). Nek appends `time u v p` per step.

## Extracting the period

FFT the transverse velocity `v` (column 3 of `1cyl.his`) over the
saturated window (drop the re-saturation transient), parabolic-interpolate
the peak — see the snippet in `../210_baseflow_newton/upo/README.md`.

Result: **St = 0.1904, T = 5.251** (86-period FFT window).

## Outputs consumed downstream

- Latest saturated snapshot `1cyl0.f0000N` → copy to
  `../210_baseflow_newton/upo/rstcyl0.f00001` (the UPO seed).
- `T = 5.251` → set as `endTime` in the UPO `.par`.
