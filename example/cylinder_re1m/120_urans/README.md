# cylinder_re1m / 120_urans

**Stage**: bare forward k-tau march (`userParam01 = 0`).
**Re**: 1e6, 2D k-tau, mesh `lelg=1480`.

This operator does not shed. Four marches:

- Release `base_converged.f00001` at `dt=0.001`, through t=70: mid-wake v stays near 1e-3, u settles at 0.63.
- Same release at `dt=0.005`, through t=132: same steady wake.
- Cold start, `dt=0.005`, to t=40: mid-wake v stays near 1e-6.
- A 0.05 Gaussian v-seed at (x,y)=(5,0) on the settled wake decays by t=40, both with hpfrt (`|v|` max 0.003) and with `filtering = none` (`|v|` max 0.001).

A June run with v about ±0.87, dated 2026-06-24, was before the 2D mesh alignment. It is not this operator, and that figure is not kept. The current binary still contains that one-time seed; do not treat a fresh run as unperturbed until it is rebuilt without it.

The unforced march leaves the SFD field and sits on a steady wake. Mid-wake u is about 0.61, not the SFD value. `steady_wake.f00001` is the t=40 field after the seed decayed onto that wake (the cold-start t=40 file was overwritten by the seeded runs). There is no resolved shedding stress to average. The large growth rates are the spectrum of the imperfect SFD base, not the spectrum of this steady wake, and not a mean-flow spectrum. This march is not the mean-flow example: a mean flow is the average of a fluctuating saturated run, and R_ij is that fluctuation correlation.

## Reynolds-stress force

This case is the write example. At the last step `1cyl.usr` calls `reynolds_commit`. That writes `RS1` (`R_uu`, `R_vv`, `R_ww`), `RS2` (`R_uv`, `R_vw`, `R_wu`), and `FRS1cyl0.f00001` (`f = -div(R)`). It does not turn the force on.

The run that should hold that mean calls this from `nekStab_usrchk`, and nothing else:

```fortran
subroutine nekStab_usrchk
   use nekstab_reynolds, only: reynolds_enable
   call reynolds_enable
end subroutine
```

`nekStab_init` then loads `FRS1cyl0.f00001` once and keeps it. `310_stability_direct` does not call `reynolds_enable`. That stage is the SFD base, not this mean.
