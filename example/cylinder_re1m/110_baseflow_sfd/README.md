# cylinder_re1m / 110_baseflow_sfd

**Stage**: SFD base flow (`userParam01=1.1`, Casacuberta when `userParam04 < 0`).
**Re**: 1e6, 2D k-tau RANS, same mesh as `000_dns/` and `310_stability_direct/` (`lelg=1480`).

The steady k-tau solution is unstable. A bare march sheds. SFD holds the march on that state so `310_stability_direct/coupled/` and `quasilaminar/` have a mean to linearise. The restart field those stages read, `base_converged.f00001`, is a local Nek field (`*.f?????` is gitignored).

`userParam04 = -0.192` is the measured wake Strouhal number used as the SFD cutoff. `dt = 0.001` is fixed. The k-tau source is stiff at larger steps.
