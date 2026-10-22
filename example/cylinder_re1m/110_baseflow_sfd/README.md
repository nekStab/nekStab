# cylinder_re1m / 110_baseflow_sfd

**Stage**: SFD base flow (`userParam01=1.1`, Casacuberta when `userParam04 < 0`).
**Re**: 1e6, 2D k-tau RANS, same mesh as `000_dns/` and `310_stability_direct/` (`lelg=1480`).
The steady k-tau solution is unstable under the Jacobian, but the bare forward
march on this mesh does not shed -- see `../120_urans/`. SFD still provides a
controlled way to hold the march on the SFD branch so `310_stability_direct/
coupled/` and `quasilaminar/` have a base to linearise. The restart field
those stages read, `base_converged.f00001`, is a local Nek field
(`*.f?????` is gitignored). This field is a base flow, not a mean flow, and
the Reynolds stress of an SFD run is the residual fluctuation correlation,
not k, tau, or eddy viscosity.

`userParam04 = -0.192` is the measured wake Strouhal number used as the SFD
cutoff (from the earlier shedding run, not this 2D mesh). `dt = 0.001` is
fixed; the k-tau source is stiff at larger steps.

## Newton probe from this field (not a stage)
Newton-GMRES (`userParam01 = 2`, `userParam07 = 5`) was started from
`base_converged.f00001` with the same mesh, dt and k-tau model. The first
forward residual was 0.54. Five Arnoldi columns left it there, at rate 0.9999.
The linear residual then jumped to about 100 and stayed there. The next
Newton residual was 0.542, and the field written by Newton matches the SFD
velocity exactly. Newton therefore does not reach a fixed point of the coupled
k-tau map from this field, and the SFD field stays the base for the stability
stages. The stage was removed from the release because it ends in this failure.
`newton_vs_sfd.png` compares the Newton output with the SFD field.
