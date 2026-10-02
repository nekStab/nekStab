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
