# cylinder_re1m / 310_stability_direct/findiff

**Stage**: direct stability via finite-difference Frechet (RANS; `iffindiff=.true.` in .usr).
**Re**: 1e6 (user-directed).
**Status**: files in place; BLOCKED on a RANS base flow (startFrom). Set up SFD/RANS baseflow first, then wire startFrom=BF before running.
**Note**: confirm the 1cyl.re2 mesh resolves the Re=1e6 RANS boundary layer before the production run.
