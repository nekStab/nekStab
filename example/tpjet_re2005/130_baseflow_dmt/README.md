# tpjet DMT Baseflow Validation

Validation case for Dynamic Mode Tracking on the tpjet baseflow at Re=1900.

The case starts from `BF_tpjet0.f00001`, targets `St_D = 0.6` with
`omc = 3.7699`, and runs DMT mode through `userParam01 = 1.3`.

Mesh, SIZE, and user routines are linked from `../baseflow/newton`; the Slurm
script writes `SESSION.NAME` dynamically so Nek reads this directory's `case.par`.
