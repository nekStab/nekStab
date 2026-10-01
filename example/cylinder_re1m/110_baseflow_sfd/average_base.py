#!/usr/bin/env python3
"""Time-average the converged-mean SFD snapshots into a steady RANS base flow.

The Re=1e6 k-tau RANS settles to a steady velocity mean with persistent local
wake unsteadiness (SFD residual plateaus ~0.16). Averaging the converged-mean
window (t=100..180) averages out the local fluctuations to give a clean steady
base flow for the 310 FD-Frechet stability stage (mean-flow stability).
Writes BF_1cyl0.f00001.
"""
import copy
import numpy as np
from pymech import readnek, writenek

files = ["1cyl0.f0000%d" % i for i in range(5, 10)]  # t=100,120,140,160,180
flds = [readnek(f) for f in files]
n = len(flds)
print("averaging %d snapshots: %s" % (n, ", ".join(files)))
print("times:", [round(f.time, 2) for f in flds])

avg = copy.deepcopy(flds[0])
for ie in range(len(avg.elem)):
    e = avg.elem[ie]
    e.vel[:]  = sum(f.elem[ie].vel  for f in flds) / n
    e.pres[:] = sum(f.elem[ie].pres for f in flds) / n
    if e.temp.size:
        e.temp[:] = sum(f.elem[ie].temp for f in flds) / n
    if e.scal.size:
        e.scal[:] = sum(f.elem[ie].scal for f in flds) / n

avg.time = float(flds[-1].time)
writenek("BF_1cyl0.f00001", avg)
print("wrote BF_1cyl0.f00001 (time=%.3f)" % avg.time)
