#!/usr/bin/env python3
"""Write the initial field of the DNS check: orbit snapshot + eps * Floquet mode.

Usage: make_ic.py <orbit snapshot> <mode field> <eps> <output>

The velocity of the output is the orbit velocity plus eps times the velocity of
the mode. The mode is the real part of the leading Floquet eigenvector. Its
field has unit norm. Pressure and time stamp are those of the orbit snapshot.
Fields are read and written with pymech, at the precision of the orbit snapshot.
"""
import sys

import numpy as np
from pymech.neksuite import readnek, writenek

orbit, mode, eps, out = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]
base = readnek(orbit)
pert = readnek(mode)
rms_base = np.sqrt(np.mean([np.sum(e.vel ** 2, axis=0).mean() for e in base.elem]))
rms_pert = np.sqrt(np.mean([np.sum(e.vel ** 2, axis=0).mean() for e in pert.elem]))
for eb, ep in zip(base.elem, pert.elem):
    eb.vel = eb.vel + eps * ep.vel
writenek(out, base)
print(f'rms(orbit) = {rms_base:.4e}  rms(mode) = {rms_pert:.4e}  eps = {eps:g}'
      f'  relative perturbation = {eps * rms_pert / rms_base:.2e}')
