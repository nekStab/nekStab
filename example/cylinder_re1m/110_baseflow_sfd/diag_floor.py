#!/usr/bin/env python3
"""What is the 5e-4 residual floor made of? Decompose by field + region.

(1) Per-field L2 change between the last two SFD snapshots (f00004->f00005,
    consecutive in time near t=100). Tells us whether velocity (physical) or
    the k/tau scalars (limiter noise) dominate the residual floor.
(2) v statistics split near-body (x<2) vs wake (x>3): a steady base has wake
    v ~ 0; residual shedding would light up the wake.
"""
import numpy as np
from pymech import readnek

def fields(path):
    f = readnek(path)
    g = lambda sel: np.concatenate([sel(e).ravel() for e in f.elem])
    return dict(
        t=f.time,
        x=g(lambda e: e.pos[0]), y=g(lambda e: e.pos[1]),
        u=g(lambda e: e.vel[0]), v=g(lambda e: e.vel[1]),
        p=g(lambda e: e.pres[0]),
        k=g(lambda e: e.scal[0]), tau=g(lambda e: e.scal[1]),
    )

a = fields("1cyl0.f00004")
b = fields("1cyl0.f00005")
dt_snap = b["t"] - a["t"]
n = len(a["u"])
print(f"snapshot times: {a['t']:.4f} -> {b['t']:.4f}  (dt_snap={dt_snap:.4f})\n")

print("per-field L2 change ||b-a||/sqrt(N)  and  relative to field rms:")
for key in ["u", "v", "p", "k", "tau"]:
    d = b[key] - a[key]
    l2 = np.sqrt(np.mean(d**2))
    rms = np.sqrt(np.mean(b[key]**2)) or 1.0
    print(f"  {key:4s}: abs={l2:.3e}   rel={l2/rms:.3e}   (field rms={rms:.3e})")

print("\nv split by region (steady base => wake v ~ 0):")
x, v = b["x"], b["v"]
nb = (x > -2) & (x < 2)
wk = (x > 3) & (x < 25) & (np.abs(b["y"]) < 3)
print(f"  near-body |x|<2 : |v|max={np.abs(v[nb]).max():.3e}  rms={np.sqrt(np.mean(v[nb]**2)):.3e}")
print(f"  wake x in [3,25]: |v|max={np.abs(v[wk]).max():.3e}  rms={np.sqrt(np.mean(v[wk]**2)):.3e}")
print(f"  wake mean v     : {v[wk].mean():+.3e}  (asymmetry; 0 = top-bottom symmetric)")

# negative k extent (limiter undershoot)
k = b["k"]
print(f"\nnegative-k cells: {(k<0).sum()}/{n} ({100*(k<0).mean():.2f}%)  min k={k.min():.3e}")
