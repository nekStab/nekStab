#!/usr/bin/env python3
"""Extract the Re=175 shedding period from the wake-probe time series (1cyl.his).
PSD peak for f_peak + sharpness; cycle-counting over the long window for a
sub-bin-precise period to seed the Re=180 Newton UPO."""
import numpy as np

# --- load .his: line1=npts, line2=coords, then rows: t, u, v, (w/p) ---
raw = np.loadtxt('1cyl.his', skiprows=2)
t, u, v = raw[:, 0], raw[:, 1], raw[:, 2]

# drop the Re180->Re175 relaxation transient (first 300 t.u.)
t0 = t[0] + 300.0
m = t >= t0
t, v = t[m], v[m]
print(f"window: t=[{t[0]:.1f}, {t[-1]:.1f}]  span={t[-1]-t[0]:.1f} t.u.  N={len(t)}")

# --- uniform resample (variableDt -> non-uniform) ---
dt = np.median(np.diff(t))
tu = np.arange(t[0], t[-1], dt)
vu = np.interp(tu, t, v)
vu = vu - vu.mean()

# --- PSD (FFT, Hann window, zero-pad x4 for finer peak) ---
w = np.hanning(len(vu))
nfft = 1 << int(np.ceil(np.log2(len(vu) * 4)))
V = np.fft.rfft(vu * w, n=nfft)
f = np.fft.rfftfreq(nfft, d=dt)
P = np.abs(V) ** 2
k = np.argmax(P[1:]) + 1
fpk = f[k]
# parabolic sub-bin refine
a, b, c = P[k-1], P[k], P[k+1]
dk = 0.5 * (a - c) / (a - 2*b + c) if (a - 2*b + c) != 0 else 0.0
fpk_ref = (k + dk) * (f[1] - f[0])
# sharpness: -3dB width in bins
half = P[k] / 2
lo = k; hi = k
while lo > 0 and P[lo] > half: lo -= 1
while hi < len(P)-1 and P[hi] > half: hi += 1
q = fpk / ((hi - lo) * (f[1]-f[0]) + 1e-30)
print(f"PSD peak: f={fpk_ref:.6f}  St={fpk_ref:.4f}  T_fft={1/fpk_ref:.5f}  Q~{q:.0f}  (binDf={f[1]-f[0]:.2e})")

# --- cycle-counting for sub-bin precision: upward zero-crossings ---
s = np.signbit(vu)
up = np.where((~s[1:]) & (s[:-1]))[0]   # - -> + crossings
# linear-interp crossing times
tc = tu[up] - vu[up] * (tu[up+1]-tu[up]) / (vu[up+1]-vu[up])
ncyc = len(tc) - 1
T_cyc = (tc[-1] - tc[0]) / ncyc
print(f"cycle-count: {ncyc} periods over {tc[-1]-tc[0]:.2f} t.u.  ->  T={T_cyc:.5f}  f={1/T_cyc:.6f}")
print(f"\n==> PERIOD for Re=180 UPO seed: T = {T_cyc:.5f}  (St={1/T_cyc:.4f})")
