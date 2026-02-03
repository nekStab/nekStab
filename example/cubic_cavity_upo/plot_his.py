#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import rfft, rfftfreq
from scipy.interpolate import interp1d
import os
import glob

PARAMS = {
    'text.usetex': False,
    'font.size': 11,
    'legend.fontsize': 11,
    'legend.handlelength': 2.5,
    'agg.path.chunksize': 100000
}
plt.rcParams.update(PARAMS)

format = 'png'
adjust = 'tight'
quality = 400
fig_width = 4.3
fig_height = 16 * fig_width / 9

def fft_analysis(x, dt):
    N = len(x)
    X = rfft(x)
    freqs = rfftfreq(N, dt)
    psd = np.abs(X)**2 * dt / N
    peak_freq = freqs[np.argmax(psd)]
    print(f'Peak at f = {peak_freq:.6f} with T = {1.0/peak_freq}')
    return X, freqs, psd, peak_freq

def process_file(filename, tskip=0, tmax=0):
    print(f'Opening file {filename}')
    
    with open(filename, 'r') as file:
        nps = int(file.readline().strip())
        print(f'Number of probes found: {nps}')
        coords = np.array([next(file).split() for _ in range(nps)], dtype=float)
    
    data = np.loadtxt(filename, skiprows=nps + 1)
    is_3d = data.shape[1] == 5 * nps
    data = data.reshape(-1, nps, 5 if is_3d else 4)
    
    for prb in range(nps):
        process_probe(data, coords, prb, is_3d, tskip, tmax)

def process_probe(data, coords, prb, is_3d, tskip, tmax):
    print(f'Processing probe number {prb + 1}')
    
    t = data[:, prb, 0]
    u = data[:, prb, 1]
    v = data[:, prb, 2]
    w = data[:, prb, 3] if is_3d else None
    p = data[:, prb, 4] if is_3d else data[:, prb, 3]
    
    t, u, v, w = crop_and_interpolate(t, u, v, w, tskip, tmax)
    
    plot_time_series_and_fft(t, v, coords[prb], prb)
    plot_phase_spaces(u, v, w, prb, is_3d)

def crop_and_interpolate(t, u, v, w, tmin, tmax):
    tmin = max(t.min(), tmin) if tmin > 0 else t.min()
    tmax = min(t.max(), tmax) if tmax > 0 else t.max()
    print(f'Time series range: {tmin:.6f} to {tmax:.6f}')
    
    mask = (t >= tmin) & (t <= tmax)
    t, u, v = t[mask], u[mask], v[mask]
    if w is not None:
        w = w[mask]
    
    dt_max = np.diff(t).max()
    tn = np.arange(t[0], t[-1], dt_max)
    un = interp1d(t, u, kind='slinear')(tn)
    vn = interp1d(t, v, kind='slinear')(tn)
    wn = interp1d(t, w, kind='slinear')(tn) if w is not None else None
    
    return tn, un, vn, wn

def plot_time_series_and_fft(t, v, coords, prb):
    fig, (ax1, ax2) = plt.subplots(2, figsize=(fig_height, fig_width), sharex=False)
    
    v_mean = np.mean(v)
    v -= v_mean
    print(f'(v) mean, min, max = {v_mean:.6f}, {v.min():.6f}, {v.max():.6f}')
    
    ax1.set_xlabel('t')
    ax1.set_title(f'x,y,z = {coords[0]:.2f}, {coords[1]:.2f}, {coords[2]:.2f}')
    ax1.scatter(t, v, c='r', s=0.001)
    
    dt = t[1] - t[0]
    X, freqs, psd, peak_freq = fft_analysis(v, dt)
    
    ax2.semilogy(freqs, psd, c='r', lw=0.6)
    ax2.axvline(x=peak_freq, c='r', lw=0.5, ls='--', label=f'St = {peak_freq:.4f}')
    ax2.set_xlabel('St')
    ax2.set_xlim(0.01, 0.4)
    ax2.set_ylim(1e-10, None)
    ax2.legend(loc='upper right')
    
    plt.tight_layout()
    filename = f'his{prb+1}_fft.{format}'
    plt.savefig(filename, format=format, dpi=quality, bbox_inches=adjust)
    plt.close()
    print(f'Saved {filename}')

def plot_phase_spaces(u, v, w, prb, is_3d):
    if is_3d:
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(fig_height, fig_width/1.9))
        
        ax1.set_xlabel('u')
        ax1.set_ylabel('v')
        ax1.plot(u, v, c='r', ls='-', lw=0.8)
        
        ax2.set_xlabel('u')
        ax2.set_ylabel('w')
        ax2.plot(u, w, c='g', ls='-', lw=0.8)
        
        ax3.set_xlabel('v')
        ax3.set_ylabel('w')
        ax3.plot(v, w, c='b', ls='-', lw=0.8)
        
        plt.tight_layout()
        filename = f'his{prb+1}_phase_space_3d.{format}'
    else:
        fig, ax = plt.subplots(figsize=(fig_height, fig_width))
        
        ax.set_xlabel('u')
        ax.set_ylabel('v')
        ax.plot(u, v, c='r', ls='-', lw=0.8)
        
        plt.tight_layout()
        filename = f'his{prb+1}_phase_space_2d.{format}'
    
    plt.savefig(filename, format=format, dpi=quality, bbox_inches=adjust)
    plt.close()
    print(f'Saved {filename}')

if __name__ == "__main__":
    
    #files = ['1cyl.his']
    files = glob.glob('*.his')[0]
    
    print(f"Processing files: {files}")
    
    tskip = 300  # You can modify this value if needed
    tmax = 0   # You can modify this value if needed
    
    process_file(files, tskip, tmax)
