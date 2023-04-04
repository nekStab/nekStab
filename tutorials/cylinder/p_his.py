#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import numpy as np
from scipy import sin, cos, pi, exp, tanh, log, interpolate
from scipy.fft import rfftfreq, rfft
import os, os.path, csv, time, pickle, glob
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.signal import savgol_filter

params = {'text.usetex': False,
          'font.size': 11,
          'legend.fontsize': 11,
          'legend.handlelength': 2.5,
          'agg.path.chunksize':100000}
plt.rcParams.update(params)

formt = 'png'
ajust = 'tight'
qual = 400
fig_width = 4.3
fig_height = 16*fig_width/9
        
def fft(u,dt):
    ufft = rfft(u)
    ufreq = rfftfreq(len(u), d=dt)
    prfreq = ufreq[np.argmax(np.abs(ufft))]
    print('Peak at f=',round(prfreq,6))
    return ufft,ufreq,prfreq

files = ['1cyl.his']
tskps = [0] # set to 0 to plot all
tmaxs = [0]

for i, filen in enumerate(files):
    print(f'Opening file {filen}')
    
    with open(filen, 'r') as file:
        nps = int(file.readline().rstrip())
        print(f'Number of probes found: {nps}')

        coords = np.zeros((nps, 3))
        for n, line in zip(range(nps), file):
            coords[n] = line.split()

    try:
        data = np.loadtxt(filen, skiprows=nps + 1).reshape(-1, nps, 5)
        is_3d = True
    except:
        data = np.loadtxt(filen, skiprows=nps + 1).reshape(-1, nps, 4)
        is_3d = False
            
    for prb in range(0,nps,1):
        print(f'Probe number {prb + 1}')

        t=data[:,prb,0][1:]
        u=data[:,prb,1][1:]
        v=data[:,prb,2][1:]
        if is_3d:
            w=data[:,prb,3][1:]
            p=data[:,prb,4][1:]
        else:
            p=data[:,prb,3][1:]

        # Adjust time series (crop and interpolate for constant time step)
        tmin, tmax = t.min(), t.max()
        print(f'Original time series: {tmin} {tmax}')

        tmin = max(tmin, float(tskps[i])) if float(tskps[i]) > 0 else tmin
        tmax = min(tmax, float(tmaxs[i])) if float(tmaxs[i]) > 0 else tmax
        print(f'Cropped time series: {tmin} {tmax}')

        itmin, itmax = np.where(t >= tmin)[0][0], np.where(t >= tmax)[0][0]
        to = t[itmin:itmax]
        dtt = np.zeros(len(to), dtype=np.float64)
        for xx in range(1, len(dtt) - 1):
            dtt[xx] = to[xx + 1] - to[xx]

        dt = dtt.max()
        tn = np.linspace(t[itmin], t[itmax], int((t[itmax] - t[itmin]) / dt), endpoint=False)

        # Figure 
        fig, axs = plt.subplots(2, sharex=False)
        fig.set_size_inches(fig_height, fig_width)
        
        # Plot time series (signal on upper part of the figure)
        axs[0].set_xlabel(r'$t$')
        axs[0].set_title(f'x,y,z={coords[prb][0]},{coords[prb][1]},{coords[prb][2]}')
        
        # We focus on the vertical velocity signal
        nm='v' # name of the signal
        cor='r' # color of the signal
        vo = v[itmin:itmax]
        qur = np.mean(vo) # compute mean
        vo -= qur # remove mean
        print('(',nm,') mean, min, max=',round(qur,6),vo.min(),vo.max())
        vn = interpolate.interp1d(to,vo,kind='slinear',fill_value="extrapolate")(tn) #'zero','slinear','quadratic','cubic'=spline interpolation of 0th,1st,2nd,3rd;
        
        axs[0].scatter(tn,vn,c=cor,ls='-',s=0.001) # plot signal

        ufft,ufreq,prfreq = fft(vn,dt) # compute fft

        axs[1].semilogy(ufreq,np.abs(ufft),c=cor,lw=0.6)#,label=nm+' St='+("{0:.4f}".format(round(prfreq,4))))
        axs[1].axvline(x=prfreq         ,c=cor,lw=0.5,ls='--',label='$St=$'+("{0:.4f}".format(round(prfreq,4))))
        axs[1].set_xlabel(r'$St$')
        axs[1].set_xlim(0., 0.4)
        axs[1].set_ylim(1e1,1e5) 

        plt.legend(loc='upper right')
        fname = 'his'+str(prb+1)+'_fft.'+formt
        print('Saving ',fname); print()
        plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust)
        plt.close(); plt.clf()
        
        # Phase space plot
        fig = plt.figure(figsize=(fig_height, fig_width))
        plt.xlabel(r'$v$')
        plt.ylabel(r'$\dot{v}$')
        plt.plot(v[2:-2],np.gradient(v)[2:-2]/dt,c='r',ls='-',lw=0.8)

        fname = 'his'+str(prb+1)+'_phase_space.'+formt
        print('Saving ',fname)
        plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust)