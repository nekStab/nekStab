#!/usr/bin/env python
#from pymech import neksuite as ns
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import numpy as np
from scipy import sin, cos, pi, exp, tanh, log, interpolate
from scipy.fft import rfftfreq, rfft
import os, os.path, csv, time, pickle
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.signal import savgol_filter


params = {'text.usetex': False,'font.size': 8,'legend.fontsize': 8,'legend.handlelength': 2.5};plt.rcParams.update(params)
transp = False
formt = 'png'
ajust = 'tight'
qual = 400
siz1 = 3.5
siz2 = 2.45
fig_width=3.5
fig_height=2.4

ratios = ['',]
tskps = [2.,]
tmaxs = [0.,]


fld = ''
filename = "flow_rate.dat"

for ratio in range(0,len(ratios)):
    filen = fld + ratios[ratio] + '' + filename
    print('Opening file',filen)
    file = open(filen, 'r')
    data = np.loadtxt(filen).reshape(-1, 2)

    fig, axs = plt.subplots(2, sharex=False)
    fig.set_size_inches(siz1*3, siz2*2)
    axs[0].set_xlabel(r'$t$')
    axs[0].set_ylabel(r'$u$')

    t=data[:,0];u=data[:,1]
    tmin=t.min();tmax=t.max()
    #axs[0].plot(t,u,c='k',ls='-',lw=0.5)

    print('Time series from:',tmin,tmax)

    tmn = float(tskps[ratio])
    tmx = float(tmaxs[ratio])

    if tmx > tmax or tmx < tmax: # sanity check for mininum
        tmx = tmax
        print('Forcing tmx to ',tmx)

    if tmn < tmin or tmn > tmin: # sanity check for maximum
        tmn = tmin
        print('Forcing tmn to ',tmn)

    itmin = np.where(t > tmn)[0][0]
    print(itmin,t[itmin],t[0])

    itmax = np.where(t >= tmx)[0][0]
    print(itmax,t[itmax],t[-1])

    to = t[itmin:itmax]

    dtt = np.zeros(len(to), dtype = np.float64)
    for i in range(1,len(dtt)-1): #compute dt array
        dtt[i] = to[i+1] - to[i]
    
    print(to.shape,dtt.shape)
    dt = dtt.max() #np.round(dtt.max(),decimals=7) # get maximum dt of array to be safe

    print('t old:',to.shape,to.min(),to.max(),' dt=',dt)

    tn = np.linspace(t[itmin], t[itmax], int((t[itmax]-t[itmin])/dt), endpoint=False)
    # endpoint = False so new interval is inside old -- we do not want to extrapolate!
    print('t new:',tn.shape,tn.min(),tn.max())

    print('err (old-new) min max=',to.min()-tn.min(),to.max()-tn.max())

    axs[1].set_xlabel(r'$St$')
    uo = u[itmin:itmax]; nm='u'; cor='k'
    qur = np.mean(uo) # compute mean
    uo -= qur # remove mean

    un = interpolate.interp1d(to,uo,kind='slinear',fill_value="extrapolate")(tn)
    #'zero','slinear','quadratic','cubic'=spline interpolation of 0th,1st,2nd,3rd;

    #qur = np.mean(un) # compute mean
    #un -= qur # remove mean

    axs[0].scatter(tn,un,c='r',ls='-',s=0.7)
    ufft = rfft(un)
    ufreq = rfftfreq(len(un), d=dt)
    ufmax = np.max(ufreq)
    prfreq = ufreq[np.argmax(np.abs(ufft))]
    print('(',nm,') mean, min, max=',round(qur,6),uo.min(),uo.max())
    print('Peak at f=',round(prfreq,6))
    print('        w=',round(2*pi*prfreq,4))
    print('        T=',round(1/prfreq,6))

    axs[1].semilogy(ufreq,np.abs(ufft) ,c=cor,lw=0.7,label=nm+' St='+("{0:.4f}".format(round(prfreq,4))))
    plt.axvline(x=prfreq    ,lw=0.4, color='k', ls=':')
    axs[1].set_xlim(0., 10.0)
    #plt.ylim(1e-5,1e5) 

    plt.legend(loc='best')
    fname = ratios[ratio]+'flow_rate.'+formt
    print('Saving ',fname)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust)

    t = tn
    u = un

    #t = np.linspace(np.min(t), np.max(t), len(p.t), endpoint=False)
    #u = interpolate.interp1d(p.t,p.u,fill_value="extrapolate")(p.t)
    u = savgol_filter(u, 3, 1, deriv=0,mode='nearest')
    ud = savgol_filter(u, 7, 3, deriv=1,mode='nearest')
    udd = savgol_filter(u, 11, 6, deriv=2,mode='nearest')
    dd = 400
    t = t[dd:-dd];u = u[dd:-dd];ud = ud[dd:-dd];udd = udd[dd:-dd]

    x=ud ; y=u ; t=t
    fig, axs = plt.subplots(1, 1);fig.set_size_inches(siz1, siz2)
    plt.xlabel(r'$\dot{v}$');plt.ylabel(r'$v$')
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    norm = plt.Normalize(t.min(), t.max())
    lc = LineCollection(segments, cmap='viridis', norm=norm);lc.set_array(t)
    line = axs.add_collection(lc);fig.colorbar(line,label=r'$t$');
    lc.set_linewidth(0.5);axs.set_xlim(1.1*x.min(), 1.1*x.max());axs.set_ylim(1.1*y.min(), 1.1*y.max())
    fname = 'projection_vd_v_cb.'+formt; print('Saving ',fname)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);plt.close();plt.clf()

    x=udd ; y=u ; t=t
    fig, axs = plt.subplots(1, 1);fig.set_size_inches(siz1, siz2)
    plt.xlabel(r'$\ddot{v}$');plt.ylabel(r'$v$')
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    norm = plt.Normalize(t.min(), t.max())
    lc = LineCollection(segments, cmap='viridis', norm=norm);lc.set_array(t)
    line = axs.add_collection(lc);fig.colorbar(line,label=r'$t$');
    lc.set_linewidth(0.5);axs.set_xlim(1.1*x.min(), 1.1*x.max());axs.set_ylim(1.2*y.min(), 1.1*y.max())
    fname = 'projection_vdd_v_cb.'+formt; print('Saving ',fname)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);plt.close();plt.clf()


    x=udd ; y=ud ; t=t
    fig, axs = plt.subplots(1, 1);fig.set_size_inches(siz1, siz2)
    plt.xlabel(r'$\ddot{v}$');plt.ylabel(r'$\dot{v}$')
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    norm = plt.Normalize(t.min(), t.max())
    lc = LineCollection(segments, cmap='viridis', norm=norm);lc.set_array(t)
    line = axs.add_collection(lc);fig.colorbar(line,label=r'$t$');
    lc.set_linewidth(0.5);axs.set_xlim(1.1*x.min(), 1.1*x.max());axs.set_ylim(1.2*y.min(), 1.1*y.max())
    fname = 'projection_vdd_vd_cb.'+formt; print('Saving ',fname)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);plt.close();plt.clf()

    fig = plt.figure();fig.set_size_inches(siz1*2, siz2*2);ax = fig.gca(projection='3d')
    ax.set_xlabel(r'$v$');ax.set_ylabel(r'$\dot{v}$');ax.set_zlabel(r'$\ddot{v}$')
    #ax.view_init(60, 35);
    factor=1.1
    vmin = np.min(u)*factor;     vmax = np.max(u)*factor
    vdmin = np.min(ud)*factor;   vdmax = np.max(ud)*factor
    vddmin = np.min(udd)*factor; vddmax = np.max(udd)*factor

    ax.plot(u,ud,vddmin, c='r',ls='dotted',lw=0.4)
    ax.plot(vmin*np.ones(len(t)),ud,udd, c='g',ls='dotted',lw=0.4)
    ax.plot(u,vdmax*np.ones(len(t)), udd,c='b',ls='dotted',lw=0.4)
    ax.plot(u,ud,udd,c='k',ls='solid',lw=0.8)
    ax.set_xlim(vmin,vmax);ax.set_ylim(vdmin,vdmax);ax.set_zlim(vddmin,vddmax)
    fname = 'attractor.'+formt;print('Saving fl ',fname)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);plt.close(); plt.clf()
