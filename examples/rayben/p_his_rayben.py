#!/usr/bin/env python
from pymech import neksuite as ns
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import numpy as np
import scipy as sp
from scipy import sin, cos, pi, exp, tanh, log
from scipy import interpolate
from scipy.fftpack import fftfreq, fft
import os, os.path, csv, time, pickle

params = {'text.usetex': False,'font.size': 8,'legend.fontsize': 8,'legend.handlelength': 2.5};plt.rcParams.update(params)
transp = False
formt = 'png'
ajust = 'tight'
qual = 400
siz1 = 3.5
siz2 = 2.45

def fft_probe(prb,fld,tsk):
    pskip = 1
    print()
    print('FFT Probe',prb,'at x,y,z='+str(x[prb-1])+','+str(y[prb-1])+','+str(z[prb-1]))
    prb-=1
    t=data[:,prb,0]
    u=data[:,prb,1]
    v=data[:,prb,2]
    w=data[:,prb,3]
    p=data[:,prb,4]
    tm=data[:,prb,5]
    if tsk == 0:
        tsk =  t[0]
    itsk = np.where(t > tsk)[0][0]
    to = t[itsk:-1]# - tsk
    dtt = np.zeros(len(to), dtype = np.float64)
    for i in range(1,len(dtt)-1):
        dtt[i] = to[i+1] - to[i]
    dt = np.round(np.max(dtt),decimals=5)
    print('to:',to.shape,to[0],to[-1],' dt=',dt)
    fig, (ax1, ax2) = plt.subplots(2, 1)#, sharex=True)
    fig.set_size_inches(siz1*2, siz2*1.5)
    ax1.set_xlabel(r'$t$')
    ax1.set_title(r'x,y,z='+str(x[prb])+','+str(y[prb])+','+str(z[prb]))
    ax2.set_xlabel(r'$f$')
    tmin = t[0]
    tmax = t[-1]
    #ax1.set_xlim(tmin, tmax)
    #ax1.set_ylim(-0.001, 0.001)
    ax2.set_xlim(0., 1.0)
    #ax2.set_ylim(1.e-2, 1.e1)
    for f in range(0,fld):
        f += 1
        if f == 1:
            cor='k';nm='u';uo=u[itsk:-1]
            tmax = np.round(tmax)
            tn = np.linspace(tsk, tmax, int((tmax-tsk)/dt), endpoint=False)
            print('tn:',tn.shape,tn[0],tn[-1])
        if f == 2:
           cor='g';nm='v';uo=v[itsk:-1]
        if f == 3:
           cor='b';nm='w';uo=w[itsk:-1]
        #if f == 4:
        #   cor='magenta';nm='p';uo=p[itsk:-1]
        if f == 4:
           cor='r';nm='t';uo=tm[itsk:-1]           
        qur = np.mean(uo)
        uo -= qur
        un = interpolate.interp1d(to,uo,fill_value="extrapolate")(tn)
        ax1.plot(tn[::pskip],un[::pskip],c=cor,ls='-',lw=0.5)#,label=r'$x=2$')
        ufft = fft(un)
        ufreq = fftfreq(len(un), d=dt)
        ufmax = np.max(ufreq)
        prfreq = ufreq[np.argmax(np.abs(ufft))]
        print(nm,' mean=',round(qur,6),'peak freq,omega=',round(prfreq,6),round(2*pi*prfreq,4))
        ax2.semilogy(ufreq[ufreq > 0], np.abs(ufft)[ufreq > 0],c=cor,ls='-',lw=0.5,label='f_'+nm+'='+("{0:.4f}".format(round(prfreq,4))))
    plt.legend(loc='upper right',fontsize=8)
    fname = filename+'_hist_p'+'{:03d}'.format(prb+1)+'.'+formt
    print('Saving ',fname)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust)
    plt.close()
    plt.clf()

########################################
filename = "ray.his"
file = open(filename, 'r')
print('Opening file:',filename)
nps = int(file.readline().rstrip())
print('Number of probes found:',nps)
x=np.zeros(nps);y=np.zeros(nps);z=np.zeros(nps)
for n, line in zip(range(nps), file):
    x[n], y[n], z[n] = line.split(" ")
    #print(n,' at x,y,z=',x[n],y[n],z[n])
file.close()
data = np.loadtxt(filename, skiprows=nps+1).reshape(-1, nps, 6)
for probe in range(1,nps,1):
    fft_probe(prb=probe,fld=4,tsk=20.0)
