import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
#from matplotlib import style as style
#from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
#                               AutoMinorLocator)
#print(plt.style.available)
import matplotlib as mpl
import numpy as np
import scipy as sp
#import pandas as pd

#from scipy.integrate import simps
from scipy import sin, cos, pi, exp, tanh, log
#from scipy.special import erf, erfc
#from scipy.optimize import curve_fit
#from scipy import stats
#from scipy.signal import butter, filtfilt, hilbert, chirp
from scipy import interpolate
from scipy.fftpack import fftfreq, fft
#from scipy.signal import find_peaks

import os
import os.path
import csv
import time
import pickle

params = {'text.usetex': False,
          'font.size': 12,
          #'font.family': 'serif',
          #'text.latex.unicode': True,
          'legend.fontsize': 10,
          'legend.handlelength': 2.5}
plt.rcParams.update(params)
#':' dotted,'-.' dashdot, '--' dashed, '-' solid
#b: blue, g: green, r: red,c: cyan,m: magenta,y: yellow,k: black

transp = False
formt = 'png'
ajust = 'tight'
qual = 300
siz1 = 3.5
siz2 = 2.45

class His3D(object):
    def __init__(self, probe, nprobes, flname):
        self.probe = probe
        self.nprobes = nprobes
        data = np.genfromtxt(flname,skip_header=1+nprobes)
        data = np.asarray(data)
        print(data.shape)
        self.t = data[::nprobes,0]
        self.u = data[::nprobes,1]
        self.v = data[::nprobes,2]
        self.w = data[::nprobes,3]
        try:
          self.p = data[::nprobes,4]
        except:
          self.p = 0.
        del data

pskip = 1

def fft_probe(fl,tsk,nps,prb,fld,ttl):
    print('FFT Probe ',prb)
    d1=His3D(prb,nps,fl)
    itsk = np.where(d1.t > tsk)[0][0]
    to = d1.t[itsk:-1]# - tsk
    print('to:',to.shape,to[0],to[-1])
    dtt = np.zeros(len(to), dtype = np.float64)
    for i in range(1,len(dtt)-1):
        dtt[i] = to[i+1] - to[i]
    dt = np.round(np.max(dtt),decimals=5)
    print('dt=',dt)
    fig, (ax1, ax2) = plt.subplots(2, 1)#, sharex=True)
    fig.set_size_inches(siz1*2, siz2*1.5)
    ax1.set_xlabel(r'$t$')
    ax1.set_title(ttl)
    ax2.set_xlabel(r'$f$')
    tmin = d1.t[0]
    tmax = d1.t[-1]
    #ax1.set_xlim(tmin, tmax)
    #ax1.set_ylim(-0.001, 0.001)
    ax2.set_xlim(0., 0.4)
    #ax2.set_ylim(1.e-2, 1.e1)
    for f in range(0,fld):
        f += 1
        if f == 1:
            cor='k';nm='u';uo=d1.v[itsk:-1]
            tmax = np.round(tmax)
            tn = np.linspace(tsk, tmax, int((tmax-tsk)/dt), endpoint=False)
            print('tn:',tn.shape,tn[0],tn[-1])
        if f == 2:
           cor='g';nm='v';uo=d1.v[itsk:-1]
        if f == 3:
           cor='b';nm='w';uo=d1.w[itsk:-1]
        if f == 4:
           cor='r';nm='p';uo=d1.p[itsk:-1]
        qur = np.mean(uo)
        print(nm,' mean=',qur)
        uo -= qur
        un = interpolate.interp1d(to,uo,fill_value="extrapolate")(tn)
        ax1.plot(tn[::pskip],un[::pskip],c=cor,ls='-',lw=0.5)#,label=r'$x=2$')
        ufft = fft(un)
        ufreq = fftfreq(len(un), d=dt)
        ufmax = np.max(ufreq)
        prfreq = ufreq[np.argmax(np.abs(ufft))]
        print('pfreq=',prfreq)
        ax2.semilogy(ufreq[ufreq > 0], np.abs(ufft)[ufreq > 0],c=cor,ls='-',lw=0.5,label='f_'+nm+'='+("{0:.4f}".format(round(prfreq,4))))
    plt.legend(loc='upper right',fontsize=8)
    fname = fl+'_hist_p'+str(prb)+'.'+formt
    print('Saving fl ',fname)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust)
    plt.close()
    plt.clf()

fft_probe(fl='1cyl.his',tsk=850.,nps=1,prb=1,fld=1,ttl=r'$x,y,z=0.25,0.25,0.$')

