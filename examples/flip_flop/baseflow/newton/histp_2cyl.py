import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import style as style
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
#print(plt.style.available)
import numpy as np
import scipy as sp
import pandas as pd

from scipy.integrate import simps
from scipy import sin, cos, pi, exp, tanh, log
from scipy.special import erf, erfc
from scipy.optimize import curve_fit
from scipy import stats
from scipy.signal import butter, filtfilt, hilbert, chirp
from scipy import interpolate
from scipy.fftpack import fftfreq, fft
from scipy.signal import find_peaks

import os
import os.path
import csv
import time
import pickle

params = {'text.usetex': False,
          'font.size': 10,
          'legend.fontsize': 8,
          'legend.handlelength': 2.5}
plt.rcParams.update(params)
#':' dotted,'-.' dashdot, '--' dashed, '-' solid
#b: blue, g: green, r: red,c: cyan,m: magenta,y: yellow,k: black

transp = False
formt = 'png'
ajust = 'tight'
qual = 900
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
        self.a = data[::nprobes,1]
        self.b = data[::nprobes,2]
        self.c = data[::nprobes,3]
        #self.d = data[::nprobes,4]
        del data

class Residu(object):
    def __init__(self, filename):
        data = np.transpose(np.loadtxt(filename))
        self.t = data[0]
        self.p = data[1]
        self.r = data[2]
        self.rt = data[3]

class fft_signal(object):
    def __init__(self, filename):
        data = np.transpose(np.loadtxt(filename))
        self.i = data[0]
        self.t = data[1]
        self.u = data[2]

class strouhal(object):
    def __init__(self,flname):
        data = np.loadtxt(flname, unpack=True)
        self.t   = data[0]
        self.st  = data[1]
        self.dt  = data[2]
        del data

class poincare(object):
    def __init__(self,flname):
        #data = np.transpose(np.loadtxt(flname))
        data = np.loadtxt(flname, unpack=True)
        self.t   = data[0]
        self.u   = data[1]
        self.ud  = data[2]
        self.udd = data[3]
        del data

class period(object):
    def __init__(self,flname):
        data = np.loadtxt(flname, unpack=True)
        self.t   = data[0]
        self.period   = data[1]
        self.mean  = data[3]
        self.l2 = data[4]
        del data

s = strouhal('strouhal.dat')
fig=plt.figure();fig.set_size_inches(siz2, siz2)
plt.axhline(y=0.114, ls=':',lw=0.5,c='r',label=r'$St=0.114$')
plt.ylabel(r'$St$');plt.xlabel(r'$t$');plt.scatter(s.t[s.st != 0],s.st[s.st != 0],c='k',s=0.01)#;plt.legend(loc='upper left',fontsize=6)
fname='strouhal.'+formt;plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()

fig=plt.figure();fig.set_size_inches(siz2, siz2)
plt.ylabel(r'$dt$');plt.xlabel(r'$t$');plt.scatter(s.t,s.dt,c='k',s=0.01)#;plt.legend(loc='upper left',fontsize=6)
fname='strouhal_dt.'+formt;plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()

f = fft_signal('fft_signal.dat')
fig=plt.figure();fig.set_size_inches(siz2, siz2);plt.xlabel(r'$t$');
plt.scatter(f.t,f.u, c='k',s=0.01)
fname='fft_signal.'+formt;#plt.legend(loc='best',fontsize=6);
plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()

f = Residu('residu_tdf.dat')
fig=plt.figure();fig.set_size_inches(siz2, siz2);plt.xlabel(r'$t$');plt.yscale('log')#plt.xscale('log')
#plt.xlim(100,200);plt.ylim(1e-333)
plt.axhline(y=1.0e-2, xmin=0, xmax=1,ls=':',lw=0.5)
plt.axhline(y=1.0e-3, xmin=0, xmax=1,ls=':',lw=0.5)
plt.axhline(y=1.0e-4, xmin=0, xmax=1,ls=':',lw=0.5)
plt.axhline(y=1.0e-5, xmin=0, xmax=1,ls=':',lw=0.5)
plt.axhline(y=1.0e-6, xmin=0, xmax=1,ls=':',lw=0.5)
plt.axhline(y=1.0e-7, xmin=0, xmax=1,ls=':',lw=0.5)
plt.axhline(y=1.0e-8, xmin=0, xmax=1,ls=':',lw=0.5)
plt.axhline(y=1.0e-9, xmin=0, xmax=1,ls=':',lw=0.5)
plt.scatter(f.t,f.r,c='k',s=0.005)
fname='residu_tdf.'+formt;#plt.legend(loc='best',fontsize=6);
plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()


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

    poi=poincare('poincare.dat')
    poif=poincare('poincaref.dat')

    ign = int(0.20*len(poi.t))
    ig2 = int(0.30*len(poi.t))

    fig = plt.figure();fig.set_size_inches(siz1, siz2)
    plt.xlabel(r'$\dot{v}$');plt.ylabel(r'$v$')
    plt.plot(poi.ud[ign:ig2],poi.u[ign:ig2],c='k',ls='dotted',lw=0.1)
    plt.plot(poi.ud[ig2:-1],poi.u[ig2:-1],c='magenta',ls='dotted',lw=0.5)
    plt.axhline(y=0.,ls=':',lw=1.,c='k');plt.axvline(x=0.,ls=':',lw=1.,c='k')
    fname = 'histp'+str(prb)+'_projection_vd_v.'+formt; print('Saving fl ',fname)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);plt.close();plt.clf()

    fig = plt.figure();fig.set_size_inches(siz1, siz2)
    plt.xlabel(r'$\dot{v}$');plt.ylabel(r'$v$')
    plt.plot(poif.ud[ign:ig2],poif.u[ign:ig2],c='k',ls='dotted',lw=0.1)
    plt.plot(poif.ud[ig2:-1],poif.u[ig2:-1],c='magenta',ls='dotted',lw=0.5)
    plt.axhline(y=0.,ls=':',lw=1.,c='k');plt.axvline(x=0.,ls=':',lw=1.,c='k')
    fname = 'histp'+str(prb)+'_projection_vd_v_f.'+formt; print('Saving fl ',fname)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);plt.close();plt.clf()


    fig = plt.figure();fig.set_size_inches(siz1, siz2)
    plt.xlabel(r'$\ddot{v}$');plt.ylabel(r'$v$')
    plt.plot(poi.udd[ign:ig2],poi.u[ign:ig2],c='k',ls='dotted',lw=0.1)
    plt.plot(poi.udd[ig2:-1],poi.u[ig2:-1],c='b',ls='dotted',lw=0.5)
    plt.axhline(y=0.,ls=':',lw=1.,c='k');plt.axvline(x=0.,ls=':',lw=1.,c='k')
    fname = 'histp'+str(prb)+'_projection_vdd_v.'+formt;print('Saving fl ',fname)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);plt.close(); plt.clf()

    fig = plt.figure();fig.set_size_inches(siz1, siz2)
    plt.xlabel(r'$\ddot{v}$');plt.ylabel(r'$v$')
    plt.plot(poif.udd[ign:ig2],poif.u[ign:ig2],c='k',ls='dotted',lw=0.1)
    plt.plot(poif.udd[ig2:-1],poif.u[ig2:-1],c='b',ls='dotted',lw=0.5)
    plt.axhline(y=0.,ls=':',lw=1.,c='k');plt.axvline(x=0.,ls=':',lw=1.,c='k')
    fname = 'histp'+str(prb)+'_projection_vdd_v_f.'+formt;print('Saving fl ',fname)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);plt.close(); plt.clf()


    fig = plt.figure();fig.set_size_inches(siz1, siz2)
    plt.xlabel(r'$\ddot{v}$');plt.ylabel(r'$\dot{v}$')
    plt.plot(poi.udd[ign:ig2],poi.ud[ign:ig2],c='k',ls='dotted',lw=0.1)
    plt.plot(poi.udd[ig2:-1],poi.ud[ig2:-1],c='g',ls='dotted',lw=0.5)
    plt.axhline(y=0.,ls=':',lw=1.,c='k');plt.axvline(x=0.,ls=':',lw=1.,c='k')
    fname = 'histp'+str(prb)+'_projection_vdd_vd.'+formt;print('Saving fl ',fname)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);plt.close(); plt.clf()

    fig = plt.figure();fig.set_size_inches(siz1, siz2)
    plt.xlabel(r'$\ddot{v}$');plt.ylabel(r'$\dot{v}$')
    plt.plot(poif.udd[ign:ig2],poif.ud[ign:ig2],c='k',ls='dotted',lw=0.1)
    plt.plot(poif.udd[ig2:-1],poif.ud[ig2:-1],c='g',ls='dotted',lw=0.5)
    plt.axhline(y=0.,ls=':',lw=1.,c='k');plt.axvline(x=0.,ls=':',lw=1.,c='k')
    fname = 'histp'+str(prb)+'_projection_vdd_vd_f.'+formt;print('Saving fl ',fname)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);plt.close(); plt.clf()

    fig = plt.figure();fig.set_size_inches(siz1, siz2)
    plt.xlabel(r'$\ddot{v}$');plt.ylabel(r'$v$')
    plt.plot(poi.u[ign:ig2][poi.ud[ign:ig2] == 0.0], poi.udd[ign:ig2][poi.ud[ign:ig2] == 0.0],c='b',ls='dotted',lw=0.5)
    #plt.plot(poi.udd[ign:ig2],poi.u[ign:ig2],c='k',ls='dotted',lw=0.1)
    #plt.plot(poi.udd[ig2:-1],poi.u[ig2:-1],c='b',ls='dotted',lw=0.5)
    plt.axhline(y=0.,ls=':',lw=1.,c='k');plt.axvline(x=0.,ls=':',lw=1.,c='k')
    fname = 'histp'+str(prb)+'_section_vdd_v.'+formt;print('Saving fl ',fname)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);plt.close(); plt.clf()


    fig = plt.figure();fig.set_size_inches(siz1*2, siz2*2);ax = fig.gca(projection='3d')
    ax.set_xlabel(r'$v$');ax.set_ylabel(r'$\dot{v}$');ax.set_zlabel(r'$\ddot{v}$')
    #ax.view_init(60, 35)
    factor=1.5
    vmin = np.min(poi.u[ig2:-1])*factor;     vmax = np.max(poi.u[ig2:-1])*factor
    vdmin = np.min(poi.ud[ig2:-1])*factor;   vdmax = np.max(poi.ud[ig2:-1])*factor
    vddmin = np.min(poi.udd[ig2:-1])*factor; vddmax = np.max(poi.udd[ig2:-1])*factor
    ax.plot(poi.u[ign:ig2],poi.ud[ign:ig2],vddmin,                             c='k',ls='dotted',lw=0.1)
    ax.plot(poi.u[ig2:-1],poi.ud[ig2:-1],vddmin,                             c='magenta',ls='dotted',lw=0.5)
    ax.plot(vmin*np.ones(len(poi.t[ign:ig2])),poi.ud[ign:ig2],poi.udd[ign:ig2], c='k',ls='dotted',lw=0.1)
    ax.plot(vmin*np.ones(len(poi.t[ig2:-1])),poi.ud[ig2:-1],poi.udd[ig2:-1], c='g',ls='dotted',lw=0.5)
    ax.plot(poi.u[ign:ig2],vdmax*np.ones(len(poi.t[ign:ig2])), poi.udd[ign:ig2],c='k',ls='dotted',lw=0.1)
    ax.plot(poi.u[ig2:-1],vdmax*np.ones(len(poi.t[ig2:-1])), poi.udd[ig2:-1],c='b',ls='dotted',lw=0.5)
    ax.plot(poi.u[ign:ig2],poi.ud[ign:ig2],poi.udd[ign:ig2],c='k',ls='dotted',lw=0.2)
    ax.plot(poi.u[ig2:-1],poi.ud[ig2:-1],poi.udd[ig2:-1],c='r',ls='dotted',lw=0.4)
    ax.set_xlim(vmin,vmax);ax.set_ylim(vdmin,vdmax);ax.set_zlim(vddmin,vddmax)
    fname = 'histp'+str(prb)+'_attractor.'+formt;print('Saving fl ',fname)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);plt.close(); plt.clf()

    fig = plt.figure();fig.set_size_inches(siz1*2, siz2*2);ax = fig.gca(projection='3d')
    ax.set_xlabel(r'$v$');ax.set_ylabel(r'$\dot{v}$');ax.set_zlabel(r'$\ddot{v}$')
    #ax.view_init(60, 35)
    factor=1.5
    vmin = np.min(poif.u[ig2:-1])*factor;     vmax = np.max(poif.u[ig2:-1])*factor
    vdmin = np.min(poif.ud[ig2:-1])*factor;   vdmax = np.max(poif.ud[ig2:-1])*factor
    vddmin = np.min(poif.udd[ig2:-1])*factor; vddmax = np.max(poif.udd[ig2:-1])*factor
    ax.plot(poif.u[ign:ig2],poif.ud[ign:ig2],vddmin,                             c='k',ls='dotted',lw=0.1)
    ax.plot(poif.u[ig2:-1],poif.ud[ig2:-1],vddmin,                             c='magenta',ls='dotted',lw=0.5)
    ax.plot(vmin*np.ones(len(poif.t[ign:ig2])),poif.ud[ign:ig2],poif.udd[ign:ig2], c='k',ls='dotted',lw=0.1)
    ax.plot(vmin*np.ones(len(poif.t[ig2:-1])),poif.ud[ig2:-1],poif.udd[ig2:-1], c='g',ls='dotted',lw=0.5)
    ax.plot(poif.u[ign:ig2],vdmax*np.ones(len(poif.t[ign:ig2])), poif.udd[ign:ig2],c='k',ls='dotted',lw=0.1)
    ax.plot(poif.u[ig2:-1],vdmax*np.ones(len(poif.t[ig2:-1])), poif.udd[ig2:-1],c='b',ls='dotted',lw=0.5)
    ax.plot(poif.u[ign:ig2],poif.ud[ign:ig2],poif.udd[ign:ig2],c='k',ls='dotted',lw=0.2)
    ax.plot(poif.u[ig2:-1],poif.ud[ig2:-1],poif.udd[ig2:-1],c='r',ls='dotted',lw=0.4)
    ax.set_xlim(vmin,vmax);ax.set_ylim(vdmin,vdmax);ax.set_zlim(vddmin,vddmax)
    fname = 'histp'+str(prb)+'_attractor_f.'+formt;print('Saving fl ',fname)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);plt.close(); plt.clf()

    
    fig = plt.figure();fig.set_size_inches(siz1, siz2); plt.xlim(2000., 2050.)
    #plt.plot(d1.t[itsk:-1,d1.b[itsk:-1],c='k',ls='-',lw=1.5,label='his')
    plt.plot(poi.t,    poi.u,c='k',ls='-',lw=0.65,label='probe')

    deri1= np.gradient(poi.u, poi.t)
    deri2 = np.gradient(deri1, poi.t)

    a = 0.1; b=2*pi*0.11;c=2*pi*0.02
    #v = (a*(sin(b*poi.t) + sin(c*poi.t)))
    #vd1 = a*b*cos(b*poi.t) + a*c*cos(c*poi.t)
    vd2 = -a*(b**2)*sin(b*poi.t) - a*(c**2)*sin(c*poi.t)

    #plt.plot(poi.t,    v,c='k',ls='dotted',lw=0.65,label='probe')
    #plt.plot( poi.t, vd1, c='brown',ls='solid',lw=0.5,label='analitica')
    plt.plot( poi.t,deri1   ,c='r',ls='-',lw=1.,label='python')
    plt.plot( poi.t,poi.ud  ,c='g',ls='-',lw=1.,label='first')
    plt.plot( poi.t,poi.udd ,c='b',ls='-',lw=1.,label='second')
    #plt.plot(poif.t, poif.ud, c='b',ls='dotted',lw=0.8,label='first f')
    plt.legend(loc='upper right',fontsize=8)
    fname = 'histp'+str(prb)+'_check_deriv1.'+formt;print('Saving fl ',fname)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);plt.close(); plt.clf()


    fig = plt.figure();fig.set_size_inches(siz1, siz2); plt.xlim(2000., 2050.)
    plt.plot( poi.t, vd2, c='brown',ls='solid',lw=0.5,label='analitica')
    #plt.plot( poi.t, deri2 ,c='r',ls='solid',lw=0.25,label='python')
    plt.plot( poi.t,  poi.udd ,c='g',ls='dotted',lw=0.3,label='second')
    #plt.plot(poif.t, poif.udd,c='b',ls='dotted',lw=0.3,label='second f')
    plt.legend(loc='upper right',fontsize=8)
    fname = 'histp'+str(prb)+'_check_deriv2.'+formt;print('Saving fl ',fname)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);plt.close(); plt.clf()






    #https://matplotlib.org/gallery/lines_bars_and_markers/multicolored_line.html
   
    fig = plt.figure(); fig.set_size_inches(siz1, siz2)
    plt.xlabel(r'$v$');plt.ylabel(r'$u$')
    plt.plot(d1.b[ign:ig2],d1.a[ign:ig2],c='k',ls='dotted',lw=0.05)
    plt.plot(d1.b[ig2:-1],d1.a[ig2:-1],c='r',ls='dotted',lw=0.2); plt.ylim(1.08,1.15)
    #plt.axhline(y=1.,ls=':',lw=1.);plt.axvline(x=0.,ls=':',lw=1.)
    fname = 'histp'+str(prb)+'_map.'+formt;print('Saving fl ',fname)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);plt.close();plt.clf()
 
    fig, (ax1, ax2) = plt.subplots(2, 1)#, sharex=True)
    fig.set_size_inches(siz1*2, siz2*1.5)
    ax1.set_xlabel(r'$t$')
    #ax1.set_yscale('log')
    #ax1.set_xscale('log')
    ax1.set_title(ttl)
    ax2.set_xlabel(r'$f$')
    tmin = d1.t[0]
    tmax = d1.t[-1]
    #ax1.set_xlim(2000, 2100)
    #ax1.set_ylim(0., 0.0015)
    ax2.set_xlim(0., 0.3)
    ax2.set_ylim(1., 1.0e4)
    #ax1.set_xlim(0., 40.)
    for f in range(0,fld):
        f += 1
        if f == 1:
            cor='g';nm='v';
            uo=d1.b[itsk:-1]
            tmax = np.round(tmax)
            tn = np.linspace(tsk, tmax, int((tmax-tsk)/dt), endpoint=False)
            print('tn:',tn.shape,tn[0],tn[-1])
        if f == 2:
            cor='r';nm='v filt';
            uo=poif.u[itsk:-1]
            to=poif.t[itsk:-1]
            dtt = np.zeros(len(to), dtype = np.float64)
            for i in range(1,len(dtt)-1):
                dtt[i] = to[i+1] - to[i]
            dt = np.round(np.max(dtt),decimals=5)
            tmax = np.round(tmax)
            tn = np.linspace(tsk, tmax, int((tmax-tsk)/dt), endpoint=False)
        qur = np.mean(uo)
        print(nm,' mean=',qur)
        uo -= qur
        un = interpolate.interp1d(to,uo,fill_value="extrapolate")(tn)
        ax1.plot(tn[::pskip],un[::pskip],c=cor,ls='-',lw=0.5,label=nm)
        ax1.axhline(y=np.max(un), xmin=0, xmax=1,ls=':',lw=1.)
        ax1.axhline(y=np.min(un), xmin=0, xmax=1,ls=':',lw=1.)
        ufft = fft(un)
        ufreq = fftfreq(len(un), d=dt)
        ufmax = np.max(ufreq)
        prfreq = ufreq[np.argmax(np.abs(ufft))]
        print('pfreq=',prfreq)
        peaks = find_peaks(np.abs(ufft)[ufreq > 0], height=900)[0]
        print(peaks,ufreq[peaks[:]])
        ax2.plot(np.abs(ufft)[peaks],peaks,'x')
        ax2.semilogy(ufreq[ufreq > 0], np.abs(ufft)[ufreq > 0],c=cor,ls='-',lw=0.5)#,label='f_'+nm+'='+("{0:.4f}".format(round(prfreq,4))))
    ax2.axvline(x=0.11250011,ls=':',lw=1.,c='red',label='f=0.1125')
    ax2.axvline(x=0.01900002,ls=':',lw=1.,c='magenta',label='f=0.0190')
    ax1.legend(loc='upper right',fontsize=8)
    ax2.legend(loc='upper right',fontsize=8)
    fname = 'histp'+str(prb)+'.'+formt;print('Saving fl ',fname)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);plt.close();plt.clf()

fft_probe(fl='2cyl.his',tsk=2000.,nps=1,prb=1,fld=2,ttl=r'$Re=67, g=0.7 \quad x,y,z=1,0,0$')
#fft_probe(fl='tpjet.his',tsk=0.,nps=1,prb=1,fld=2,ttl=r'$Re=2000 \quad x,y,z=5,1,0$')
