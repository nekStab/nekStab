#!/usr/bin/env python
#from pymech import neksuite as ns
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
fig_width=3.5
fig_height=2.4
dtmult = 1
def fft_probe(prb,fld,tsk):
    pskip = 1
    print()
    print('FFT Probe',prb,'at x,y,z='+str(x[prb-1])+','+str(y[prb-1])+','+str(z[prb-1]))
    prb-=1
    t=data[:,prb,0]
    u=data[:,prb,1]
    v=data[:,prb,2]
    w=data[:,prb,3]
    #p=data[:,prb,4]
    if tsk == 0:
        tsk =  t[0]
    itsk = np.where(t > tsk)[0][0]
    to = t[itsk:-1]# - tsk
    dtt = np.zeros(len(to), dtype = np.float64)
    for i in range(1,len(dtt)-1):
        dtt[i] = to[i+1] - to[i]
    dt = np.round(np.max(dtt),decimals=5)*dtmult
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
    ax2.set_xlim(0., 2.0)
    #ax2.set_ylim(1.e-2, 1.e1)
    for f in range(0,fld):
        f += 1
        if f == 1:
            cor='k';nm='v';uo=v[itsk:-1]
            tmax = np.round(tmax)
            tn = np.linspace(tsk, tmax, int((tmax-tsk)/dt), endpoint=False)
            print('tn:',tn.shape,tn[0],tn[-1])
        if f == 2:
           cor='g';nm='v';uo=v[itsk:-1]
        if f == 3:
           cor='b';nm='w';uo=w[itsk:-1]
        qur = np.mean(uo)
        uo -= qur
        un = interpolate.interp1d(to,uo,fill_value="extrapolate")(tn)
        ax1.plot(tn[::pskip],un[::pskip],c=cor,ls='-',lw=0.2)#,label=r'$x=2$')
        ufft = fft(un)
        ufreq = fftfreq(len(un), d=dt)
        ufmax = np.max(ufreq)
        prfreq = ufreq[np.argmax(np.abs(ufft))]
        print(nm,' mean=',round(qur,6),'peak freq,omega=',round(prfreq,6),round(2*pi*prfreq,4))
        ax2.semilogy(ufreq[ufreq > 0], np.abs(ufft)[ufreq > 0],c=cor,ls='-',lw=0.2,label='f_'+nm+'='+("{0:.4f}".format(round(prfreq,4))))
    plt.legend(loc='upper right',fontsize=8)
    fname = '1cav_hist_p'+'{:03d}'.format(prb+1)+'.'+formt
    print('Saving ',fname)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust)
    plt.close()
    plt.clf()

########################################

ratios = [
'',
]
tskps = [
20,
]

fld = ''
filename = "naca4412.his"

for ratio in range(0,len(ratios)):
    filen = fld + ratios[ratio] + '' + filename
    print('Opening file',filen)
    file = open(filen, 'r')
    nps = int(file.readline().rstrip())
    print('Number of probes found:',nps)
    x=np.zeros(nps);y=np.zeros(nps);z=np.zeros(nps)
    for n, line in zip(range(nps), file):
        x[n], y[n], z[n] = line.split(" ")
    file.close()
    data = np.loadtxt(filen, skiprows=nps+1).reshape(-1, nps, 5)

    for prb in range(0,nps,1):
        print('Probe number ',prb)
        fig, axs = plt.subplots(4, sharex=True)
        fig.set_size_inches(siz1*3, siz2*2)
        axs[0].set_xlabel(r'$t$')
        axs[0].set_ylabel(r'$u$')
        axs[1].set_ylabel(r'$v$')
        axs[2].set_ylabel(r'$w$')
        axs[3].set_ylabel(r'$p$')
        axs[0].set_title(r'R=0.'+ratios[ratio]+' x,y,z='+str(x[prb])+','+str(y[prb])+','+str(z[prb]))
        t=data[:,prb,0]
        u=data[:,prb,1]
        v=data[:,prb,2]
        w=data[:,prb,3]
        p=data[:,prb,4]
    
        #if tskps[ratio] == 0:
        #   tsk =  t[0]
        #else:
        #   tsk = tskps[ratio] 
        tsk = tskps[ratio] 
        tmax = t[-1]
        
        itsk = np.where(t > tsk)[0][0]
        to = t[itsk:-1]# - tsk
        dtt = np.zeros(len(to), dtype = np.float64)
        for i in range(1,len(dtt)-1):
           dtt[i] = to[i+1] - to[i]
        dt = np.round(np.max(dtt),decimals=6)*dtmult
        print('to:',to.shape,to[0],to[-1],' dt=',dt)
        tmax = np.round(tmax)
        tn = np.linspace(tsk, tmax, int((tmax-tsk)/dt), endpoint=False)
        print('tn:',tn.shape,tn[0],tn[-1])

        axs[0].plot(t,u,c='k',ls='-',lw=1)
        axs[1].plot(t,v,c='g',ls='-',lw=1)
        axs[2].plot(t,w,c='b',ls='-',lw=1)
        axs[3].plot(t,p,c='r',ls='-',lw=1)

        #plt.legend(loc='upper right',fontsize=6)
        fname = ratios[ratio]+'probe'+'{:03d}'.format(prb+1)+'.'+formt
        print('Saving ',fname)
        plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust)
        plt.close()
        plt.clf()

        fig=plt.figure();fig.set_size_inches(fig_width*1.5, fig_height)
        plt.xlabel(r'$St$')

        uo = u[itsk:-1]; nm='u'; cor='k'
        dt = t[1] - t[0]
        qur = np.mean(uo)
        uo -= qur
        un = interpolate.interp1d(to,uo,fill_value="extrapolate")(tn)
        ufft = fft(un)
        ufreq = fftfreq(len(un), d=dt)
        ufmax = np.max(ufreq)
        prfreq = ufreq[np.argmax(np.abs(ufft))]
        print(nm,' mean=',round(qur,6),uo.min(),uo.max())
        print('peak freq,omega=',round(prfreq,6),round(2*pi*prfreq,4))
        plt.semilogy(ufreq[ufreq>0],np.abs(ufft)[ufreq>0],c=cor,lw=0.7,label=nm+' St='+("{0:.4f}".format(round(prfreq,4))))

        plt.axvline(x=prfreq    ,lw=0.4, color='k', ls=':')
        plt.axvline(x=prfreq*2  ,lw=0.2, color='k', ls=':')
        plt.axvline(x=prfreq*3  ,lw=0.2, color='k', ls=':')
        plt.axvline(x=prfreq*4  ,lw=0.2, color='k', ls=':')
        plt.axvline(x=prfreq*5  ,lw=0.2, color='k', ls=':')
        plt.axvline(x=prfreq*6  ,lw=0.2, color='k', ls=':')
        plt.axvline(x=prfreq*7  ,lw=0.2, color='k', ls=':')
        plt.axvline(x=prfreq*8  ,lw=0.2, color='k', ls=':')
        plt.axvline(x=prfreq*0.5,lw=0.2, color='k', ls=':')

        uo = v[itsk:-1]; nm='v'; cor='g'
        qur = np.mean(uo)
        uo -= qur
        un = interpolate.interp1d(to,uo,fill_value="extrapolate")(tn)
        ufft = fft(un)
        ufreq = fftfreq(len(un), d=dt)
        ufmax = np.max(ufreq)
        prfreq = ufreq[np.argmax(np.abs(ufft))]
        print(nm,' mean=',round(qur,6),uo.min(),uo.max())
        print('peak freq,omega=',round(prfreq,7),round(2*pi*prfreq,4))
        plt.semilogy(ufreq[ufreq>0],np.abs(ufft)[ufreq>0],c=cor,lw=0.6,label=nm+' St='+("{0:.4f}".format(round(prfreq,4))))

        uo = w[itsk:-1]; nm='w'; cor='b'
        qur = np.mean(uo)
        uo -= qur
        un = interpolate.interp1d(to,uo,fill_value="extrapolate")(tn)
        ufft = fft(un)
        ufreq = fftfreq(len(un), d=dt)
        ufmax = np.max(ufreq)
        prfreq = ufreq[np.argmax(np.abs(ufft))]
        print(nm,' mean=',round(qur,6),uo.min(),uo.max())
        print('peak freq,omega=',round(prfreq,6),round(2*pi*prfreq,4))
        plt.semilogy(ufreq[ufreq>0],np.abs(ufft)[ufreq>0],c=cor,lw=0.4,label=nm+' St='+("{0:.4f}".format(round(prfreq,4))))

        uo = p[itsk:-1]; nm='p'; cor='r'
        qur = np.mean(uo)
        uo -= qur
        un = interpolate.interp1d(to,uo,fill_value="extrapolate")(tn)
        ufft = fft(un)
        ufreq = fftfreq(len(un), d=dt)
        ufmax = np.max(ufreq)
        prfreq = ufreq[np.argmax(np.abs(ufft))]
        print(nm,' mean=',round(qur,6),uo.min(),uo.max())
        print('peak freq,omega=',round(prfreq,6),round(2*pi*prfreq,4))
        plt.semilogy(ufreq[ufreq>0],np.abs(ufft)[ufreq>0],c=cor,lw=0.2,label=nm+' St='+("{0:.4f}".format(round(prfreq,4))))

        plt.xlim(0., 3.0)
        plt.ylim(1e-5,1e5)

        plt.legend(loc='upper right',fontsize=6)
        fname = ratios[ratio]+'probe'+'{:03d}'.format(prb+1)+'fft.'+formt
        print('Saving ',fname)
        plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust)
        plt.close()
        plt.clf()
