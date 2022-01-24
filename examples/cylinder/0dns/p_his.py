#!/usr/bin/env python
#from pymech import neksuite as ns
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

params = {'text.usetex': False,'font.size': 8,'legend.fontsize': 8,'legend.handlelength': 2.5,
'agg.path.chunksize':100000
};plt.rcParams.update(params)
transp = False
formt = 'png'
ajust = 'tight'
qual = 900
siz2 = 5.33
siz1 = 16.0*siz2/9.0

#large-sized journals 3.3 double-column, 6.85 (single-column) wide and not higher 9.2126.
#small-sized journals, 4.685 wide and not higher than 7.677 mm.


lspec = 0.5 # width spctra lines
smarker = 0.001 # width signal markers
        
def fft(u,t,dt):
    ufft = rfft(u)
    ufreq = rfftfreq(len(u), d=dt)
    prfreq = ufreq[np.argmax(np.abs(ufft))]
    print('Peak at f=',round(prfreq,6))
    #print('        w=',round(2*pi*prfreq,4))
    #print('        T=',round(1/prfreq,6))
    return ufft,ufreq,prfreq

#files = sorted(glob.glob('*/*.his', recursive=True))
files = [
'1cyl.his',
]
tskps = [
0,
]
tmaxs = [
0,
]

i = 0
for filen in files:
    #filen += '/jcf.his'
    print('Opening file',filen)
    file = open(filen, 'r')
    fld = filen.split('/')[0]
    nps = int(file.readline().rstrip())
    print('Number of probes found:',nps)
    x=np.zeros(nps);y=np.zeros(nps);z=np.zeros(nps)
    for n, line in zip(range(nps), file):
        x[n], y[n], z[n] = line.split()
    file.close()
    try:
        data = np.loadtxt(filen, skiprows=nps+1).reshape(-1, nps, 5)
        if3d = True
    except:
        data = np.loadtxt(filen, skiprows=nps+1).reshape(-1, nps, 4)
        if3d = False

    #nps = 1
    for prb in range(0,nps,1):
        print('Probe number ',prb+1)
        #fig, axs = plt.subplots(4, sharex=True)
        fig, axs = plt.subplots(2, sharex=False)
        fig.set_size_inches(siz1, siz2)
        axs[0].set_xlabel(r'$t$')
        #axs[0].set_ylabel(r'$q$')

        axs[0].set_title(r''+' x,y,z='+str(x[prb])+','+str(y[prb])+','+str(z[prb]))
        t=data[:,prb,0][3:]
        u=data[:,prb,1][3:]
        v=data[:,prb,2][3:]
        if if3d:
            w=data[:,prb,3][3:]
            p=data[:,prb,4][3:]
        else:
            p=data[:,prb,3][3:]


        tmin=t.min();tmax=t.max()
        print('Original time series from:',tmin,tmax)


        if float(tskps[i]) > 0 and float(tskps[i]) > tmin:
            tmin = float(tskps[i])

        if float(tmaxs[i]) > 0 and float(tmaxs[i]) < tmax:
            tmax = float(tmaxs[i])
 
        print('Cropped time series from:',tmin,tmax)

        itmin = np.where(t >= tmin)[0][0]
        #print(itmin,t[itmin],t[0])

        itmax = np.where(t >= tmax)[0][0]
        #print(itmax,t[itmax],t[-1])

        to = t[itmin:itmax]

        dtt = np.zeros(len(to), dtype = np.float64)
        for xx in range(1,len(dtt)-1): #compute dt array
            dtt[xx] = to[xx+1] - to[xx]
    
        print(to.shape,dtt.shape)
        dt = dtt.max() #np.round(dtt.max(),decimals=7) # get maximum dt of array to be safe
        #print('t old:',to.shape,to.min(),to.max(),' dt=',dt)

        tn = np.linspace(t[itmin], t[itmax], int((t[itmax]-t[itmin])/dt), endpoint=False)
        # endpoint = False so new interval is inside old -- we do not want to extrapolate!
        #print('t new:',tn.shape,tn.min(),tn.max())
        #print('err (old-new) min max=',to.min()-tn.min(),to.max()-tn.max())

        axs[1].set_xlabel(r'$St$')

        uo = u[itmin:itmax]; nm='u'; cor='r'
        qur = np.mean(uo); uo -= qur; print('(',nm,') mean, min, max=',round(qur,6),uo.min(),uo.max())
        un = interpolate.interp1d(to,uo,kind='slinear',fill_value="extrapolate")(tn) #'zero','slinear','quadratic','cubic'=spline interpolation of 0th,1st,2nd,3rd;
        u = un
        #axs[0].scatter(tn,un,c=cor,ls='-',s=smarker)
        ufft,ufreq,prfreq = fft(un,tn,dt)

        #axs[1].semilogy(ufreq,np.abs(ufft),c=cor,lw=lspec,label=nm+' St='+("{0:.4f}".format(round(prfreq,4))))
        #axs[1].axvline(x=prfreq           ,c=cor,lw=lspec,ls=':')

        uo = v[itmin:itmax]; nm='v'; cor='g'
        qur = np.mean(uo); uo -= qur; print('(',nm,') mean, min, max=',round(qur,6),uo.min(),uo.max())
        un = interpolate.interp1d(to,uo,kind='slinear',fill_value="extrapolate")(tn) #'zero','slinear','quadratic','cubic'=spline interpolation of 0th,1st,2nd,3rd;
        ud = un

      
        axs[0].scatter(tn,un,c=cor,ls='-',s=smarker)
        ufft,ufreq,prfreq = fft(un,tn,dt)
        axs[1].semilogy(ufreq,np.abs(ufft),c=cor,lw=lspec)#,label=nm+' St='+("{0:.4f}".format(round(prfreq,4))))

        axs[1].axvline(x=prfreq         ,c=cor,lw=0.5,ls='--',label='$St=$'+("{0:.4f}".format(round(prfreq,4))))
   
        if if3d:
            uo = w[itmin:itmax]; nm='w'; cor='b'
        else:
            uo = p[itmin:itmax]; nm='p'; cor='b'
        qur = np.mean(uo); uo -= qur; print('(',nm,') mean, min, max=',round(qur,6),uo.min(),uo.max())
        un = interpolate.interp1d(to,uo,kind='slinear',fill_value="extrapolate")(tn) #'zero','slinear','quadratic','cubic'=spline interpolation of 0th,1st,2nd,3rd;
        udd = un
        
        #axs[0].scatter(tn,un,c=cor,ls='-',s=smarker)
        #ufft,ufreq,prfreq = fft(un,tn,dt)
        #axs[1].semilogy(ufreq,np.abs(ufft),c=cor,lw=lspec,label=nm+' St='+("{0:.4f}".format(round(prfreq,4))))
        #axs[1].axvline(x=prfreq           ,c=cor,lw=lspec,ls=':')
        
        axs[1].set_xlim(0., 0.4)
        axs[1].set_ylim(1e1,1e5) 

        #plt.ylim(1e-5,1e5) 

        plt.legend(loc='upper right')
        fname = fld+'_his'+str(prb+1)+'.'+formt
        print('Saving ',fname); print()
        plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust)
        plt.close(); plt.clf()
        '''
        t = tn
        fig = plt.figure();fig.set_size_inches(siz1*2, siz2*2);ax = fig.gca(projection='3d')
        ax.set_title(r''+ratios[i]+' x,y,z='+str(x[prb])+','+str(y[prb])+','+str(z[prb]))
        ax.set_xlabel(r'$u$');ax.set_ylabel(r'$v$');ax.set_zlabel(r'$w$')
        #ax.view_init(60, 35);
        factor=1.1
        vmin = np.min(u)*factor;     vmax = np.max(u)*factor
        vdmin = np.min(ud)*factor;   vdmax = np.max(ud)*factor
        vddmin = np.min(udd)*factor; vddmax = np.max(udd)*factor

        ax.plot(u,ud,vddmin, c='r',ls='dotted',lw=0.4)
        ax.plot(vmin*np.ones(len(t)),ud,udd, c='g',ls='dotted',lw=0.4)
        ax.plot(u,vdmax*np.ones(len(t)),udd,c='b',ls='dotted',lw=0.4)
        ax.plot(u,ud,udd,c='k',ls='solid',lw=0.8)

        ax.set_xlim(vmin,vmax);ax.set_ylim(vdmin,vdmax);ax.set_zlim(vddmin,vddmax)
        fname = fld+'_'+ratios[i]+'his'+str(prb+1)+'_attractor.'+formt
        print('Saving ',fname)
        plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);plt.close(); plt.clf()

        u = savgol_filter(ud, 3, 1, deriv=0, mode='nearest')
        ud = savgol_filter(ud, 15, 6, deriv=1, mode='nearest')
        udd = savgol_filter(ud, 23, 9, deriv=2, mode='nearest')
        dd = 6
        t = tn[dd:-dd];u = u[dd:-dd];ud = ud[dd:-dd]; udd = udd[dd:-dd]

        fig = plt.figure();fig.set_size_inches(siz1*2, siz2*2);ax = fig.gca(projection='3d')
        #ax.set_title(r'R=0.'+ratios[i]+' x,y,z='+str(x[prb])+','+str(y[prb])+','+str(z[prb]))
        ax.set_title(r'x,y,z='+str(x[prb])+','+str(y[prb])+','+str(z[prb]))
        ax.set_xlabel(r'$v$');ax.set_ylabel(r'$\dot{v}$');ax.set_zlabel(r'$\ddot{v}$')
        #ax.view_init(60, 35);
        factor=1.1
        vmin = np.min(u)*factor;     vmax = np.max(u)*factor
        vdmin = np.min(ud)*factor;   vdmax = np.max(ud)*factor
        vddmin = np.min(udd)*factor; vddmax = np.max(udd)*factor

        #ax.plot(u,ud,vddmin, c='r',ls='dotted',lw=0.4)
        #ax.plot(vmin*np.ones(len(t)),ud,udd, c='g',ls='dotted',lw=0.4)
        #ax.plot(u,vdmax*np.ones(len(t)),udd,c='b',ls='dotted',lw=0.4)
        ax.plot(u,ud,udd,c='k',ls='solid',lw=0.8)

        ax.set_xlim(vmin,vmax);ax.set_ylim(vdmin,vdmax);ax.set_zlim(vddmin,vddmax)
        fname = fld+'_'+ratios[i]+'his'+str(prb+1)+'_attractor_rec.'+formt
        print('Saving ',fname)
        plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);plt.close(); plt.clf()
        '''
    i += 1


