#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.fft import rfftfreq, rfft
params = {'text.usetex': False,
          'font.size': 10,
          'legend.fontsize': 6,
          'legend.handlelength': 1.3,
          'agg.path.chunksize':100000}
plt.rcParams.update(params)

formt = 'png'
ajust = 'tight'
qual = 400
fig_width = 3.4
fig_height = 4.3
        
def fft(u,dt):
    ufft = rfft(u)
    ufreq = rfftfreq(len(u), d=dt)
    prfreq = ufreq[np.argmax(np.abs(ufft))]
    print('Peak at f=',round(prfreq,6))
    return ufft,ufreq,prfreq

files = ['1cyl.his']

# set the skipping time
tskps = [20] # set to 0 to plot all

# set the maximum time 
tmaxs = [0] # set to 0 to plot all

for i, filen in enumerate(files):
    print(f'Opening file {filen}')
    
    with open(filen, 'r') as file:
        nps = int(file.readline().rstrip())
        print(f'Number of probes found: {nps}')

        coords = np.zeros((nps, 3))
        for n, line in zip(range(nps), file):
            coords[n] = line.split()

    #try:
    data = np.loadtxt(filen, skiprows=nps + 1).reshape(-1, nps, 8)
    #     is_3d = True
    # except:
    #     data = np.loadtxt(filen, skiprows=nps + 1).reshape(-1, nps, 7)
    #     is_3d = False
            
    for prb in range(0,nps,1):
        print(f'Probe number {prb + 1}')

        t=data[:,prb,0]
        u=data[:,prb,1]
        v=data[:,prb,2]
        # if is_3d:
        w=data[:,prb,3]
        p=data[:,prb,4]
        te=data[:,prb,5]
        s1=data[:,prb,6]
        s2=data[:,prb,7]
        # else:
        #     p=data[:,prb,3][1:]
        #     te=data[:,prb,4][1:]
        #     s1=data[:,prb,5][1:]
        #     s2=data[:,prb,6][1:]


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
        

        nm='u' # name of the signal
        cor='k' # color of the signal
        vo = u[itmin:itmax]
        qur = np.mean(vo) # compute mean
        vo -= qur # remove mean
        print('(',nm,') mean, min, max=',round(qur,6),vo.min(),vo.max())
        vn = interpolate.interp1d(to,vo,kind='slinear',fill_value="extrapolate")(tn) #'zero','slinear','quadratic','cubic'=spline interpolation of 0th,1st,2nd,3rd;
        axs[0].scatter(tn,vn,c=cor,ls='-',s=0.02) # plot signal
        ufft,ufreq,prfreq = fft(vn,dt) # compute fft
        # axs[1].plot(ufreq,np.abs(ufft),c=cor,lw=0.6)#,label=nm+' St='+("{0:.4f}".format(round(prfreq,4))))
        # axs[1].axvline(x=1*prfreq,c=cor,lw=0.5,ls='--',label=' $St=$'+("{0:.4f}".format(round(1*prfreq,4))))

        # We focus on the vertical velocity signal
        nm='v' # name of the signal
        cor='r' # color of the signal
        vo = v[itmin:itmax]
        qur = np.mean(vo) # compute mean
        vo -= qur # remove mean
        print('(',nm,') mean, min, max=',round(qur,6),vo.min(),vo.max())
        vn = interpolate.interp1d(to,vo,kind='slinear',fill_value="extrapolate")(tn) #'zero','slinear','quadratic','cubic'=spline interpolation of 0th,1st,2nd,3rd;
        axs[0].scatter(tn,vn,c=cor,ls='-',s=0.02) # plot signal
        ufft,ufreq,prfreq = fft(vn,dt) # compute fft
        axs[1].plot(ufreq,np.abs(ufft),c=cor,lw=0.6)#,label=nm+' St='+("{0:.4f}".format(round(prfreq,4))))
        axs[1].axvline(x=1*prfreq,c=cor,lw=0.5,ls='--',label=' $St=$'+("{0:.4f}".format(round(1*prfreq,4))))
        axs[1].axvline(x=3*prfreq,c=cor,lw=0.5,ls='--',label='$3St=$'+("{0:.4f}".format(round(3*prfreq,4))))
        axs[1].axvline(x=5*prfreq,c=cor,lw=0.5,ls='--',label='$5St=$'+("{0:.4f}".format(round(5*prfreq,4))))
        axs[1].axvline(x=7*prfreq,c=cor,lw=0.5,ls='--',label='$7St=$'+("{0:.4f}".format(round(7*prfreq,4))))


        nm='w' # name of the signal
        cor='g' # color of the signal
        vo = w[itmin:itmax]
        qur = np.mean(vo) # compute mean
        vo -= qur # remove mean
        print('(',nm,') mean, min, max=',round(qur,6),vo.min(),vo.max())
        vn = interpolate.interp1d(to,vo,kind='slinear',fill_value="extrapolate")(tn) #'zero','slinear','quadratic','cubic'=spline interpolation of 0th,1st,2nd,3rd;
        axs[0].scatter(tn,vn,c=cor,ls='-',s=0.02) # plot signal
        ufft,ufreq,prfreq = fft(vn,dt) # compute fft
        axs[1].plot(ufreq,np.abs(ufft),c=cor,lw=0.6)#,label=nm+' St='+("{0:.4f}".format(round(prfreq,4))))
        axs[1].axvline(x=1*prfreq,c=cor,lw=0.5,ls='--',label=' $St=$'+("{0:.4f}".format(round(1*prfreq,4))))
        axs[1].axvline(x=2*prfreq,c=cor,lw=0.5,ls='--',label=' $2St=$'+("{0:.4f}".format(round(2*prfreq,4))))
        axs[1].axvline(x=4*prfreq,c=cor,lw=0.5,ls='--',label=' $4St=$'+("{0:.4f}".format(round(4*prfreq,4))))

        axs[1].set_xscale('log')
        axs[1].set_yscale('log')
        axs[1].set_xlabel(r'$St$')
        axs[1].set_xlim(0.04, 4)
        axs[1].set_ylim(1e-1,1e3) 

        plt.legend(loc='upper left')
        fname = 'his'+str(prb+1)+'_fft.'+formt
        print('Saving ',fname); print()
        plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust)
        plt.close(); plt.clf()
        
        
        
        # log 
        fig = plt.figure(figsize=(fig_height, fig_width))
        plt.xlabel(r'$t$')
        plt.ylabel(r'$v$')
        plt.yscale('log')
        #plt.xscale('log')
        plt.plot(tn,vn,c='r',ls='-',lw=0.8)
        
        from scipy.signal import hilbert
        vn_hilbert = np.abs(hilbert(vn))
        plt.plot(tn[1:-1], vn_hilbert[1:-1], c='b', ls='--', lw=0.8)
        
        plt.xlim(0,200)
        plt.ylim(1e-4,1)

        fname = 'his'+str(prb+1)+'_log.'+formt
        print('Saving ',fname)
        plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust)
        
        
        # Phase space plot
        fig = plt.figure(figsize=(fig_height, fig_width))
        plt.xlabel(r'$v$')
        plt.ylabel(r'$\dot{v}$')
        plt.plot(v[2:-2],np.gradient(v)[2:-2]/dt,c='r',ls='-',lw=0.8)

        fname = 'his'+str(prb+1)+'_phase_space.'+formt
        print('Saving ',fname)
        plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust)
