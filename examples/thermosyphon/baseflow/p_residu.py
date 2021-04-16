#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
from numpy import pi
import scipy as sp
import os
import os.path
import csv
import pickle
params = {'text.usetex': False,
          'font.size': 8,
          'legend.fontsize': 8,
          'legend.handlelength': 2.5,}
plt.rcParams.update(params)
plt.style.use('seaborn-white')

formt = 'png'
ajust = 'tight'
qual = 500
fig_width = 3.5
fig_height = 2.45

class res(object):
    def __init__(self, filename):
        print('Reading '+filename)
        data = np.transpose(np.genfromtxt(filename))
        self.t  = data[0]
        self.r  = data[1]
        self.rt = data[2]
        del data

class res1(object):
    def __init__(self, filename):
        print('Reading '+filename)
        data = np.genfromtxt(filename)
        self.r = data
        del data

def plot_rs(ax, filename, sized = 0, color='gray', label=None):
    file = res(filename)
    plt.plot(file.t,file.r,c=color,lw=0.18+sized,ls='-',label=label)
    plt.axhline(y=file.r.min(), c=color,lw=0.18+sized,ls='--')
    plt.axvline(x=file.t.max(), c=color,lw=0.18+sized,ls='--')

    return

########################################

if __name__ == '__main__':

    fig=plt.figure();fig.set_size_inches(fig_width, fig_height)
    plt.yscale('log');plt.xlabel(r'$t$')
    plt.axhline(y=1e-02, lw=0.1, c='k', ls='dotted')
    plt.axhline(y=1e-04, lw=0.1, c='k', ls='dotted')
    plt.axhline(y=1e-06, lw=0.1, c='k', ls='dotted')
    plt.axhline(y=1e-08, lw=0.1, c='k', ls='dotted')
    plt.axhline(y=1e-10, lw=0.1, c='k', ls='dotted')
    plt.axhline(y=1e-12, lw=0.1, c='k', ls='dotted')
    plt.axhline(y=1e-13, lw=0.1, c='k', ls='dotted')

    plot_rs(plt, 'residu.dat',   0.2, 'r', r'Ra=400')

    plt.legend(loc='best',fontsize=6);
    fname='residu.'+formt
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')


    fig=plt.figure();fig.set_size_inches(fig_width, fig_height)
    plt.yscale('log');plt.xlabel(r'$nsteps$')
    plt.xscale('log')
    plt.axhline(y=1e-9, lw=0.1, c='k', ls='dotted')

    file = res1('residu_newton.dat')
    plt.plot(file.r,c='m',lw=0.5,label=r'NEWTON Ra=400')
    file = res1('residu_gmres.dat')
    plt.plot(file.r,c='g',lw=0.5,label=r'GMRES Ra=400')
    file = res1('residu_arnoldi.dat')
    plt.plot(file.r,c='b',lw=0.5,label=r'ARNOLDI Ra=400')

    plt.legend(loc='best',fontsize=6)
    fname='residu_newton.'+formt
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')