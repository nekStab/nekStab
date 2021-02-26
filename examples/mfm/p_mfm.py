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
          'legend.handlelength': 2.5}
plt.rcParams.update(params)
plt.style.use('seaborn-white')

#':' dotted,'-.' dashdot, '--' dashed, '-' solid
#b: blue, g: green, r: red,c: cyan,m: magenta,y: yellow,k: black

formt = 'png'
ajust = 'tight'
qual = 1500
fig_width = 3.5
fig_height = 2.45

class res(object):
    def __init__(self, filename):
        #data = np.transpose(np.loadtxt(filename))
        print('Reading '+filename)
        data = np.transpose(np.genfromtxt(filename))
        self.t  = data[0]
        self.r  = data[1]
        self.rt = data[2]
        del data
def plot_rs(ax, filename, sized = 0, color='gray', label=None):
    file = res(filename)
    plt.plot(file.t[1:-1],file.r[1:-1],c=color,lw=0.18+sized,ls='-',label=label)
    plt.axhline(y=file.r[-1], c=color,lw=0.1+sized,ls='--')
    #plt.axvline(x=file.t.max(), c=color,lw=0.18+sized,ls='--')
    return
def plot_rs2(ax, filename, sized = 0, color='gray', label=None):
    file = res(filename)
    plt.plot(file.t[1:-1],file.rt[1:-1],c=color,lw=0.18+sized,ls='-',label=label)
    plt.axhline(y=file.rt[-1], c=color,lw=0.1+sized,ls='--')
    #plt.axvline(x=file.t.max(), c=color,lw=0.18+sized,ls='--')
            
    return

########################################

if __name__ == '__main__':
    fig=plt.figure();fig.set_size_inches(fig_width, fig_height)
    plt.ylabel(r'$A$');#plt.ylim(-0.1,0.005)#plt.yscale('log');
    plt.xlabel(r'$t$')#;plt.xlim(-0.2,0.2)

    #plt.axhline(y=1e-02, lw=0.1, c='k', ls='dotted')

    plot_rs2(plt, 'mfm.dat',   0.2, 'm', r'Re=1000')

    plt.legend(loc='best',fontsize=6);
    fname='mfm_A.'+formt
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')
    fig=plt.figure();fig.set_size_inches(fig_width, fig_height)
    plt.ylabel(r'$K$');#plt.ylim(-0.1,0.005)
    plt.yscale('log');
    plt.xlabel(r'$t$')#;plt.xlim(-0.2,0.2)
    #plt.axhline(y=1e-02, lw=0.1, c='k', ls='dotted')

    plot_rs(plt, 'mfm.dat',   0.2, 'm', r'Re=1000')

    plt.legend(loc='best',fontsize=6);
    fname='mfm_K.'+formt
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')
