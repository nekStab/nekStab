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
params = {'text.usetex': True,
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
        data = np.transpose(np.genfromtxt(filename))
        self.i  = data[0]
        self.r  = data[1]
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
    plt.yscale('log');plt.xlabel(r'$i$');plt.ylabel(r'$\epsilon$')
    #plt.xscale('log');
    plt.title(r'$2D-Cylinder,Re=80$',fontsize=8)
    plt.axhline(y=1e-11, lw=0.1, c='k', ls='dotted')

    f = res1('residu_arnoldi.dat')
    plt.plot(f.i,f.r,c='b',lw=0.0,ls='--', marker='o',markersize=0.4,label=r'ARNOLDI')
    for i in range(len(f.i)):
        print(i+1,f.r[i],f.i[i])
        if i> 1 and f.i[i] == 1.0:
            print()
            #plt.text(float(i), f.r[i], str(int(f.i[i])), color="k", fontsize=3)
            #plt.text(float(i-1), f.r[i-1], str(int(f.i[i-1])), color="k", fontsize=3)
        #else:
        #    plt.text(float(i+1), f.r[i], str(int(f.i[i])), color="b", fontsize=2)



    f = res1('residu_gmres.dat')
    plt.plot(f.i+1,f.r,c='g',lw=0.5,ls=':', marker='s',markersize=0.5,label=r'GMRES')
    try:
        for i in range(len(f.i)):
            print(i,f.r[i],f.i[i])
            if i> 1 and f.i[i] == 1.0:
                print()
                #plt.text(float(i), f.r[i], str(int(f.i[i])), color="k", fontsize=3)
                #plt.text(float(i-1), f.r[i-1], str(int(f.i[i-1])), color="k", fontsize=3)
            # plt.text(float(i), f.r[i], str(int(f.i[i])), color="k", fontsize=3)
            #else:
            #    plt.text(float(i+1), f.r[i], str(int(f.i[i])), color="g", fontsize=2)
    except:
        pass


    f = res1('residu_newton.dat')
    plt.plot(f.i,f.r,c='m',lw=0.5,ls='--', marker='d',markersize=0.8,label=r'NEWTON')
    #try:
    #    for i in range(len(f.i)):
    #        #print(i+1,f.r[i],f.i[i])
    #        #plt.text(float(i+1), f.r[i], str(int(f.i[i])), color="m", fontsize=2)
    #except:
    #    pass

    plt.legend(loc='upper right',fontsize=6)
    fname='residu_newton.'+formt
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')
