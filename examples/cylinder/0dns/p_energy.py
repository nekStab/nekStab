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
params = {'text.usetex': False,'font.size': 8,'legend.fontsize': 8,'legend.handlelength': 2.5,
'agg.path.chunksize':100000
};plt.rcParams.update(params)
transp = False
formt = 'png'
ajust = 'tight'
qual = 900
siz2 = 5.33
siz1 = 16.0*siz2/9.0

class ener(object):
    def __init__(self, filename):
        print('Reading '+filename)
        data = np.transpose(np.genfromtxt(filename))
        self.t  = data[0]
        self.u = data[1]
        self.v = data[2]
        self.w = data[3]
        self.et = data[4]
        del data

########################################

if __name__ == '__main__':

    fig=plt.figure();fig.set_size_inches((siz1, siz2))
    plt.xlabel(r'$t$')
    # plt.yscale('log')
    # plt.ylim(0.1915,0.21)
    # plt.xlim(0.,1000)

    f = ener('global_energy.dat');label=r'Re=50'
    plt.plot(f.t,f.u,c='k',lw=0.4,ls='-',label=label+' u')
    plt.plot(f.t,f.v,c='r',lw=0.4,ls='-',label=label+' v')

    plt.legend(loc='best',fontsize=6)
    fname='cube_global_energy.'+formt
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')
