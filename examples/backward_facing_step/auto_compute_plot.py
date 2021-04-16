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
        #data = np.transpose(np.loadtxt(filename))
        print('Reading '+filename)
        data = np.transpose(np.genfromtxt(filename))
        self.t  = data[0]
        self.e  = data[1]
        del data

########################################

if __name__ == '__main__':

    fig=plt.figure();fig.set_size_inches(fig_width, fig_height)
    #plt.yscale('log')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$E$')

    f = res('ref_barkley2008_fig5.dat')
    plt.plot(f.t,f.e,c='k',lw=0.2,ls='-',label='Barkley (2008)')


    plt.legend(loc='best',fontsize=6);
    fname='envelope.'+formt
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')
