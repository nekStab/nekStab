#!/usr/bin/env python
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
import os, os.path, csv, pickle
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

class time(object):
    def __init__(self, filename):
        print('Reading '+filename)
        data = np.transpose(np.genfromtxt(filename))
        self.k  = data[0]
        self.s  = data[1]
        self.t = data[2]
        del data

########################################

if __name__ == '__main__':

    fig=plt.figure();fig.set_size_inches(fig_width, fig_height)

    plt.yscale('log')
    plt.xlabel(r'$m\tau$')
    plt.ylabel(r'$min$')

    f = time('time.dat')
    for i in range(len(f.k)):
        print(f.k[i],f.s[i],f.k[i]*f.s[i],f.t[i]/60.)
        plt.scatter(f.k[i]*f.s[i],f.t[i]/60.)#,c='b',s='0.5', marker='o')
        plt.annotate('('+str(int(f.k[i]))+','+str(f.s[i])+')', (f.k[i]*f.s[i], f.t[i]/60.),fontsize=7)

    #plt.legend(loc='upper right',fontsize=6)
    fname='timerel.'+formt
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')
    

    fig=plt.figure();fig.set_size_inches(fig_width, fig_height)
    ax = fig.add_subplot(111,projection='3d')

    ax.set_zscale('log')
    ax.set_xlabel(r'$m$')
    ax.set_ylabel(r'$\tau$')
    ax.set_zlabel(r'$s$')

    for i in range(len(f.k)):
        ax.scatter(int(f.k[i]),f.s[i],f.t[i])#,c='b',s='0.5', marker='o')

    #plt.legend(loc='upper right',fontsize=6)
    fname='time.'+formt
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')
