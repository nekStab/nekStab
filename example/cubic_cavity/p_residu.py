#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
params = {'text.usetex': True,
          'font.size': 11,
          'legend.fontsize': 11,
          'legend.handlelength': 1.5,}
plt.rcParams.update(params)

formt = 'png'
ajust = 'tight'
qual = 300
fig_width = 3.2
fig_height = 2.3

class Res:
    def __init__(self, filename):
        print(f"Reading {filename}")
        data = np.transpose(np.genfromtxt(filename))
        if data.shape[0] == 2:
            self.i, self.r = data[0], data[1]
        if data.shape[0] == 3:
            self.i, self.k, self.r = data[0], data[1], data[2]
            
def plot_rs(ax, fname, lw = 1, c='gray', label=None):
    d = Res(fname)
    plt.plot(d.i,d.r,c=c,lw=lw,ls='-',label=label)
    return

########################################

if __name__ == '__main__':

    fig=plt.figure();fig.set_size_inches(fig_width, fig_height)
    plt.yscale('log')
    plt.ylabel(r'Residual')
    plt.xlabel(r'Iterations')
    plt.axhline(y=1e-09, lw=0.1, c='k', ls='dotted')
    
    plot_rs(plt,'residu.dat',lw=0.8,c='m',label=r'Case 1')

    plt.legend(loc='best')
    fname='residu.'+formt
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')

    fig=plt.figure();fig.set_size_inches(fig_width, fig_height)
    plt.yscale('log')
    plt.ylabel(r'Residual')
    plt.xlabel(r'Iteration')
    plt.axhline(y=1e-09, lw=0.1, c='k', ls='dotted')
    
    plot_rs(plt,'residu_arnoldi.dat', 0.8, 'k', r'Arnoldi')
    plot_rs(plt,'residu_newton.dat', 0.8, 'm', r'Newton')

    plt.legend(loc='best')
    fname='residu_newton.'+formt
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')
