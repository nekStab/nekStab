#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
formt = 'png'
ajust = 'tight'
qual = 500
fig_width = 4.3
fig_height = 3.4
        
class residual_reader:
    def __init__(self, filename):
        print(f"Reading {filename}")
        self.calls, self.ttime, self.r = np.transpose(np.genfromtxt(filename))

if __name__ == '__main__':

    fig=plt.figure()
    fig.set_size_inches(fig_width, fig_height)
    
    plt.xscale('log')
    plt.yscale('log')
    
    plt.xlabel(r'linearized calls')
    plt.ylabel(r'residual')
    
    try:
        d = residual_reader('residu.dat')
        plt.plot(d.calls,d.r,c='m',lw=2,ls='-', marker='o',markersize=1)
    except:
        print('First Newton iteartion not finished!')
        pass
    
    fname='residu_newton.'+formt
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust)
    print('Saving '+fname)
    plt.close()
    print('------------------------------------------')