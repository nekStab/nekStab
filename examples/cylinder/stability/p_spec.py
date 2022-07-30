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
          'legend.handlelength': 1.,}
plt.rcParams.update(params)
plt.style.use('seaborn-white')

#':' dotted,'-.' dashdot, '--' dashed, '-' solid
#b: blue, g: green, r: red,c: cyan,m: magenta,y: yellow,k: black

formt = 'png'
ajust = 'tight'
qual = 500
fig_width = 3.5
fig_height = 2.45

class SpectreH(object):
    def __init__(self, filename):
        #data = np.transpose(np.loadtxt(filename))
        print('Reading '+filename)
        data = np.transpose(np.genfromtxt(filename))
        self.vpr = data[0]
        self.vpi = data[1]
        self.r   = data[2]
        del data

class SpectreNS(object):
     def __init__(self, filename):
        print('Reading '+filename)
        #data = np.transpose(np.loadtxt(filename))
        data = np.transpose(np.genfromtxt(filename))
        self.vpr  = data[0]
        self.vpi  = data[1]
        del data

def plot_H(ax, vpr, vpi, r, sized = 0, color='gray',symb = 'o', label=None):
    iflabel = False
    theta = np.linspace(0.0,2.*np.pi,400);ax.plot(np.cos(theta),np.sin(theta), lw=0.2, color='r', ls='-')
    for k in range(len(vpr)):
        if r[k] < tolerance:
            mod = np.sqrt( (vpr[k])**2 + (vpi[k])**2 )
            if mod == 1:
                print('Time derivative of the baseflow found=',mod,vpr[k],vpi[k])
                plt.scatter(vpr[k],vpi[k], s=8+sized, alpha=0.8,marker='x',facecolors=color,edgecolors=color,linewidth=0.55)
            elif mod > 1:
                plt.scatter(vpr[k],vpi[k], s=5+sized, alpha=0.8,marker=symb,facecolors=color,edgecolors='k',linewidth=0.18)
        else:
            if iflabel == False:
                plt.scatter(vpr[k],vpi[k], s=5+sized, alpha=1,marker=symb,facecolors=color,edgecolors='k',linewidth=0.18,label=label)
                iflabel=True
            else:
                plt.scatter(vpr[k],vpi[k], s=5+sized, alpha=0.4,marker=symb,facecolors=color,edgecolors='k',linewidth=0.18)
            
    ax.axhline(y=0., xmin=0, xmax=1, lw=0.2, color='k', ls=':')
    ax.axvline(x=0., ymin=0, ymax=1, lw=0.2, color='k', ls=':')
    return

def plot_NS(ax, vpr, vpi, freq=False, sized = 0, color='gray', symb = 'o', label=None):
    b = 1; iflabel=False
    if freq:
        b = (2.*np.pi)
    for k in range(len(f.vpr)):
        if f.vpr[k] > 0:
            print('Mode: ',(k+1),f.vpr[k],f.vpi[k])
            print('      f=',f.vpi[k]/b)
            ax.scatter(f.vpi[k]/b,f.vpr[k], alpha=1, s=7+sized,marker=symb, facecolors=color,edgecolors='k',linewidth=0.3)
        else:
            if iflabel == False:
                ax.scatter(f.vpi[k]/b,f.vpr[k], s=5+sized, alpha=0.8,marker=symb,facecolors=color,edgecolors='k',linewidth=0.18,label=label)
                iflabel = True
            else:
                ax.scatter(f.vpi[k]/b,f.vpr[k], s=5+sized, alpha=0.8,marker=symb,facecolors=color,edgecolors='k',linewidth=0.18)
                
    ax.axhline(y=0., xmin=0, xmax=1, lw=0.2, color='k', ls=':')
    ax.axvline(x=0., ymin=0, ymax=1, lw=0.2, color='k', ls=':')
    return

########################################

tolerance = 1.0e-6
if __name__ == '__main__':

    fig=plt.figure();fig.set_size_inches(fig_height, fig_height)
    plt.axis('equal');xrz = [-1,0,1];xlabels = ['-1','0','1']
    plt.xticks(xrz,xlabels);plt.xlim(-1.1,1.1);plt.xlabel(r'$\Re (\mu)$')
    plt.yticks(xrz,xlabels);plt.ylim(-1.1,1.1);plt.ylabel(r'$\Im (\mu)$')

    f = SpectreH('direct/Spectre_Hd.dat')
    plot_H(plt, f.vpr, f.vpi, f.r, 8, 'k', 'o', r'$Re=50$')
    f = SpectreH('adjoint/Spectre_Ha.dat')
    plot_H(plt, f.vpr, f.vpi, f.r, 2, 'r', 'o', r'$Re=50^{\dagger}$')
    
    fname='Spectre_H.'+formt
    plt.legend(loc='best',fontsize=6)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')

    fig=plt.figure();fig.set_size_inches(fig_width, fig_height)
    plt.xlabel(r'$\Re (\lambda)$');#plt.xlim(-1,1)
    plt.ylabel(r'$\Im (\lambda)$');# plt.ylim(-0.1,0.02)

    plt.scatter(0.73562,0.01264,s=4,marker='*',label=r'Marquet et al (2019)')
    
    f = SpectreNS('direct/Spectre_NSd_conv.dat')
    plot_NS(plt, f.vpr, f.vpi, False, 8, 'k', 'o', r'$Re=50$')
    f = SpectreNS('adjoint/Spectre_NSa_conv.dat')
    plot_NS(plt, f.vpr, f.vpi, False, 2, 'r', 'o', r'$Re=50^{\dagger}$')
    
    fname='Spectre_NS.'+formt
    plt.legend(loc='best',fontsize=6)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')
