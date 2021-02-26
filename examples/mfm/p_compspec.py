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
          'legend.handlelength': 1.,
          'font.family' : 'lmodern'}
plt.rcParams['text.latex.preamble']=[r'\usepackage{lmodern}']
plt.rcParams.update(params)
plt.style.use('seaborn-white')

#':' dotted,'-.' dashdot, '--' dashed, '-' solid
#b: blue, g: green, r: red,c: cyan,m: magenta,y: yellow,k: black

formt = 'png'
ajust = 'tight'
qual = 900
fig_width = 3.5 #5.33
fig_height = 2.45 #5.33/3

class SpectreH(object):
    def __init__(self, filename):
        #data = np.transpose(np.loadtxt(filename))
        print('Reading '+filename)
        data = np.transpose(np.genfromtxt(filename))
        self.vp_real  = data[0]
        self.vp_imag  = data[1]
        self.residu   = data[2]
        del data

class SpectreNS(object):
     def __init__(self, filename):
        print('Reading '+filename)
        #data = np.transpose(np.loadtxt(filename))
        data = np.transpose(np.genfromtxt(filename))
        self.vp_real  = data[0]
        self.vp_imag  = data[1]
        del data

def plot_H(ax, vp_real, vp_imag, residu, sized = 0, color='gray', label=None):
    iflabel = False
    theta = np.linspace(0.0,2.*np.pi,400);ax.plot(np.cos(theta),np.sin(theta), lw=0.2, color='r', ls='-')
    for k in range(len(vp_real)):
        if residu[k] < tolerance:
            mod = np.sqrt( (vp_real[k])**2 + (vp_imag[k])**2 )
            if mod == 1:
                print('Time derivative of the baseflow found=',mod,vp_real[k],vp_imag[k])
                plt.scatter(vp_real[k],vp_imag[k], s=8+sized, alpha=0.8,marker='x',facecolors=color,edgecolors=color,linewidth=0.55)
            elif mod > 1:
                print('Mode: ',(k+1),'|',mod,'|')
                print('   Re=',vp_real[k])
                print('angle: ',np.angle(vp_imag[k]))
                print('    f=',(vp_imag[k]/(2.*np.pi)))
                print('--------------------------------------------------------------------')
                if iflabel == False:
                    plt.scatter(vp_real[k],vp_imag[k], s=5+sized, alpha=0.8,marker='o',facecolors='none',edgecolors='k',linewidth=0.18,label=label)
                    iflabel=True
                else:
                    plt.scatter(vp_real[k],vp_imag[k], s=5+sized, alpha=0.8,marker='o',facecolors='none',edgecolors='k',linewidth=0.18)
            else:
                if iflabel == False:
                    plt.scatter(vp_real[k],vp_imag[k], s=5+sized, alpha=0.8,marker='o',facecolors=color,edgecolors='k',linewidth=0.18,label=label)
                    iflabel=True
                else:
                    plt.scatter(vp_real[k],vp_imag[k], s=5+sized, alpha=0.8,marker='o',facecolors=color,edgecolors='k',linewidth=0.18)
        else:
            plt.scatter(vp_real[k],vp_imag[k], s=1, alpha=0.5,color=color)
    ax.axhline(y=0., xmin=0, xmax=1, lw=0.2, color='k', ls=':')
    ax.axvline(x=0., ymin=0, ymax=1, lw=0.2, color='k', ls=':')
    return

def plot_NS(ax, vp_real, vp_imag, freq=False, sized = 0, color='gray',label=None):
    b = 1; iflabel=False;
    if freq:
        b = (2.*np.pi)
    for k in range(len(file.vp_real)):
        if file.vp_real[k] > 0:
            print('Mode: ',(k+1),file.vp_real[k],file.vp_imag[k])
            print('      f=',file.vp_imag[k]/b)
            ax.scatter(file.vp_imag[k]/b,file.vp_real[k], alpha=0.8, s=7+sized,marker='o', facecolors=color,edgecolors='k',linewidth=0.3)
        else:
            if iflabel == False:
                ax.scatter(file.vp_imag[k]/b,file.vp_real[k], alpha=0.8, s=5+sized,marker='o',facecolors=color,edgecolors='k',linewidth=0.2,label=label)
                iflabel = True
            else:
                ax.scatter(file.vp_imag[k]/b,file.vp_real[k], alpha=0.8, s=5+sized,marker='o',facecolors=color,edgecolors='k',linewidth=0.2)
    ax.axhline(y=0., xmin=0, xmax=1, lw=0.2, color='k', ls=':')
    ax.axvline(x=0., ymin=0, ymax=1, lw=0.2, color='k', ls=':')
    return

########################################

tolerance = 1.0e-6


if __name__ == '__main__':

    fig=plt.figure();fig.set_size_inches(fig_height, fig_height)
    plt.axis('equal');xrz = [-1,0,1];xlabels = ['-1','0','1']
    plt.xticks(xrz,xlabels);plt.xlim(-1.5,1.5);plt.xlabel(r'$\Re (\mu)$')
    plt.yticks(xrz,xlabels);plt.ylim(-1.5,1.5);plt.ylabel(r'$\Im (\mu)$')

    file = SpectreH('Spectre_H.dat')
    plot_H(plt, file.vp_real, file.vp_imag, file.residu, 0, 'r', r'Re=400')

    fname='Spectre_H.'+formt
    plt.legend(loc='best',fontsize=6);
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')

    fig=plt.figure();fig.set_size_inches(fig_width, fig_height)
    plt.xlabel(r'$\Re (\lambda)$');#plt.xlim(-1,1)
    plt.ylabel(r'$\Im (\lambda)$');#plt.ylim(-0.2,0.2)

    # plt.axvline(x=0, lw=0.4, color='g', ls='dashdot',label='S1')

    file = SpectreNS('Spectre_NS_conv.dat')
    plot_NS(plt, file.vp_real, file.vp_imag, False, 0, 'r', r'Re=400')


    fname='Spectre_NS.'+formt
    plt.legend(loc='best',fontsize=4);
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')


    fig=plt.figure();fig.set_size_inches(fig_width, fig_height)
    plt.ylabel(r'$\sigma$');#plt.ylim(-0.1,0.005)
    plt.xlabel(r'$f$');#plt.xlim(-0.2,0.2)

    # plt.axvline(x=0, lw=0.4, color='g', ls='dashdot',label='S1')
    # plt.axvline(x=0.12, lw=0.4, color='m', ls='dashed',label='T1')

    file = SpectreNS('Spectre_NS_conv.dat')
    plot_NS(plt, file.vp_real, file.vp_imag, True, 0, 'r', r'Re=400')

    plt.legend(loc='best',fontsize=4);
    fname='Spectre_NSf.'+formt
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')
