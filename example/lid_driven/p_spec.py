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

class Spectre:
    def __init__(self, filename):
        print(f"Reading {filename}")
        data = np.transpose(np.genfromtxt(filename))
        self.re, self.im = data[0], data[1]
        if data.shape[0] == 3:
            self.res = data[2]

def plot_H(ax, re, im, r, sized = 0, color='gray',symb = 'o', label=None):
    iflabel = False
    theta = np.linspace(0.0,2.*np.pi,400);ax.plot(np.cos(theta),np.sin(theta), lw=0.2, color='r', ls='-')
    for k in range(len(re)):
        if r[k] < tolerance:
            mod = np.sqrt( (re[k])**2 + (im[k])**2 )
            if mod == 1:
                print('Time derivative of the baseflow found=',mod,re[k],im[k])
                plt.scatter(re[k],im[k], s=8+sized, alpha=0.8,marker='x',facecolors=color,edgecolors=color,lw=0.55)
            elif mod > 1:
                plt.scatter(re[k],im[k], s=5+sized, alpha=0.8,marker=symb,facecolors=color,edgecolors='k',lw=0.18)
        else:
            if iflabel == False:
                plt.scatter(re[k],im[k], s=5+sized, alpha=1,marker=symb,facecolors=color,edgecolors='k',lw=0.18,label=label)
                iflabel=True
            else:
                plt.scatter(re[k],im[k], s=5+sized, alpha=0.4,marker=symb,facecolors=color,edgecolors='k',lw=0.18)
            
    ax.axhline(y=0, lw=0.2, c='k', ls=':')
    ax.axvline(x=0, lw=0.2, c='k', ls=':')
    return

def plot_NS(ax, re, im, freq=False, sized = 0, color='gray', symb = 'o', label=None):
    b = 1; iflabel=False
    if freq:
        b = (2.*np.pi)
    for k in range(len(f.re)):
        if f.re[k] > 0:
            print('Mode: ',(k+1),f.re[k],f.im[k])
            print('      f=',f.im[k]/b)
            ax.scatter(f.im[k]/b,f.re[k], alpha=1, s=7+sized,marker=symb, facecolors=color,edgecolors='k',lw=0.3)
        else:
            if iflabel == False:
                ax.scatter(f.im[k]/b,f.re[k], s=5+sized, alpha=0.8,marker=symb,facecolors=color,edgecolors='k',lw=0.18,label=label)
                iflabel = True
            else:
                ax.scatter(f.im[k]/b,f.re[k], s=5+sized, alpha=0.8,marker=symb,facecolors=color,edgecolors='k',lw=0.18)
                
    ax.axhline(y=0., lw=0.2, c='k', ls=':')
    ax.axvline(x=0., lw=0.2, c='k', ls=':')
    return

########################################

tolerance = 1.0e-6
if __name__ == '__main__':

    fig=plt.figure();fig.set_size_inches(fig_width, fig_height)
    plt.ylabel(r'$\Re (\lambda)$');#plt.xlim(-1,1)
    plt.xlabel(r'$\Im (\lambda)=\omega$');# plt.ylim(-0.1,0.02)
    
    f = Spectre('Spectre_NSd_conv.dat')
    plot_NS(plt, f.re, f.im, False, 8, 'k', 'o')
 
    fname='Spectre_NS.'+formt
    #plt.legend(loc='best')
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')
    
    
    fig=plt.figure();fig.set_size_inches(fig_width, fig_height)
    plt.ylabel(r'$\Re (\lambda)$');#plt.xlim(-1,1)
    plt.xlabel(r'$f=\omega/2\pi$');# plt.ylim(-0.1,0.02)
    
    f = Spectre('Spectre_NSd_conv.dat')
    plot_NS(plt, f.re, f.im, True, 8, 'k', 'o',label=r'$Re=3600$')
    plt.axvline(x=0.34, lw=0.9, c='r', ls=':',label=r'$f=0.34$')
 
    fname='Spectre_NSf.'+formt
    plt.legend(loc='best')
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')
    
    
    fig=plt.figure();fig.set_size_inches(fig_height, fig_height)
    plt.axis('equal');xrz = [-1,0,1];xlabels = ['-1','0','1']
    plt.xticks(xrz,xlabels);plt.xlim(-1.1,1.1);plt.xlabel(r'$\Re (\mu)$')
    plt.yticks(xrz,xlabels);plt.ylim(-1.1,1.1);plt.ylabel(r'$\Im (\mu)$')

    f = Spectre('Spectre_Hd.dat')
    plot_H(plt, f.re, f.im, f.res, 8, 'k', 'o', r'$Re=50$')

    fname='Spectre_H.'+formt
    plt.legend(loc='best')
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')
