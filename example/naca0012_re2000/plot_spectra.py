#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
params = {'text.usetex': False,
          'font.size': 8,
          'legend.fontsize': 8,
          'legend.handlelength': 1.,}
plt.rcParams.update(params)

formt = 'png'
ajust = 'tight'
qual = 500
fig_width = 4.3
fig_height = 3.4

class Spectre(object):
    def __init__(self, filename):
        print('Reading '+filename)
        data = np.transpose(np.genfromtxt(filename))
        self.R = data[0]
        self.I = data[1]
        if data.shape[0] > 2:
            self.r = data[2]
        else:
            self.r = None

def plot_H(ax, R, I, r, sized = 0, color='gray',symb = 'o', label=None):
    iflabel = False
    theta = np.linspace(0.0,2.*np.pi,400);ax.plot(np.cos(theta),np.sin(theta), lw=0.2, color='r', ls='-')
    for k in range(len(R)):
        if r[k] < tolerance:
            mod = np.sqrt( (R[k])**2 + (I[k])**2 )
            if mod == 1:
                print('Time derivative of the baseflow found=',mod,R[k],I[k])
                plt.scatter(R[k],I[k], s=8+sized, alpha=0.8,marker='x',facecolors=color,edgecolors=color,linewidth=0.55)
            elif mod > 1:
                plt.scatter(R[k],I[k], s=5+sized, alpha=0.8,marker=symb,facecolors=color,edgecolors='k',linewidth=0.18)
        else:
            if iflabel == False:
                plt.scatter(R[k],I[k], s=5+sized, alpha=1,marker=symb,facecolors=color,edgecolors='k',linewidth=0.18,label=label)
                iflabel=True
            else:
                plt.scatter(R[k],I[k], s=5+sized, alpha=0.4,marker=symb,facecolors=color,edgecolors='k',linewidth=0.18)
            
    ax.axhline(y=0., xmin=0, xmax=1, lw=0.2, color='k', ls=':')
    ax.axvline(x=0., ymin=0, ymax=1, lw=0.2, color='k', ls=':')
    return

def plot_NS(ax, R, I, r, freq=False, sized = 0, color='gray', symb = 'o', label=None):
    b = 1; iflabel=False
    if freq:
        b = (2.*np.pi)
    for k in range(len(f.R)):
        if k == 0:
            ax.axhline(y=f.R[k], lw=0.2, color='r', ls=':')
            ax.axvline(x=f.I[k]/b, lw=0.2, color='r', ls=':')
        if r[k] > tolerance:
            ax.scatter(f.I[k]/b,f.R[k], s=5+sized, alpha=0.4,marker=symb,facecolors=color,edgecolors='k',linewidth=0.18)
        else:
            if f.R[k] > 0:
                print('Mode: ',(k+1))
                print(' sigma=',f.R[k])
                print(' omega=',f.I[k])
                print('     f=',round(f.I[k]/b,7))
                if iflabel == False:
                    ax.scatter(f.I[k]/b,f.R[k], alpha=1, s=7+sized,marker=symb, facecolors=color,edgecolors='k',linewidth=0.3,label=label)
                    iflabel = True
                else:
                    ax.scatter(f.I[k]/b,f.R[k], alpha=1, s=7+sized,marker=symb, facecolors=color,edgecolors='k',linewidth=0.3)
                
    ax.axhline(y=0., xmin=0, xmax=1, lw=0.2, color='k', ls=':')
    ax.axvline(x=0., ymin=0, ymax=1, lw=0.2, color='k', ls=':')
    return

########################################

tolerance = 1.0e-6
if __name__ == '__main__':

    fig=plt.figure();fig.set_size_inches(fig_width, fig_height)
    plt.axis('equal');xrz = [-1,0,1];xlabels = ['-1','0','1']
    plt.xticks(xrz,xlabels);plt.xlim(-1.1,1.1);plt.xlabel(r'$\Re (\mu)$')
    plt.yticks(xrz,xlabels);plt.ylim(-1.1,1.1);plt.ylabel(r'$\Im (\mu)$')

    f = Spectre('Spectre_Hd.dat')
    plot_H(plt, f.R, f.I, f.r, 8, 'k', 'o', r'$Re=1500$')

    fname='Spectre_H.'+formt
    plt.legend(loc='best',fontsize=6)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')

    fig=plt.figure();fig.set_size_inches(fig_width, fig_height)
    plt.ylabel(r'$\sigma$');#plt.xlim(-1,1)
    plt.xlabel(r'$f=\omega/2\pi$'); plt.ylim(-0.4,0.3)

    f = Spectre('Spectre_NSd.dat')
    plot_NS(plt, f.R, f.I, f.r, True, 8, 'k', 'o', r'$Re=2500$')

    fname='Spectre_NS.'+formt
    plt.legend(loc='best',fontsize=6)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')
