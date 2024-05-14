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
        num_columns = data.shape[0]

        if num_columns == 2:
            self.R = data[0]
            self.I = data[1]
            self.r = None
        elif num_columns == 3:
            self.R = data[0]
            self.I = data[1]
            self.r = data[2]
        else:
            raise ValueError(f"Expected 2 or 3 columns in file, but got {num_columns}")

def plot_H(ax, R, I, r, sized = 0, color='gray',symb = 'o', label=None):
    iflabel = False
    theta = np.linspace(0.0,2.*np.pi,400);ax.plot(np.cos(theta),np.sin(theta), lw=0.2, color='r', ls='-')
    for k in range(len(R)):
        if r[k] < tolerance:
            plt.scatter(R[k],I[k], s=5+sized, alpha=1,marker=symb,facecolors=color,edgecolors='k',linewidth=0.18)
        else:
            if iflabel == False:
                plt.scatter(R[k],I[k], s=5+sized, alpha=1,marker=symb,facecolors=color,edgecolors='k',linewidth=0.18,label=label)
                iflabel=True
            else:
                plt.scatter(R[k],I[k], s=5+sized, alpha=1,marker=symb,facecolors=color,edgecolors='k',linewidth=0.18)
            
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

def plot_S(ax, R, I, sized = 0, color='gray',symb = 'o', label=None):
    iflabel = False
    for k in range(len(R)):
        #if I[k] < tolerance:
        if iflabel == False:
            plt.scatter(k+1,R[k], s=5+sized, alpha=1,marker=symb,facecolors=color,edgecolors='k',linewidth=0.18,label=label)
            iflabel=True
        else:
            plt.scatter(k+1,R[k], s=5+sized, alpha=1,marker=symb,facecolors=color,edgecolors='k',linewidth=0.18)
    return
########################################

tolerance = 1.0e-6
if __name__ == '__main__':

    spectre_files_H = [
        ('/home/rfrantz/nekStab/examples/cylinder/stability/direct/Spectre_Hd.dat', 22, 'k', 'o', r'$Re=50$ (legacy)'),
        ('Spectrum_Hd.dat', 12, 'g', 'o', r'$Re=50$'),
        ('Spectrum_Hf.dat', 12, 'm', 'd', r'$Re=50$ FD=2'),
        ('Spectrum_Ha.dat', 20, 'r', '+', r'$Re=50^{\dagger}$')
    ]

    spectre_files_NS = [
        ('/home/rfrantz/nekStab/examples/cylinder/stability/direct/Spectre_NSd.dat', 22, 'k', 'o', r'$Re=50$ (legacy)'),
        ('Spectrum_NSd.dat', 12, 'g', 'o', r'$Re=50$'),
        ('Spectrum_NSf.dat', 12, 'm', 'd', r'$Re=50$ FD=2'),
        ('Spectrum_NSa.dat', 20, 'r', '+', r'$Re=50^{\dagger}$')
    ]

    spectre_files_S = [
        ('Spectrum_Sp.dat', 12, 'g', 'o', r'$Re=50$ (transient growth)'),
        ('Spectrum_Sr.dat', 12, 'r', 'o', r'$Re=50$ (resolvent)'),
    ]

    fig=plt.figure();fig.set_size_inches(fig_width, fig_height)
    plt.axis('equal');xrz = [-1,0,1];xlabels = ['-1','0','1']
    plt.xticks(xrz,xlabels);plt.xlim(-1.1,1.1);plt.xlabel(r'$\Re (\mu)$')
    plt.yticks(xrz,xlabels);plt.ylim(-1.1,1.1);plt.ylabel(r'$\Im (\mu)$')
    for file, size, color, marker, label in spectre_files_H:
        try:
            f = Spectre(file)
            plot_H(plt, f.R, f.I, f.r, size, color, marker, label)
        except Exception as e:
            print(f"Error processing file {file}: {e}")
    fname='Spectre_H.'+formt
    plt.legend(loc='best',fontsize=6)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')

    fig=plt.figure();fig.set_size_inches(fig_width, fig_height)
    plt.ylabel(r'$\sigma$');#plt.xlim(-1,1)
    plt.xlabel(r'$f=\omega/2\pi$'); plt.ylim(-0.4,0.3)
    for file, size, color, marker, label in spectre_files_NS:
        try:
            f = Spectre(file)
            plot_NS(plt, f.R, f.I, f.r, True, size, color, marker, label)
        except Exception as e:
            print(f"Error processing file {file}: {e}")
    fname='Spectre_NS.'+formt
    plt.legend(loc='best',fontsize=6)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')



    fig=plt.figure();fig.set_size_inches(fig_width, fig_height)
    plt.ylabel(r'$\sigma$');#plt.xlim(-1,1)
    plt.yscale('log');plt.xscale('log')
    plt.xlabel(r'$i$')#; plt.ylim(-0.4,0.3)
    for file, size, color, marker, label in spectre_files_S:
        try:
            f = Spectre(file)
            plot_S(plt, f.R, f.I, size, color, marker, label)
        except Exception as e:
            print(f"Error processing file {file}: {e}")

    fname='Spectre_S.'+formt
    plt.legend(loc='best',fontsize=6)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')
