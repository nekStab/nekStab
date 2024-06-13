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
            self.r[self.r < 1e-12] = 1e-12

            # Filter out data points where the residual is larger than 0.1
            mask = self.r <= 0.1
            self.R = self.R[mask]
            self.I = self.I[mask]
            self.r = self.r[mask]
        else:
            self.r = None

def plot_H(ax, R, I, r, sized = 0, color='gray',symb = 'o', label=None):
    theta = np.linspace(0.0,2.*np.pi,400);ax.plot(np.cos(theta),np.sin(theta), lw=0.2, color='r', ls='-')
    
    # Normalize residuals to [0, 1]
    r_normalized = 1 - (10*r) #  - np.max(r)) / (np.max(r) - np.min(r))

    alpha_min = 0.1
    # # Scale to [alpha_min, 1]
    r_normalized = alpha_min + r_normalized * (1 - alpha_min)

    # Apply logarithmic transformation

    print('r_normalized:', r_normalized)
    
    # Create an empty list to store scatter plot objects
    scatters = []
    for k in range(len(R)):
        if symb == 'x':
            scatters.append(plt.scatter(R[k],I[k], s=5+sized, alpha=r_normalized[k],marker=symb,facecolors=color,linewidth=0.18))
        else:
            scatters.append(plt.scatter(R[k],I[k], s=5+sized, alpha=r_normalized[k],marker=symb,facecolors=color,edgecolors='k',linewidth=0.18))
            
    # Add label to the first scatter plot object
    scatters[0].set_label(label)
    
    ax.axhline(y=0., xmin=0, xmax=1, lw=0.2, color='k', ls=':')
    ax.axvline(x=0., ymin=0, ymax=1, lw=0.2, color='k', ls=':')
    return

def plot_NS(ax, R, I, r, freq=False, sized = 0, color='gray', symb = 'o', label=None):
    b = 1
    if freq:
        b = (2.*np.pi)
    
    # Normalize residuals to [0, 1]
    r_normalized =  (r - np.min(r)) / (np.max(r) - np.min(r))
    
    # Create an empty list to store scatter plot objects
    scatters = []
    for k in range(len(R)):
        if symb == 'x':
            scatters.append(ax.scatter(I[k]/b,R[k], s=5+sized, alpha=r_normalized[k],marker=symb,facecolors=color,linewidth=0.18))
        else:
            scatters.append(ax.scatter(I[k]/b,R[k], s=5+sized, alpha=r_normalized[k],marker=symb,facecolors=color,edgecolors='k',linewidth=0.18))
        if R[k] > 0:
            print('Mode: ',(k+1))
            print(' sigma=',R[k])
            print(' omega=',I[k])
            print('     f=',round(I[k]/b,7))
                
    # Add label to the first scatter plot object
    scatters[0].set_label(label)
    
    ax.axhline(y=0., xmin=0, xmax=1, lw=0.2, color='k', ls=':')
    ax.axvline(x=0., ymin=0, ymax=1, lw=0.2, color='k', ls=':')
    return

########################################

if __name__ == '__main__':

    fig=plt.figure();fig.set_size_inches(fig_width, fig_height)
    plt.axis('equal');xrz = [-1,0,1];xlabels = ['-1','0','1']
    plt.xticks(xrz,xlabels);plt.xlim(-1.1,1.1);plt.xlabel(r'$\Re (\mu)$')
    plt.yticks(xrz,xlabels);plt.ylim(-1.1,1.1);plt.ylabel(r'$\Im (\mu)$')

    f = Spectre('Spectre_Hd_reference.dat')
    plot_H(plt, f.R, f.I, f.r, 20, 'k', 'x', r'$Re=2000$')

    f = Spectre('Spectre_Hd.dat')
    plot_H(plt, f.R, f.I, f.r, 8, 'k', 'o', r'$Re=2000$')

    print("Most unstable Floquet multiplier:", f.R[0])

    reference_value = -1.17 # Shaabani 2019 JFM

    percentage_difference = ((f.R[0] - reference_value) / reference_value) * 100
    print("Percentage difference from the reference value: {:.2f}%".format(percentage_difference))

    fname='Spectre_H.'+formt
    plt.legend(loc='best',fontsize=6)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')

    fig=plt.figure();fig.set_size_inches(fig_width, fig_height)
    plt.ylabel(r'$\sigma$');#plt.xlim(-1,1)
    plt.xlabel(r'$f=\omega/2\pi$')#; plt.ylim(-0.4,0.3)

    f = Spectre('Spectre_NSd.dat')
    plot_NS(plt, f.R, f.I, f.r, True, 8, 'k', 'o', r'$Re=2000$')

    fname='Spectre_NS.'+formt
    plt.legend(loc='best',fontsize=6)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
    print('------------------------------------------')
