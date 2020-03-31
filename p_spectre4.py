import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import scipy as sp
import os
import os.path
import csv
import pickle
params = {'text.usetex': False,
          'font.size': 10,
          'legend.fontsize': 8,
          'legend.handlelength': 2.5}
plt.rcParams.update(params)
#':' dotted,'-.' dashdot, '--' dashed, '-' solid
#b: blue, g: green, r: red,c: cyan,m: magenta,y: yellow,k: black

transp = False
formt = 'png'
ajust = 'tight'
qual = 1500
siz1 = 3.5
siz2 = 2.45

class Spectre(object):
    def __init__(self, filename):
        #data = np.transpose(np.loadtxt(filename))
        data = np.transpose(np.genfromtxt(filename))
        self.vp_real  = data[0]
        self.vp_imag  = data[1]
        self.residu   = data[2]
        del data

class Spectre_c(object):
    def __init__(self, filename):
        #data = np.transpose(np.loadtxt(filename))
        data = np.transpose(np.genfromtxt(filename))
        self.vp_real  = data[0]
        self.vp_imag  = data[1]
        del data


########################################

tolerance = 1.0e-6

fname='Spectre_H.'+formt
#file = Spectre('H0098')
file = Spectre('Spectre_H.dat')
print('Reading Spectre_H.dat')

fig=plt.figure()
fig.set_size_inches(siz2, siz2)
plt.axis('equal')
plt.xlabel(r'$\Re (\mu)$')
plt.ylabel(r'$\Im (\mu)$')
plt.xlim(-1.5,1.5)
plt.ylim(-1.5,1.5)
xrz = [-1,0,1]
xlabels = ['-1','0','1']
plt.xticks(xrz,xlabels)
plt.yticks(xrz,xlabels)

plt.axhline(y=0., xmin=0, xmax=1, lw=0.2, color='k', ls=':')
plt.axvline(x=0., ymin=0, ymax=1, lw=0.2, color='k', ls=':')

theta = np.linspace(0.0,2.*np.pi,400)
plt.plot(np.cos(theta),np.sin(theta), lw=0.2, color='r', ls='-')

for k in range(len(file.vp_real)):
    if file.residu[k] < tolerance:
        mod = np.sqrt( (file.vp_real[k])**2 + (file.vp_imag[k])**2 )
        if mod == 1:
            print('Time derivative of the baseflow found=',mod,file.vp_real[k],file.vp_imag[k])
            plt.scatter(file.vp_real[k],file.vp_imag[k], s=8,marker='x',facecolors='gray',edgecolors='gray',linewidth=0.6)
        elif mod > 1:
            print('Mode: ',(k+1),'|',mod,'|')
            print('   Re=',file.vp_real[k])
            print('angle: ',np.angle(file.vp_imag[k]))
            print('    f=',(file.vp_imag[k]/(2.*np.pi)))
            print('--------------------------------------------------------------------')
            plt.scatter(file.vp_real[k],file.vp_imag[k], s=5,marker='o',facecolors='none',edgecolors='k',linewidth=0.2)
        else:
            plt.scatter(file.vp_real[k],file.vp_imag[k], s=5,marker='o',facecolors='gray',edgecolors='k',linewidth=0.2)
    else:
        plt.scatter(file.vp_real[k],file.vp_imag[k], s=1, alpha=0.5,color='gray')

#plt.legend(loc='best',fontsize=6);
plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()
print('------------------------------------------')

fname='Spectre_NS.'+formt
fig=plt.figure()
fig.set_size_inches(siz1, siz2)
#plt.xlabel(r'$\sigma$')
#plt.ylabel(r'$\omega$')
plt.xlabel(r'$\Re (\lambda)$')
plt.ylabel(r'$\Im (\lambda)$')
#plt.ylim(-0.25,0.05)
#plt.xlim(-1.,1.)

plt.axhline(y=0., xmin=0, xmax=1, lw=0.2, color='k', ls=':')
plt.axvline(x=0., ymin=0, ymax=1, lw=0.2, color='k', ls=':')

print('Reading Spectre_NS.dat')
file = Spectre('Spectre_NS.dat')
for k in range(len(file.vp_real)):
         if file.vp_real[k] > 0:
             print('Mode: ',(k+1),file.vp_real[k],file.vp_imag[k],file.vp_imag[k]/(2.*np.pi))
             plt.scatter(file.vp_real[k],file.vp_imag[k], s=5,marker='o',facecolors='none',edgecolors='k',linewidth=0.2)
         else:
             #print('Mode: ',(k+1),file.vp_real[k],file.vp_imag[k])
             plt.scatter(file.vp_real[k],file.vp_imag[k], s=5,marker='o',facecolors='gray',edgecolors='k',linewidth=0.2)

#plt.legend(loc='best',fontsize=6);
plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()


fig=plt.figure()
fig.set_size_inches(siz1, siz2)
plt.ylabel(r'$\sigma$')
plt.xlabel(r'$f$')

#plt.ylim(-0.005,0.002)
#plt.xlim(0,0.009)

plt.axhline(y=0., xmin=0, xmax=1, lw=0.2, color='k', ls=':')
plt.axvline(x=0., ymin=0, ymax=1, lw=0.2, color='k', ls=':')
#plt.axvline(x=0.00821362, ymin=0, ymax=1, lw=0.6, color='r', ls=':',label='f')
#plt.axvline(x=0.00821362/2, ymin=0, ymax=1, lw=0.4, color='g', ls=':',label='f/2')
#plt.axvline(x=0.00821362/4., ymin=0, ymax=1, lw=0.2, color='b', ls=':',label='f/4')

print('Reading Spectre_NS.dat')
file = Spectre('Spectre_NS.dat')
for k in range(len(file.vp_real)):
         if file.vp_real[k] > 0:
             print('Mode: ',(k+1),file.vp_real[k],file.vp_imag[k])
             print('      f=',file.vp_imag[k]/(2.*np.pi))
             plt.scatter(file.vp_imag[k]/(2.*np.pi),file.vp_real[k], s=5,marker='o',facecolors='none',edgecolors='k',linewidth=0.2)
         else:
             #print('Mode: ',(k+1),file.vp_real[k],file.vp_imag[k])
             plt.scatter(file.vp_imag[k]/(2.*np.pi),file.vp_real[k], s=5,marker='o',facecolors='gray',edgecolors='k',linewidth=0.2)

plt.legend(loc='best',fontsize=6);
fname = 'Spectre_NS2f'
plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust);print('Saving '+fname);plt.close()


'''
for i in range(2, 300, 1):
    try:
        fname = 'H'+str(i).zfill(4)
        print('Reading '+fname)
        file = Spectre(fname)
        fig=plt.figure()
        fig.set_size_inches(siz2, siz2)
        plt.axis('equal')
        plt.xlabel(r'$\Re (\mu)$')
        plt.ylabel(r'$\Im (\mu)$')
        xrz = [-1,0,1]
        xlabels = ['-1','0','1']
        plt.xticks(xrz,xlabels)
        plt.yticks(xrz,xlabels)
        plt.axhline(y=0., xmin=0, xmax=1, lw=0.2, color='k', ls=':')
        plt.axvline(x=0., ymin=0, ymax=1, lw=0.2, color='k', ls=':')
        theta = np.linspace(0.0,2.*np.pi,400)
        plt.plot(np.cos(theta),np.sin(theta), lw=0.2, color='r', ls='-')
        saveyes = 0
        for k in range(i):
            if file.residu[k] < tolerance:
                saveyes = 1
                mod = np.sqrt( file.vp_real[k]**2 + file.vp_imag[k]**2 )
                if mod == 1:
                    #print('Time derivative of the baseflow found=',mod,file.vp_real[k],file.vp_imag[k])
                    plt.scatter(file.vp_real[k],file.vp_imag[k], s=8,marker='x',facecolors='gray',edgecolors='gray',linewidth=0.6)
                elif mod > 1:
                    #print('Mode: ',(k+1), mod, file.vp_real[k],file.vp_imag[k])
                    plt.scatter(file.vp_real[k],file.vp_imag[k], s=5,marker='o',facecolors='none',edgecolors='k',linewidth=0.2)
                else:
                    plt.scatter(file.vp_real[k],file.vp_imag[k], s=5,marker='o',facecolors='gray',edgecolors='k',linewidth=0.2)
            else:
                plt.scatter(file.vp_real[k],file.vp_imag[k], s=0, alpha=0,color='gray')
        #plt.legend(loc='best',fontsize=6);
        if saveyes == 1:
            print('Saving '+fname+'.'+formt)
            plt.savefig(fname+'.'+formt,format=formt,dpi=qual,bbox_inches=ajust);plt.close()
    except:
        #print('Skipping '+fname)
        pass
'''
