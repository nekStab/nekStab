#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.colors as mcolors
import numpy as np
from cycler import cycler
from scipy.optimize import curve_fit
from scipy import interpolate
from numpy import pi, interp
import scipy as sp
from math import isclose
import os
import os.path
import csv
import pickle

params = {
    "text.usetex": False,
    "legend.handlelength": 1.0,
}
plt.rcParams.update(params)
formt = "png"
ajust = "tight"
qual = 300
fig_width = 3.45
fig_height = 3

class file(object):
    def __init__(self, filename):
        # data = np.transpose(np.loadtxt(filename))
        print("Reading " + filename)
        data = np.transpose(np.genfromtxt(filename))
        self.x = data[0]
        self.y = data[1]
        del data

if __name__ == "__main__":


    fig = plt.figure()
    fig.set_size_inches(fig_width, fig_height)
    plt.ylabel(r"$St$")
    plt.xlabel(r"$\alpha_c$")
    plt.xlim(0,100); #plt.ylim(0,0.8)

    f = file("naca4412_1.csv")
    poly = np.poly1d(np.polyfit(f.x,f.y, 9))
    print(poly)
    plt.scatter(f.x, f.y, c='k',s=8)
    # def St(x):
    #     return  -6.192e-20*x**12  + 3.422e-17*x**11  - 8.226e-15*x**10  + 1.122e-12*x**9 - 9.456e-11*x**8 + 4.952e-09*x**7 - 1.455e-07*x**6 + 9.969e-07*x**5 + 9.282e-05*x**4 - 0.003851*x**3 + 0.07199*x**2 - 0.7218*x**1 + 3.675

    def func(x, a0,a1,a2,a3,a4,a5,a6,a7,a8,a9):
        #return a * np.exp(-b * x) + c
        return a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4 + a5*x**5 + a6*x**6 + a7*x**7 + a8*x**8 + a9*x**9 

    popt, pcov = curve_fit(func,f.x, f.y)
    plt.plot(f.x, func(f.x, *popt), c='r',ls='--',lw=0.7)
    print(popt)


    def St(x):
        return 2.65677490e+00 -4.11598506e-01*x +  3.52670003e-02*x**2 -1.83464151e-03*x**3 + 6.05343182e-05*x**4 -1.28872800e-06*x**5 +1.75989973e-08*x**6 -1.48597380e-10*x**7 + 7.05235396e-13*x**8 -1.43724492e-15*x**9 
    aoac = np.linspace(0, 100, 200)
    plt.plot(aoac, St(aoac), c='g',ls='--',lw=1.7)


    fname = "naca4412_St_alpha_c." + formt
    plt.savefig(fname, format=formt, dpi=qual, bbox_inches=ajust)
    print("Saving " + fname)
    plt.close()
    print("------------------------------------------")

    fig = plt.figure()
    fig.set_size_inches(fig_width, fig_height)
    plt.ylabel(r"$\alpha_c$")
    plt.xlabel(r"$Re_c$")
    #plt.xlim(0,700); #plt.ylim(0,100)

    f = file("naca4412_2.csv")    
    plt.scatter(f.x, f.y, c='k',s=8)

    print(np.poly1d(np.polyfit(f.y,f.x, 9)))
    def func(x, a0,a1,a2,a3,a4,a5,a6,a7,a8,a9):
        #return a * np.exp(-b * x) + c
        return a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4 + a5*x**5 + a6*x**6 + a7*x**7 + a8*x**8 + a9*x**9 

    popt, pcov = curve_fit(func,f.y, f.x)
    plt.plot(func(f.y, *popt),f.y,  'r-')#,label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
    print(popt)
    def rec(x):
        return 3.83207915e+03 -6.96661369e+02*x  +5.99070559e+01*x**2 -3.00364213e+00*x**3  +9.47417747e-02*x**4 -1.93536665e-03*x**5 +2.55607403e-05*x**6 -2.10608078e-07*x**7 + 9.84008556e-10*x**8 -1.99034184e-12*x**9
    aoa = np.linspace(0, 100, 200)
    plt.plot(rec(aoa),aoa, c='y',ls='--',lw=2)

    print(St(8),rec(8))
    print(St(60),rec(60))


    fname = "naca4412_alpha_c_Re_c." + formt
    plt.savefig(fname, format=formt, dpi=qual, bbox_inches=ajust)
    print("Saving " + fname)
    plt.close()
    print("------------------------------------------")
