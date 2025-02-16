#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.fft import rfftfreq, rfft
params = {'text.usetex': False,
          'font.size': 11,
          'legend.fontsize': 11,
          'legend.handlelength': 2.5,
          'agg.path.chunksize':100000}
plt.rcParams.update(params)

formt = 'png'
ajust = 'tight'
qual = 400
fig_width = 3.4
fig_height = 4.3
        
def main():

    fig = plt.figure(figsize=(fig_height, fig_width))

    # Load data from file
    data = np.loadtxt('ubar.dat')
    t = data[:, 0]
    ubar = data[:, 1]

    # Plot data
    plt.xlabel('t')
    plt.ylabel('ubar')
    plt.plot(t, ubar, c='r', ls='-', lw=0.8)
    plt.axhline(y=ubar[-1], color='b', linestyle=':')
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\bar{u}$')

    diff_percent = ((ubar[-1] - ubar[0]) / ubar[0]) * 100
    print(f"Ubar percentage difference: {diff_percent}%")

    fname = 'ubar.'+formt
    print('Saving ',fname)
    plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust)

if __name__ == "__main__":
    main()