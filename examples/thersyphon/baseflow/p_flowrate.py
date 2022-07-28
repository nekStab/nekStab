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
import os, sys, traceback
import os.path
import csv
import pickle

params = {
    "text.usetex": False,
    "legend.handlelength": 1.0,
}
plt.rcParams.update(params)
new_prop_cycle = cycler("color", ["k", "r"]) * cycler("linewidth", [1.0, 1.5, 2.0])

plt.rc("axes", prop_cycle=new_prop_cycle)
plt.style.use("seaborn-white")

#':' :,'-.' dashdot, '--' --, '-' -
# b: blue, g: green, r: red,c: cyan,m: magenta,y: yellow,k: black

my_color = "gray"  # ['red', 'black', 'blue', 'brown', 'green']
formt = "png"
ajust = "tight"
qual = 300
fig_width = 8
fig_height = 6
mksize = 15  # marker size
ftsize = 6  # maker number font size
lgsize = 4  # legend size
syset = -3  # offset y label

# https://stackoverflow.com/questions/16006572/plotting-different-colors-in-matplotlib/45916274
# default color cycle is ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"], the Vega category10 palette.


def extrap(x, fp, xp):
    return (xp[-2] - xp[-1]) / (fp[-2] - fp[-1]) * (x - fp[-1]) + xp[-1]


class SpectreH(object):
    def __init__(self, filename):
        # data = np.transpose(np.loadtxt(filename))
        print("Reading " + filename)
        data = np.transpose(np.genfromtxt(filename))
        self.r = data[0]
        self.i = data[1]
        self.residu = data[2]
        del data


class SNS(object):
    def __init__(self, filename):
        print("Reading " + filename)
        # data = np.transpose(np.loadtxt(filename))
        data = np.transpose(np.genfromtxt(filename))
        self.r = data[0]
        self.i = data[1]
        del data

class FR(object):
    def __init__(self, filename):
        data = np.transpose(np.genfromtxt(filename))
        self.t = data[0]
        self.fr = data[1]
        del data
        #print("Read " + filename)


class SNS2(object):
    def __init__(self, filename):
        print("Reading " + filename)
        # data = np.loadtxt(filename)
        data = np.transpose(np.genfromtxt(filename))
        self.r = data[0]
        self.i = data[1]
        del data


def plot_H(ax, r, i, residu, sized=0, my_color="gray", label=None):
    iflabel = False
    theta = np.linspace(0.0, 2.0 * np.pi, 400)
    ax.plot(np.cos(theta), np.sin(theta), lw=0.2, color="r", ls="-")
    for k in range(len(r)):
        if residu[k] < tolerance:
            mod = np.sqrt((r[k]) ** 2 + (i[k]) ** 2)
            if mod == 1:
                print("Time derivative of the baseflow found=", mod, r[k], i[k])
                plt.scatter(
                    r[k],
                    i[k],
                    s=mksize + 3 + sized,
                    alpha=0.8,
                    marker="x",
                    facecolors=my_color,
                    edgecolors=my_color,
                    linewidth=0.55,
                )
            elif mod > 1:
                print("Mode: ", (k + 1), "|", mod, "|")
                print("   Re=", r[k])
                print("angle: ", np.angle(i[k]))
                print("    f=", (i[k] / (2.0 * np.pi)))
                print(
                    "--------------------------------------------------------------------"
                )
                if iflabel == False:
                    plt.scatter(
                        r[k],
                        i[k],
                        s=mksize + sized,
                        alpha=0.8,
                        marker="o",
                        facecolors="none",
                        edgecolors="k",
                        linewidth=0.18,
                        label=label,
                    )
                    iflabel = True
                else:
                    plt.scatter(
                        r[k],
                        i[k],
                        s=mksize + sized,
                        alpha=0.8,
                        marker="o",
                        facecolors="none",
                        edgecolors="k",
                        linewidth=0.18,
                    )
            else:
                if iflabel == False:
                    plt.scatter(
                        r[k],
                        i[k],
                        s=mksize + sized,
                        alpha=0.8,
                        marker="o",
                        facecolors=color,
                        edgecolors="k",
                        linewidth=0.18,
                        label=label,
                    )
                    iflabel = True
                else:
                    plt.scatter(
                        r[k],
                        i[k],
                        s=mksize + sized,
                        alpha=0.8,
                        marker="o",
                        facecolors=color,
                        edgecolors="k",
                        linewidth=0.18,
                    )
        else:
            plt.scatter(r[k], i[k], s=1, alpha=0.5, color=color)
    ax.axhline(y=0.0, xmin=0, xmax=1, lw=0.2, color="k", ls=":")
    ax.axvline(x=0.0, ymin=0, ymax=1, lw=0.2, color="k", ls=":")
    return


def pNS(ax, r, i, bds, freq=False, sized=0, my_color="gray", label=None):
    b = 1
    iflabel = False
    if freq:
        b = 2.0 * np.pi
    for k in range(len(r)):
        tr = 0
        ta = 1e-10
        if r[k] > 0:
            print("Mode: ", (k + 1), r[k], i[k])
            print("      f=", i[k] / b)
            if isclose(abs(r[k]), 0.0, rel_tol=tr, abs_tol=ta) and isclose(
                abs(i[k]), 0.0, rel_tol=tr, abs_tol=ta
            ):
                print("skipping plot of spurious mode 1")
                pass
            else:
                ax.scatter(
                    i[k] / b,
                    r[k],
                    alpha=0.8,
                    s=mksize + 3 + sized,
                    marker="o",
                    facecolors=my_color,
                    edgecolors="k",
                    linewidth=0.3,
                )
        else:
            if iflabel == False:
                if isclose(abs(r[k]), 0.0, rel_tol=tr, abs_tol=ta) and isclose(
                    abs(i[k]), 0.0, rel_tol=tr, abs_tol=ta
                ):
                    print("skipping plot of spurious mode 2")
                    iflabel = False
                    pass
                else:
                    ax.scatter(
                        i[k] / b,
                        r[k],
                        alpha=0.8,
                        s=mksize + sized,
                        marker="o",
                        facecolors=my_color,
                        edgecolors="k",
                        linewidth=0.2,
                        label=label,
                    )
                    iflabel = True
            else:
                if isclose(abs(r[k]), 0.0, rel_tol=tr, abs_tol=ta) and isclose(
                    abs(i[k]), 0.0, rel_tol=tr, abs_tol=ta
                ):
                    print("skipping plot of spurious mode 3")
                    pass
                else:
                    ax.scatter(
                        i[k] / b,
                        r[k],
                        alpha=0.8,
                        s=mksize + sized,
                        marker="o",
                        facecolors=my_color,
                        edgecolors="k",
                        linewidth=0.2,
                    )
    ax.axhline(y=0.0, xmin=0, xmax=1, lw=0.2, color="k", ls=":")
    ax.axvline(x=0.0, ymin=0, ymax=1, lw=0.2, color="k", ls=":")
    return


def pNSn(ax, r, i, bds, freq=False, sized=0, my_color="gray", label=None):
    b = 1
    iflabel = False
    if freq:
        b = 2.0 * np.pi
    for k in range(len(r)):
        if abs(r[k]) < 1e-7 and abs(i[k]) < 1e-7:
            print("skipping spurious mode annotate", k + 1)
            # ax.scatter(file.i[k]/b,file.r[k], alpha=0.8, s=5+sized,marker='x', facecolors=color,edgecolors='k',linewidth=0.3)
        else:
            if r[k] > 0:
                ax.scatter(
                    i[k] / b,
                    r[k],
                    alpha=0.8,
                    s=mksize + 3 + sized,
                    marker="o",
                    facecolors=my_color,
                    edgecolors="k",
                    linewidth=0.3,
                )
                ax.annotate(
                    str(k + 1),
                    (i[k] / b, r[k]),
                    color=my_color,
                    fontsize=ftsize + float(0.5 * sized),
                )

            else:
                if iflabel == False:
                    ax.scatter(
                        i[k] / b,
                        r[k],
                        alpha=0.8,
                        s=mksize + sized,
                        marker="o",
                        facecolors=my_color,
                        edgecolors="k",
                        linewidth=0.2,
                        label=label,
                    )
                    ax.annotate(
                        str(k + 1),
                        (i[k] / b, r[k]),
                        color=my_color,
                        fontsize=ftsize + float(0.5 * sized),
                    )
                    iflabel = True
                else:
                    if bds[0] <= i[k] / b <= bds[1] and bds[2] <= r[k] <= bds[3]:
                        ax.scatter(
                            i[k] / b,
                            r[k],
                            alpha=0.8,
                            s=mksize + sized,
                            marker="o",
                            facecolors=my_color,
                            edgecolors="k",
                            linewidth=0.2,
                        )
                        ax.annotate(
                            str(k + 1),
                            (i[k] / b, r[k]),
                            color=my_color,
                            fontsize=ftsize + float(0.5 * sized),
                        )

    ax.axvline(x=0.0, lw=0.1, color="k", ls="--")
    return


if __name__ == "__main__":

    fig = plt.figure()
    fig.set_size_inches(fig_width, fig_height)
    plt.ylabel(r"Flow rate", labelpad=syset)
    plt.xlabel(r"$t$")    


    data=FR("flow_rate.dat")
    plt.plot(data.t, data.fr,ls='-',lw=1,c='k',label='Ra=20k')
    
    fname = "flow_rate." + formt    
    plt.legend(loc="center right", fontsize=8)
    plt.savefig(fname, format=formt, dpi=qual, bbox_inches=ajust)
    print("Saving " + fname)
    plt.close()
    print("------------------------------------------")
