#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import os

params = {
    "text.usetex": False,
    "font.size": 8,
    "legend.fontsize": 8,
    "legend.handlelength": 1.0,
}
plt.rcParams.update(params)

formt = "png"
ajust = "tight"
qual = 500
fig_width = 4.3
fig_height = 3.4

# Color cycle for multiple spectra
colors = ['r', 'b', 'g', 'g', 'm', 'c', 'y']
markers = ['o', 's', '^', '<', 'p', 'v', '>']

class Spectre(object):
    def __init__(self, filename):
        print("Reading " + filename)
        data = np.transpose(np.genfromtxt(filename))
        self.R = data[0]
        self.I = data[1]
        if data.shape[0] > 2:
            self.r = data[2]
        else:
            self.r = np.zeros_like(self.R)  # Default to zeros if no residuals
        self.name = os.path.splitext(os.path.basename(filename))[0]

def plot_H(ax, R, I, r, sized=0, color="gray", symb="o", label=None):
    iflabel = False
    theta = np.linspace(0.0, 2.0 * np.pi, 400)
    ax.plot(np.cos(theta), np.sin(theta), lw=0.2, color="r", ls="-")
    for k in range(len(R)):
        if r is None or r[k] < 1.0e-6:
            mod = np.sqrt((R[k]) ** 2 + (I[k]) ** 2)
            if mod == 1:
                print("Time derivative of the baseflow found=", mod, R[k], I[k])
                plt.scatter(R[k], I[k], s=8 + sized, alpha=0.8, marker="x", facecolors=color, edgecolors=color, linewidth=0.55)
            elif mod > 1:
                plt.scatter(R[k], I[k], s=5 + sized, alpha=0.8, marker=symb, facecolors=color, edgecolors="k", linewidth=0.18)
        else:
            if iflabel == False:
                plt.scatter(R[k], I[k], s=5 + sized, alpha=1, marker=symb, facecolors=color, edgecolors="k", linewidth=0.18, label=label)
                iflabel = True
            else:
                plt.scatter(R[k], I[k], s=5 + sized, alpha=0.4, marker=symb, facecolors=color, edgecolors="k", linewidth=0.18)

    ax.axhline(y=0.0, xmin=0, xmax=1, lw=0.2, color="k", ls=":")
    ax.axvline(x=0.0, ymin=0, ymax=1, lw=0.2, color="k", ls=":")
    return

def plot_NS(ax, R, I, r, freq=False, sized=0, color="gray", symb="o", label=None):
    b = 1
    iflabel = False
    if freq:
        b = 2.0 * np.pi
    for k in range(len(R)):
        if k == 0:
            ax.axhline(y=R[k], lw=0.2, color="r", ls=":")
            ax.axvline(x=I[k] / b, lw=0.2, color="r", ls=":")
        if r is None or r[k] > 1.0e-6:
            ax.scatter(I[k] / b, R[k], s=5 + sized, alpha=0.4, marker=symb, facecolors=color, edgecolors="k", linewidth=0.18)
        else:
            if R[k] > 0:
                print("Mode: ", (k + 1))
                print(" sigma=", R[k])
                print(" omega=", I[k])
                print("     f=", round(I[k] / b, 7))
                if iflabel == False:
                    ax.scatter(I[k] / b, R[k], alpha=1, s=7 + sized, marker=symb, facecolors=color, edgecolors="k", linewidth=0.3, label=label)
                    iflabel = True
                else:
                    ax.scatter(I[k] / b, R[k], alpha=1, s=7 + sized, marker=symb, facecolors=color, edgecolors="k", linewidth=0.3)

    ax.axhline(y=0.0, xmin=0, xmax=1, lw=0.2, color="k", ls=":")
    ax.axvline(x=0.0, ymin=0, ymax=1, lw=0.2, color="k", ls=":")
    return

########################################

tolerance = 1.0e-6
if __name__ == "__main__":
    # Define folders to look for spectra
    ref_dir = ["kdim50", "kdim92", "kdim200"]

    # Marquet data
    marquet_omega = np.array([
        0.0319534, 0.0367601, 0.0770401, 0.0802253, 0.120809, 0.132061, 0.162649,
        0.181628, 0.202581, 0.241326, 0.280496, 0.318067, 0.355046, 0.394003,
        0.431520, 0.467594, 0.503667, 0.541181, 0.576779, 0.612850, 0.648923,
        0.685329, 0.720664, 0.735069, 0.755106, 0.792106, 0.825283, 0.858346, 0.892409
    ])
    marquet_sigma = np.array([
        -0.0597649, -0.0789646, -0.0654494, -0.0912914, -0.0664227, -0.100849,
        -0.0671014, -0.105598, -0.0674317, -0.0670804, -0.0666034, -0.0665784,
        -0.0660692, -0.0657232, -0.0657658, -0.0658217, -0.0655798, -0.0650630,
        -0.0660087, -0.0656648, -0.0654161, -0.0672625, -0.0699638, 0.0127451,
        -0.0727645, -0.0784841, -0.0848244, -0.0909969, -0.0974422
    ])

    # Plot H spectra
    fig = plt.figure()
    fig.set_size_inches(fig_width, fig_height)
    plt.axis("equal")
    xrz = [-1, 0, 1]
    xlabels = ["-1", "0", "1"]
    plt.xticks(xrz, xlabels)
    plt.xlim(-1.1, 1.1)
    plt.xlabel(r"$\Re (\mu)$")
    plt.yticks(xrz, xlabels)
    plt.ylim(-1.1, 1.1)
    plt.ylabel(r"$\Im (\mu)$")

    # Plot reference H spectra first (background)
    plotted = False
    for i, ref_path in enumerate(ref_dir):
        ref_h = os.path.join(ref_path, "Spectre_Hd.dat")
        if os.path.exists(ref_h):
            f_ref = Spectre(ref_h)
            plot_H(plt, f_ref.R, f_ref.I, f_ref.r, 4, colors[i], markers[i], ref_path)
            plotted = True

    # Plot current directory H spectrum last (on top)
    if os.path.exists("Spectre_Hd.dat"):
        f_current = Spectre("Spectre_Hd.dat")
        plot_H(plt, f_current.R, f_current.I, f_current.r, 12, "k", "*", "current")
        plotted = True

    if plotted:
        plt.legend(loc="best", fontsize=6)
        plt.savefig("Spectre_H_comparison." + formt, format=formt, dpi=qual, bbox_inches=ajust)
        print("Saving Spectre_H_comparison." + formt)
        plt.close()
        print("------------------------------------------")

    # Plot NS spectra
    fig = plt.figure()
    fig.set_size_inches(fig_width, fig_height)
    plt.ylabel(r"$\sigma$")
    plt.xlabel(r"$f=\omega/2\pi$")
    plt.ylim(-0.4, 0.05)

    # Plot reference NS spectra first (background)
    plotted = False
    for i, ref_path in enumerate(ref_dir):
        ref_ns = os.path.join(ref_path, "Spectre_NSd.dat")
        if os.path.exists(ref_ns):
            f_ref = Spectre(ref_ns)
            plot_NS(plt, f_ref.R, f_ref.I, f_ref.r, True, 4, colors[i], markers[i], ref_path)
            plotted = True

    # Plot current directory NS spectrum last (on top)
    if os.path.exists("Spectre_NSd.dat"):
        f_current = Spectre("Spectre_NSd.dat")
        plot_NS(plt, f_current.R, f_current.I, f_current.r, True, 12, "k", "*", "current")
        plotted = True

    if plotted:
        plt.legend(loc="best", fontsize=6)
        plt.savefig("Spectre_NS_comparison." + formt, format=formt, dpi=qual, bbox_inches=ajust)
        print("Saving Spectre_NS_comparison." + formt)
        plt.close()
        print("------------------------------------------")
