#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import os
plt.rcParams.update({
    'text.usetex': False,
    'font.size': 8,
    'legend.fontsize': 8,
    'legend.handlelength': 2.5,
})
output_format = 'png'
figure_adjust = 'tight'
dpi = 500
fig_width = 3.5
fig_height = 2.45

class ResidualData:
    def __init__(self, filename):
        print(f'Reading {filename}')
        try:
            data = np.genfromtxt(filename).T
            self.x = data[0]
            self.y = data[1]
            if len(data) > 2:
                self.z = data[2]
            else:
                self.z = None
        except Exception as e:
            print(f"Error reading {filename}: {e}")
            self.x = self.y = self.z = None

def plot_residuals(input_file, output_file, x_label, y_label, title):
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.set_yscale('log')
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)

    f = ResidualData(input_file)
    if f.x is not None and f.y is not None and len(f.x) > 0 and len(f.y) > 0:
        try:
            if f.z is not None:
                ax.plot(f.x, f.z, c='b', lw=0.3, ls='-', marker='o', markersize=0.4)
            else:
                ax.plot(f.x, f.y, c='b', lw=0.3, ls='-', marker='o', markersize=0.4)
            fig.savefig(output_file, format=output_format, dpi=dpi, bbox_inches=figure_adjust)
            print(f'Saving {output_file}')
        except ValueError as e:
            print(f"Error plotting data: {e}")
    else:
        print(f"No valid data found in {input_file}")

    plt.close(fig)

def plot_all_residuals():
    files_to_plot = [
        ('residu.dat', 'residu.png', 'Linearized calls', 'Residual', 'Newton Residuals'),
        ('residu_newton.dat', 'residu_newton.png', 'Iteration', 'Residual', 'Newton Residuals'),
        ('residu_gmres.dat', 'residu_gmres.png', 'Iteration', 'Residual', 'GMRES Residuals'),
        ('residu_arnoldi.dat', 'residu_arnoldi.png', 'Iteration', 'Residual', 'Arnoldi Residuals')
    ]

    for input_file, output_file, x_label, y_label, title in files_to_plot:
        if os.path.exists(input_file):
            plot_residuals(input_file, output_file, x_label, y_label, title)
        else:
            print(f"File {input_file} not found.")

if __name__ == '__main__':
    plot_all_residuals()
