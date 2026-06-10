#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.gridspec import GridSpec

plt.rcParams.update({
    'text.usetex': False,
    'font.size': 8,
    'legend.fontsize': 8,
    'legend.handlelength': 2.5,
})
output_format = 'png'
figure_adjust = 'tight'
dpi = 1000  # For retina quality

fig_width = 4 # inches
golden_ratio = 1.618
fig_height = fig_width / golden_ratio  # ≈2.16 inches

class ResidualData:
    def __init__(self, filename):
        """Initialize ResidualData from file."""
        self.iter = None
        self.total_calls = None
        self.iter_calls = None
        self.k_sum = None
        self.time = None
        self.tol = None
        self.residual = None
        self.dtol = None
        self.initial_residual = None

        if os.path.exists(filename):
            try:
                with open(filename, 'r') as f:
                    # Read initial residual from header
                    header = f.readline()
                    if header.startswith('# Initial residual:'):
                        self.initial_residual = float(header.split(':')[1])
                        # Skip column header line
                        f.readline()
                    else:
                        # Old format, rewind file
                        f.seek(0)
                    
                    # Read the rest of the data
                    data = np.loadtxt(f)
                    if data.size > 0:  # Check if we have any data
                        if len(data.shape) == 1:  # Single row
                            data = data.reshape(1, -1)
                        if filename == 'residu_gmres.dat':
                            # For GMRES: [newton_iter, gmres_iter, k, k_sum, tol, beta2, dtol]
                            if len(data[0]) >= 7:  # Check if we have all required columns
                                self.newton_iter = data[:, 0]
                                self.gmres_iter = data[:, 1]
                                self.k = data[:, 2]
                                self.k_sum = data[:, 3]  # k_sum is column 4
                                self.total_calls = data[:, 3]  # Use k_sum as total_calls for GMRES
                                self.tol = data[:, 4]
                                self.residual = data[:, 5]
                                self.dtol = data[:, 6]
                            else:
                                print(f"File {filename} has insufficient columns")
                                self.residual = None
                                self.tol = None
                        elif filename == 'residu_arnoldi.dat':
                            # For Arnoldi: [k, arnoldi_iter_time, tol, beta2, dtol]
                            if len(data[0]) >= 5:  # Check if we have all required columns
                                self.k = data[:, 0]
                                self.time = data[:, 1]
                                self.tol = data[:, 2]
                                self.residual = data[:, 3]
                                self.dtol = data[:, 4]
                            else:
                                print(f"File {filename} has insufficient columns")
                                self.residual = None
                                self.tol = None
                        elif filename == 'residu_newton.dat':
                            # For Newton: [i, total_calls, iter_calls, k_sum, tottime, solver_tol, residual, dtol]
                            if len(data[0]) >= 8:
                                self.iter = data[:, 0]
                                self.total_calls = data[:, 1]  # Use total_calls directly from Newton file
                                self.iter_calls = data[:, 2]
                                self.k_sum = data[:, 3]  # k_sum is already accumulated in the data file
                                self.time = data[:, 4]
                                self.tol = data[:, 5]
                                self.residual = data[:, 6]
                                self.dtol = data[:, 7]
                            else:
                                print(f"File {filename} has insufficient columns")
                                self.residual = None
                                self.tol = None
                        else:
                            self.x = data[:, 0]
                            self.y = data[:, 1]
                            if len(data[0]) > 2:
                                self.z = data[:, 2]
                            else:
                                self.z = None
            except Exception as e:
                print(f"Error reading {filename}: {e}")
                self.residual = None
                self.tol = None

def plot_arnoldi_residuals(data, output_file):
    if data.residual is None:
        print("No valid Arnoldi data to plot")
        return
        
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.set_yscale('log')
    
    # Lower x-axis (Arnoldi iterations)
    ax.set_xlabel('Arnoldi Iterations')
    ax.set_ylabel(r'$\Vert r \Vert^2$')
    
    # Colors
    gmres_color = 'purple'
    newton_color = 'blue'
    
    # Plot Arnoldi residuals
    x = np.arange(len(data.residual)) if isinstance(data.residual, np.ndarray) else [0]
    ax.plot(x, data.residual if isinstance(data.residual, np.ndarray) else [data.residual], 
            c='k', lw=0.6, ls='-', marker='.', markersize=2, label='Arnoldi', zorder=3)
    
    # Plot tolerance
    if data.tol is not None:
        ax.plot(x, data.tol if isinstance(data.tol, np.ndarray) else [data.tol], 'b--', lw=0.8, label='tol', zorder=2)
    
    # Plot Newton residuals if available
    if os.path.exists('residu_newton.dat'):
        newton_data = ResidualData('residu_newton.dat')
        if newton_data.residual is not None:
            # Plot Newton residuals
            ax.plot(newton_data.k_sum, newton_data.residual,c=newton_color, lw=0.8, ls='none', marker='s', markersize=4, markerfacecolor='none', label='Newton', zorder=4)
            
            # Plot dtol
            if newton_data.dtol is not None:
                ax.axhline(y=newton_data.dtol if isinstance(newton_data.dtol, float) else newton_data.dtol[0], color='red', linestyle='--', linewidth=0.6, label='dtol', zorder=1)
    
    # Add grid
    ax.grid(True, which='both', linestyle=':', linewidth=0.5, alpha=0.5)
    ax.legend(ncol=2, fontsize=6)
    
    print(f'Saving {output_file}')
    fig.savefig(output_file, format=output_format, dpi=dpi, bbox_inches=figure_adjust)
    plt.close(fig)

def plot_newton_krylov_stages():
    # Create figure with two subplots
    fig = plt.figure(figsize=(fig_width, 2*fig_height))
    gs = GridSpec(2, 1, height_ratios=[1, 1], hspace=0.2)
    
    # Top subplot for Newton residuals
    ax1 = fig.add_subplot(gs[0])
    ax1.set_yscale('log')
    ax1.set_ylabel(r'$\Vert r \Vert^2$', labelpad=2)
    ax1.grid(True, which='both', ls='-', alpha=0.2)
    
    # Bottom subplot for Arnoldi/GMRES residuals
    ax2 = fig.add_subplot(gs[1])
    ax2.set_yscale('log')
    ax2.set_xlabel('Matrix-vector products')
    ax2.set_ylabel(r'$\Vert r \Vert^2$')
    ax2.grid(True, which='both', ls='-', alpha=0.2)
    
    # Read data files
    newton_data = ResidualData('residu_newton.dat')

    if newton_data.residual is not None:
        # Plot initial residual if available
        if newton_data.initial_residual is not None:
            ax1.scatter(0, newton_data.initial_residual, color='k', marker='o', 
                      label='Initial', zorder=3)
        
        # Plot Newton data
        x = newton_data.k_sum
        ax1.plot(x, newton_data.residual, 'bo-', label='Newton', markersize=2)
        ax1.plot(x, newton_data.dtol, 'r--', label='Target (dtol)', alpha=0.7)
        ax1.plot(x, newton_data.tol, 'g--', label='GMRES (tol)', alpha=0.7)
        
        # Add Newton iteration markers
        for i, (k_sum, res) in enumerate(zip(newton_data.k_sum, newton_data.residual)):
            ax1.annotate(f'N{i}', (k_sum, res), xytext=(2, 2), 
                        textcoords='offset points', fontsize=6)
            
        # Add vertical lines for Newton iterations
        for k_sum in newton_data.k_sum:
            ax1.axvline(x=k_sum, color='b', linestyle=':', alpha=0.2)
    
    # Plot Arnoldi and GMRES data in bottom subplot
    arnoldi_data = ResidualData('residu_arnoldi.dat')
    gmres_data = ResidualData('residu_gmres.dat')
    
    if arnoldi_data.residual is not None:
        # Use k_sum for x-axis to match Newton plot
        x_arnoldi = np.cumsum(arnoldi_data.k)  # Accumulate k values
        ax2.plot(x_arnoldi, arnoldi_data.residual, 'ro-', 
                label='Arnoldi', markersize=1)
        ax2.plot(x_arnoldi, arnoldi_data.tol, 'r--', 
                label='Arnoldi tol', alpha=0.5)
        ax2.axhline(y=newton_data.dtol[0], color='r', linestyle='--', label='Target (dtol)', alpha=0.7)
        ax2.set_xlabel('Matrix-vector products')  # Update label to match

    # Adjust plot aesthetics
    ax1.legend(loc='best', ncol=2)
    ax2.legend(loc='best', ncol=2)
    
    print('Saving residu2.' + output_format)
    plt.savefig('residu2.' + output_format, dpi=dpi, bbox_inches=figure_adjust)
    plt.close()

if __name__ == '__main__':
    # Plot GMRES residuals if file exists
    if os.path.exists('residu_gmres.dat'):
        gmres_data = ResidualData('residu_gmres.dat')

    # Plot Arnoldi residuals if file exists
    if os.path.exists('residu_arnoldi.dat'):
        arnoldi_data = ResidualData('residu_arnoldi.dat')
        if arnoldi_data.residual is not None:
            plot_arnoldi_residuals(arnoldi_data, 'residu.' + output_format)
        
    # Plot Newton-Krylov stages
    plot_newton_krylov_stages()