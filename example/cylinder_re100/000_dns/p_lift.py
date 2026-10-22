#!/usr/bin/env python
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os, os.path, csv, time, pickle, glob
from scipy import signal
from scipy.interpolate import interp1d

try:
    from scipy.fft import rfft, rfftfreq  # SciPy >= 1.4.0
except ImportError:
    from scipy.fftpack import rfft, rfftfreq  # SciPy < 1.4.0


def periodogram_rfft(x, fs):
    """Compute PSD using periodogram with real FFT."""
    freqs, psd = signal.periodogram(x, fs, scaling="density")
    return freqs, psd


def find_peaks(st, psd, threshold=0.01):
    """Return peak locations and values above a fraction of the maximum."""
    peak_indices = signal.find_peaks(psd, height=max(psd) * threshold)[0]
    return st[peak_indices], psd[peak_indices]


class lif_drag3d_just_dgxyz(object):
    def __init__(self, filename, flip=False):
        print("Reading " + filename)
        data = []
        with open(filename, "r") as f:
            for line in f:
                parts = line.strip().split()
                # Skip lines that are too short or not data
                if len(parts) < 6:
                    continue
                # Always extract time, dragx, dragy
                time = float(parts[1])
                dragx = float(parts[2])
                dragy = float(parts[5])
                # Try to get dragz if present (3D), else None
                dragz = float(parts[8]) if len(parts) > 8 else None
                data.append([time, dragx, dragy, dragz])
        d = np.transpose(data)
        self.t = d[0]
        self.dgx = d[1]
        self.dgy = d[2]
        self.dgz = d[3] if np.all([x is not None for x in d[3]]) else None
        if flip and self.dgz is not None:
            self.dgy, self.dgz = self.dgz, self.dgy

    #      do i = 1, nobj ! Write results
    #         if (nio == 0) then
    #            if (if3D .or. ifaxis) then ! 3D or axisymmetric
    #               write (737, "(i8,19E15.7)") istep, time,
    #  $   dragx(i), dragpx(i), dragvx(i), dragy(i), dragpy(i), dragvy(i), dragz(i), dragpz(i), dragvz(i),
    #  $   torqx(i), torqpx(i), torqvx(i), torqy(i), torqpy(i), torqvy(i), torqz(i), torqpz(i), torqvz(i)
    #            else ! 2D or not axisymmetric
    #            write (737, "(i8,10E15.7)") istep, time,
    #  $   dragx(i), dragpx(i), dragvx(i), dragy(i), dragpy(i), dragvy(i), torqz(i), torqpz(i), torqvz(i)
    #         end if ! if3D .or. ifaxis
    #      end if ! nio == 0
    #     end do ! i = 1, nobj


params = {"text.usetex": False, "font.size": 11, "legend.fontsize": 11, "legend.handlelength": 2.5, "agg.path.chunksize": 100000}
plt.rcParams.update(params)

formt = "png"
ajust = "tight"
qual = 400
fig_width = 4.3
fig_height = 16 * fig_width / 9

# Load and process lift_drag.dat using the loader class; remove all probe logic
lift_drag_file = "lift_drag.dat"
file_root = lift_drag_file.split(".")[0]
loader = lif_drag3d_just_dgxyz(lift_drag_file)
t = loader.t
# Compute means and remove them for analysis
dgx_mean = np.mean(loader.dgx)
dgx = loader.dgx - dgx_mean  # Zero-mean Cx
dgy_mean = np.mean(loader.dgy)
dgy = loader.dgy - dgy_mean  # Zero-mean Cy

# Compute RMS values for later use
dgy_rms = np.sqrt(np.mean(dgy**2))

# --- Advanced composite plot: signals, periodogram, phase space ---

# 1. Interpolate to uniform time grid for FFT/periodogram
diff_t = np.gradient(t)
dt = np.min(diff_t)
nt = int((t.max() - t.min()) / dt)
uniform_t = np.linspace(t.min(), t.max(), nt)
dgx_interp = interp1d(t, dgx, kind="linear")(uniform_t)
dgy_interp = interp1d(t, dgy, kind="linear")(uniform_t)

# 2. Compute periodogram
freq_cx, psd_cx = periodogram_rfft(dgx_interp, fs=1.0 / dt)
freq_cy, psd_cy = periodogram_rfft(dgy_interp, fs=1.0 / dt)
peak_freqs_cx, peak_vals_cx = find_peaks(freq_cx, psd_cx, threshold=0.01)
peak_freqs_cy, peak_vals_cy = find_peaks(freq_cy, psd_cy, threshold=0.01)

# 3. Compute phase space (signal vs. time derivative)
dgx_dot = np.gradient(dgx_interp, dt)
dgy_dot = np.gradient(dgy_interp, dt)

# 4. Set up gridspec for custom layout
fig = plt.figure(figsize=(fig_width * 2.2, fig_height * 1.2), constrained_layout=True)
gs = gridspec.GridSpec(2, 2, width_ratios=[2, 1], height_ratios=[1, 1.1], figure=fig)


# Reference values (Re=50)
Cx_mean_ref = 1.43
Cy_rms_ref = 0.035
St_ref = 0.125

# --- Top: signals ---
ax0 = fig.add_subplot(gs[0, :])
ax0.plot(t, dgx, label="Cx (zero-mean)", color="b")
ax0.plot(t, dgy, label="Cy (zero-mean)", color="g")
ax0.set_xlabel("t")
ax0.set_ylabel("Cx / Cy")
ax0.set_title("Time vs Cx and Cy (zero-mean)")
ax0.legend(loc="upper right")

# Add text box with comparison
stats_text = (
    f"Cx_mean = {dgx_mean:.5f} (ref: {Cx_mean_ref})\n"
    f"Cy_rms = {dgy_rms:.5f} (ref: {Cy_rms_ref})\n"
    f"St = {peak_freqs_cy[0]:.5f} (ref: {St_ref})"
    if len(peak_freqs_cy) > 0
    else f"Cx_mean = {dgx_mean:.5f} (ref: {Cx_mean_ref})\nCy_rms = {dgy_rms:.5f} (ref: {Cy_rms_ref})"
)

props = dict(boxstyle='round', facecolor='wheat', alpha=0.9)
ax0.text(0.02, 0.98, stats_text, transform=ax0.transAxes,
        verticalalignment='top', bbox=props, fontsize=9)

# --- Bottom left: periodogram ---
ax1 = fig.add_subplot(gs[1, 0])
ax1.semilogy(freq_cx, psd_cx, color="b", lw=1, label="Cx")
ax1.semilogy(freq_cy, psd_cy, color="g", lw=1, label="Cy")
ax1.scatter(peak_freqs_cx, peak_vals_cx, c="b", marker="x", s=40, label="Cx peaks")
ax1.scatter(peak_freqs_cy, peak_vals_cy, c="g", marker="x", s=40, label="Cy peaks")
if len(peak_freqs_cx) > 0:
    ax1.axvline(x=peak_freqs_cx[0], c="b", lw=0.7, ls="--", label=f"Cx $St$={peak_freqs_cx[0]:.4f}")
if len(peak_freqs_cy) > 0:
    ax1.axvline(x=peak_freqs_cy[0], c="g", lw=0.7, ls="--", label=f"Cy $St$={peak_freqs_cy[0]:.4f}")
# Reference St (Re=50)
ax1.axvline(x=0.125, color="r", linestyle=":", linewidth=1.2, label="Ref St=0.125 (Re=50)")
ax1.set_xlabel(r"$St$")
ax1.set_ylabel("PSD")
ax1.set_xlim(0.0, 1)
ax1.set_ylim(1e-3, max(psd_cx.max(), psd_cy.max()) * 1.5)
ax1.set_title("Periodogram (PSD)")
ax1.legend(loc="upper right", fontsize=9)

# --- Bottom right: phase space ---
ax2 = fig.add_subplot(gs[1, 1])
ax2.plot(dgx_interp[2:-2], dgx_dot[2:-2], color="b", lw=0.9, label="Cx")
ax2.plot(dgy_interp[2:-2], dgy_dot[2:-2], color="g", lw=0.9, label="Cy")
ax2.set_xlabel(r"Signal")
ax2.set_ylabel(r"$\dot{Signal}$")
ax2.set_title("Phase Space")
ax2.legend(fontsize=9)

figname = file_root
plt.savefig(figname + "." + formt, format=formt, dpi=qual, bbox_inches=ajust)
print("Saved: " + figname + "." + formt)
plt.close()



class lif_drag3d(object):
    def __init__(self, filename, flip=False):
        print("-------------------------------------------------------")
        print("Reading " + filename)
        # 2D: 10 columns (i, t, dgx, dpx, dvx, dgy, dpy, dvy, tqz, tpz, tvz)
        # 3D: 20 columns (i, t, dgx, dpx, dvx, dgy, dpy, dvy, dgz, dpz, dvz, tqx, tpx, tvx, tqy, tpy, tvy, tqz, tpz, tvz)
        # We'll use only drag components for plotting
        data = []
        last_line = None
        with open(filename, "r") as f:
            for line in f:
                parts = line.strip().split()
                line_wo_first = " ".join(parts[1:])
                if line_wo_first != last_line:
                    data.append([float(x) for x in parts])
                last_line = line_wo_first
        d = np.transpose(data)
        ncol = d.shape[0]
        print(f"File has {ncol} columns.")
        # 2D Nek5000 drag file: 11 columns, last three are torque (tqz, tpz, tvz)
        # 3D Nek5000 drag file: >=14 columns, includes z drag (dgz, dpz, dvz) before torque columns
        if ncol == 11:
            self.is3D = False
            print("Detected 2D file: plotting x, y drag components only (dgx, dpx, dvx, dgy, dpy, dvy)")
            print("Column mapping: [i, t, dgx, dpx, dvx, dgy, dpy, dvy, tqz, tpz, tvz]")
            drag_attrs = ["dgx", "dpx", "dvx", "dgy", "dpy", "dvy"]
            drag_indices = [2, 3, 4, 5, 6, 7]
        elif ncol >= 14:
            self.is3D = True
            print("Detected 3D file: plotting x, y, z drag components (dgx, dpx, dvx, dgy, dpy, dvy, dgz, dpz, dvz)")
            print("Column mapping: [i, t, dgx, dpx, dvx, dgy, dpy, dvy, dgz, dpz, dvz, tqx, tqy, tqz, ...]")
            drag_attrs = ["dgx", "dpx", "dvx", "dgy", "dpy", "dvy", "dgz", "dpz", "dvz"]
            drag_indices = [2, 3, 4, 5, 6, 7, 8, 9, 10]
        else:
            raise ValueError(f"Unrecognized file format: {ncol} columns. Please check the file.")
        # Assign time
        self.t = d[1]
        # Assign drag components (only those that exist)
        for attr, idx in zip(drag_attrs, drag_indices):
            if idx < ncol:
                setattr(self, attr, d[idx])
        # Optionally flip y/z if needed for 3D
        if flip and getattr(self, 'is3D', False):
            for y, z in [(5,6), (7,8)]:
                yattr, zattr = drag_attrs[y], drag_attrs[z]
                if hasattr(self, yattr) and hasattr(self, zattr):
                    arr_y, arr_z = getattr(self, yattr), getattr(self, zattr)
                    setattr(self, yattr, arr_z)
                    setattr(self, zattr, arr_y)

# --- New Figure: All d* components vs time (log y) ---
loader_full = lif_drag3d(lift_drag_file)
t = loader_full.t

# List all d* attributes to plot
# Only plot z components if the file is 3D (at least 11 columns)
# Use only the drag components present in the loader (set by loader based on 2D/3D)
if loader_full.is3D:
    components = [
        ("dgx", "Cx"), ("dpx", "Cpx"), ("dvx", "Cvx"),
        ("dgy", "Cy"), ("dpy", "Cpy"), ("dvy", "Cvy"),
        ("dgz", "Cz"), ("dpz", "Cpz"), ("dvz", "Cvz")
    ]
    colors = ["b", "c", "m", "g", "lime", "olive", "r", "orange", "k"]
else:
    components = [
        ("dgx", "Cx"), ("dpx", "Cpx"), ("dvx", "Cvx"),
        ("dgy", "Cy"), ("dpy", "Cpy"), ("dvy", "Cvy")
    ]
    colors = ["b", "c", "m", "g", "lime", "olive"]

print("Plotting components:", [label for comp, label in components if hasattr(loader_full, comp)])

fig = plt.figure()
fig.set_size_inches((fig_width * 1.1, fig_width * 0.8))
ax = fig.add_subplot(1, 1, 1)

for (comp, label), color in zip(components, colors):
    try:
        if hasattr(loader_full, comp):
            arr = getattr(loader_full, comp)
            if arr is not None and np.any(np.isfinite(arr)) and not np.all(arr == 0):
                arr = np.array(arr)
                # Do NOT remove mean
                if np.any(np.isfinite(arr)) and not np.all(np.abs(arr) < 1e-14):
                    ax.plot(t, np.abs(arr), label=f"|{label}|", color=color)
    except Exception as e:
        print(f"[WARN] Could not plot {label}: {e}")

ax.set_yscale('log')
ax.set_xlabel("t")
ax.set_ylabel("|component|")
ax.set_title("All components vs time (log y)")
ax.legend(loc="best", fontsize=8)
plt.tight_layout()
plt.savefig(file_root + "_components." + formt, format=formt, dpi=qual, bbox_inches=ajust)
print("Saved: " + file_root + "_components." + formt)
plt.close()

# # Figure
# fig, axs = plt.subplots(2, sharex=False)
# fig.set_size_inches(fig_height, fig_width)

# # Plot time series (signal on upper part of the figure)
# axs[0].set_xlabel(r'$t$')
# axs[0].set_title(f'x,y,z={coords[prb][0]},{coords[prb][1]},{coords[prb][2]}')

# # We focus on the vertical velocity signal
# nm='v' # name of the signal
# cor='r' # color of the signal
# vo = v[itmin:itmax]
# qur = np.mean(vo) # compute mean
# vo -= qur # remove mean
# print('(',nm,') mean, min, max=',round(qur,6),vo.min(),vo.max())
# vn = interp1d(to, vo, kind='nearest', fill_value="extrapolate")(tn)  # snap to nearest value

# axs[0].scatter(tn,vn,c=cor,ls='-',s=0.001) # plot signal

# # Compute PSD using periodogram_rfft
# ufreq, psd = periodogram_rfft(vn, fs=1.0/dt)

# # Find peaks in the PSD
# peak_freqs, peak_vals = find_peaks(ufreq, psd, threshold=0.01)

# # Plot PSD
# axs[1].semilogy(ufreq, psd, c=cor, lw=0.6)
# axs[1].scatter(peak_freqs, peak_vals, c='k', marker='x', s=40, label='peaks')
# if len(peak_freqs) > 0:
#     axs[1].axvline(x=peak_freqs[0], c=cor, lw=0.5, ls='--', label=f'$St={peak_freqs[0]:.4f}$')
# axs[1].set_xlabel(r'$St$')
# axs[1].set_xlim(0., 1)
# axs[1].set_ylim(1e-3, psd.max()*1.5)
# plt.legend(loc='upper right')
# fname = 'his'+str(prb+1)+'_fft.'+formt
# print('Saving ',fname); print()
# plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust)
# plt.close(); plt.clf()

# # Phase space plot
# fig = plt.figure(figsize=(fig_height, fig_width))
# plt.xlabel(r'$v$')
# plt.ylabel(r'$\dot{v}$')
# plt.plot(v[2:-2],np.gradient(v)[2:-2]/dt,c='r',ls='-',lw=0.8)

# fname = 'his'+str(prb+1)+'_phase_space.'+formt
# print('Saving ',fname)
# plt.savefig(fname,format=formt,dpi=qual,bbox_inches=ajust)
