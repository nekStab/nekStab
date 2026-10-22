#!/usr/bin/env python3
"""plot.py: compare the DNS of the unstable Re=1950 limit cycle with its Floquet pair.

The DNS starts on the orbit plus a small part of the leading Floquet mode
(make_ic.py). Two checks use the probe file cav.his:

(a) Stroboscopic sampling. The probe is read once per orbit period T. For a
    Floquet multiplier mu the deviation from the orbit follows mu^n, so
    its growth rate is ln|mu|/T and its frequency is arg(mu)/T. The script fits
    d_n = a + b n + Re(c mu^n) and compares sigma and omega with the Floquet
    pair of ../311_stability_direct_floquet. The term b n is the neutral
    multiplier mu = 1 (shift along the orbit) that a period known to seven
    digits leaves in the samples. The window starts after the stable modes
    have decayed. It ends at the last period for which the relative fit
    residual stays below RESID_LIMIT: past that point the deviation is no
    longer small and the single-multiplier model fails. Panel (c) shows sigma
    for every window end, so the choice of the end can be judged.
(b) Spectrum. A secondary mode with the Floquet frequency f_F = omega/(2 pi)
    puts sidebands at n f_1 +- f_F around the harmonics of the orbit frequency
    f_1 = 1/T.

OUTPUTS: plot_floquet_vs_dns.png, floquet_vs_dns.dat and a printed table.
USAGE:   python plot.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

CASE_DIR = Path(__file__).resolve().parent
FLOQUET_DIR = CASE_DIR.parent / "311_stability_direct_floquet"
PERIOD = 10.78909  # orbit period, written in the time stamp of the base flow
N_FIT_START = 20  # periods dropped at the start: stable modes still decay
N_FIT_MIN = 20  # shortest window, in periods
RESID_LIMIT = 0.10  # largest relative residual of the single-multiplier fit


def read_probes(path):
    """Return t[n] and the array x[n, point, component] for u, v, w."""
    rows = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) == 5:
            rows.append([float(x.replace("D", "E")) for x in fields])
    data = np.array(rows)
    npts = int(path.read_text().splitlines()[0])
    n = (len(data) // npts) * npts
    data = data[:n].reshape(-1, npts, 5)
    return data[:, 0, 0], data[:, :, 1:4]


def floquet_pair():
    """Leading Floquet pair (sigma, omega) from the direct Floquet stage."""
    for line in (FLOQUET_DIR / "Spectre_NSd_conv.dat").read_text().splitlines():
        fields = line.split()
        if len(fields) >= 2:
            return float(fields[0]), abs(float(fields[1]))
    raise SystemExit("no Floquet pair found")


def fit_multiplier(d):
    """Fit d_n = a + b n + Re(c mu^n) over mu = exp(s + i phi).

    A coarse grid over (s, phi) is refined once around its best point; for each
    (s, phi) the other coefficients follow from linear least squares.
    Returns s = sigma T, phi = omega T and the relative residual.
    """
    n = np.arange(len(d), dtype=float)

    def residual(s, phi):
        z = np.exp((s + 1j * phi) * n)
        a = np.column_stack([np.ones_like(n), n, z.real, z.imag])
        coef, *_ = np.linalg.lstsq(a, d, rcond=None)
        return np.linalg.norm(a @ coef - d)

    def search(s_grid, phi_grid):
        best = (np.inf, 0.0, 0.0)
        for s in s_grid:
            for phi in phi_grid:
                r = residual(s, phi)
                if r < best[0]:
                    best = (r, s, phi)
        return best

    r, s, phi = search(np.linspace(-0.2, 0.4, 61), np.linspace(0.0, np.pi, 315))
    ds, dphi = 0.6 / 60, np.pi / 314
    r, s, phi = search(
        np.linspace(s - ds, s + ds, 21), np.linspace(phi - dphi, phi + dphi, 21)
    )
    return s, phi, r / np.linalg.norm(d - d.mean())


def main():
    plt.switch_backend("Agg")
    t, x = read_probes(CASE_DIR / "cav.his")
    sigma_f, omega_f = floquet_pair()
    f1, f_f = 1.0 / PERIOD, omega_f / (2.0 * np.pi)
    nper = int(t[-1] // PERIOD)
    tn = PERIOD * np.arange(nper + 1)

    # stroboscopic deviation of every probe component; keep the one with the largest growth
    best = None
    for ip in range(x.shape[1]):
        for ic in range(3):
            xn = np.interp(tn, t, x[:, ip, ic])
            d = xn - xn[0]
            if best is None or np.abs(d).max() > np.abs(best[2]).max():
                best = (ip, ic, d)
    ip, ic, d = best
    dev = np.abs(d - np.median(d[N_FIT_START : N_FIT_START + 20]))
    # sigma, omega and residual for every window end; keep the last end below RESID_LIMIT
    ends = np.arange(N_FIT_START + N_FIT_MIN, len(d) + 1)
    sweep = [(e, *fit_multiplier(d[N_FIT_START:e])) for e in ends]
    end = ends[0]
    for e, _, _, resid in sweep:
        if resid > RESID_LIMIT:
            break
        end = e
    window = slice(N_FIT_START, end)
    s, phi, resid = next(item[1:] for item in sweep if item[0] == end)
    sigma_d, omega_d = s / PERIOD, phi / PERIOD

    print(
        f"probe {ip + 1}, component {'uvw'[ic]}; window n = {N_FIT_START}..{end - 1} "
        f"(t = {tn[N_FIT_START]:.0f}..{tn[end - 1]:.0f}); relative fit residual {resid:.2f}"
    )
    print("          sigma [1/t.u.]   omega [rad/t.u.]   f [1/t.u.]")
    print(f"Floquet   {sigma_f:+.4e}     {omega_f:.5f}          {f_f:.5f}")
    print(
        f"DNS       {sigma_d:+.4e}     {omega_d:.5f}          {omega_d / (2 * np.pi):.5f}"
    )

    # numbers for ref/reference.json: Floquet and DNS growth rate and frequency
    np.savetxt(
        CASE_DIR / "floquet_vs_dns.dat",
        [[sigma_f, sigma_d, omega_f, omega_d, abs(sigma_d / sigma_f - 1.0), abs(omega_d / omega_f - 1.0)]],
        fmt="%.6e",
        header="sigma_floquet sigma_dns omega_floquet omega_dns relerr_sigma relerr_omega",
    )

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 4))
    names = "uvw"
    ax1.semilogy(
        tn,
        dev + 1e-30,
        "-o",
        ms=3,
        mfc="none",
        lw=0.8,
        label=f"DNS, probe {ip + 1}, {names[ic]}",
    )
    tw = tn[window]
    ref = dev[N_FIT_START]
    ax1.semilogy(
        tw,
        ref * np.exp(sigma_f * (tw - tw[0])),
        "--s",
        ms=3,
        mfc="none",
        lw=1.6,
        markevery=4,
        label=f"Floquet, sigma = {sigma_f:.2e}",
    )
    ax1.semilogy(
        tw,
        ref * np.exp(sigma_d * (tw - tw[0])),
        ":^",
        ms=3,
        mfc="none",
        lw=1.0,
        markevery=4,
        label=f"DNS fit, sigma = {sigma_d:.2e}",
    )
    ax1.axvline(tn[end - 1], color="gray", lw=0.6, label="end of linear window")
    ax1.set_xlabel("time [t.u.]")
    ax1.set_ylabel("deviation from the orbit, once per period")
    ax1.set_title("(a) stroboscopic growth", loc="left", fontsize=10)
    ax1.legend(fontsize=7)

    # spectrum of the same signal over the whole run, Hann window
    sel = t > 200.0
    tt, sig = t[sel], x[sel, ip, ic]
    dt = np.median(np.diff(tt))
    tu = np.arange(tt[0], tt[-1], dt)
    sig = np.interp(tu, tt, sig)
    sig = sig - sig.mean()
    nfft = 1 << int(np.ceil(np.log2(len(sig) * 2)))
    psd = np.abs(np.fft.rfft(sig * np.hanning(len(sig)), n=nfft)) ** 2
    fr = np.fft.rfftfreq(nfft, d=dt)
    m = fr < 0.45
    ax2.semilogy(
        fr[m], psd[m], "-o", ms=2, mfc="none", lw=0.6, markevery=40, label="DNS probe"
    )
    first = True
    for k in range(0, 5):
        for sgn in (-1, 1):
            f = k * f1 + sgn * f_f
            if 0.0 < f < 0.45:
                ax2.axvline(
                    f,
                    color="tab:red",
                    lw=0.6,
                    ls="--",
                    label=r"$n f_1 \pm f_F$ (Floquet)" if first else None,
                )
                first = False
    ax2.axvline(f1, color="k", lw=0.6, label=f"$f_1$ = {f1:.4f}")
    ax2.set_xlabel("frequency [1/t.u.]")
    ax2.set_ylabel("power")
    ax2.set_title("(b) spectrum", loc="left", fontsize=10)
    ax2.legend(fontsize=7)

    # sigma of the fit for every window end, with the Floquet value as a guide
    t_end = np.array([tn[e - 1] for e, *_ in sweep])
    sig_end = np.array([r[1] / PERIOD for r in sweep])
    good = np.array([r[3] <= RESID_LIMIT for r in sweep])
    ax3.plot(
        t_end[good],
        sig_end[good],
        "-o",
        ms=3,
        lw=1.0,
        markevery=3,
        label=f"DNS fit, residual <= {RESID_LIMIT:.2f}",
    )
    ax3.plot(
        t_end[~good],
        sig_end[~good],
        "o",
        ms=3,
        mfc="none",
        lw=0,
        label=f"DNS fit, residual > {RESID_LIMIT:.2f}",
    )
    ax3.axhline(sigma_f, color="tab:orange", lw=1.6, label="Floquet")
    ax3.axvline(tn[end - 1], color="gray", lw=0.6, label="window end used in (a)")
    ax3.set_ylim(-1e-3, 8e-3)
    ax3.set_xlabel("end of the fit window [t.u.]")
    ax3.set_ylabel("sigma [1/t.u.]")
    ax3.set_title("(c) sigma versus window end", loc="left", fontsize=10)
    ax3.legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(CASE_DIR / "plot_floquet_vs_dns.png", dpi=400)
    print("Saved plot_floquet_vs_dns.png")


if __name__ == "__main__":
    main()
