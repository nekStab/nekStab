#!/usr/bin/env python3
"""
validate_eigensolver.py: Eigensolver validation for cylinder at Re=50.

Runs multiple Arnoldi/Krylov-Schur cases varying k_dim, collects spectra
and cost metrics, generates a publication-quality comparison figure, and
validates converged eigenvalues against Marquet et al. (2009).

OUTPUTS: validate_eigensolver.png
USAGE:
    python validate_eigensolver.py              # Compile + run all + plot
    python validate_eigensolver.py --plot-only  # Plot from existing results/
    python validate_eigensolver.py --dry-run    # Show commands without executing
    python validate_eigensolver.py -n 8         # Override MPI process count
    python validate_eigensolver.py --cases k40 k200  # Run subset
"""
from pathlib import Path
import argparse
import os
import re
import shutil
import subprocess
import sys
import time

import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SCRIPT = Path(__file__).resolve()
CASE_DIR = SCRIPT.parent
CASE_NAME = "1cyl"
PAR_FILE = CASE_DIR / f"{CASE_NAME}.par"
RESULTS_DIR = CASE_DIR / "results"
OUTPUT_PNG = CASE_DIR / "validate_eigensolver.png"

# Each case: (label, k_dim, schur_tgt, maxmodes)
# schur_tgt=0 and maxmodes=0 use the .usr defaults (2 and schur_tgt).
CASES = [
    ("k20",       20,   0,  0),
    ("k30",       30,   0,  0),
    ("k40",       40,   0,  0),
    ("k50_t10",   50,  10, 10),
    ("k50_t20",   50,  20, 20),
    ("k60",       60,   0,  0),
    ("k100",     100,   0,  0),
    ("k200",     200,   0,  0),
]

# Marquet et al. (2009) reference eigenvalues (NS-plane: sigma, omega)
MARQUET_SIGMA = np.array([
    -0.0597649, -0.0789646, -0.0654494, -0.0912914, -0.0664227, -0.100849,
    -0.0671014, -0.105598, -0.0674317, -0.0670804, -0.0666034, -0.0665784,
    -0.0660692, -0.0657232, -0.0657658, -0.0658217, -0.0655798, -0.0650630,
    -0.0660087, -0.0656648, -0.0654161, -0.0672625, -0.0699638,  0.0127451,
    -0.0727645, -0.0784841, -0.0848244, -0.0909969, -0.0974422,
])
MARQUET_OMEGA = np.array([
    0.0319534, 0.0367601, 0.0770401, 0.0802253, 0.120809, 0.132061, 0.162649,
    0.181628, 0.202581, 0.241326, 0.280496, 0.318067, 0.355046, 0.394003,
    0.431520, 0.467594, 0.503667, 0.541181, 0.576779, 0.612850, 0.648923,
    0.685329, 0.720664, 0.735069, 0.755106, 0.792106, 0.825283, 0.858346,
    0.892409,
])

EIGENTOL = 1.0e-6  # convergence threshold for residuals

# Plot styling — colorblind-safe
STYLE = {
    "k20":      {"color": "tab:red",    "marker": "o"},
    "k30":      {"color": "tab:orange", "marker": "s"},
    "k40":      {"color": "tab:blue",   "marker": "^"},
    "k50":      {"color": "tab:red",    "marker": "o"},
    "k50_t10":  {"color": "tab:olive",  "marker": "P"},
    "k50_t20":  {"color": "tab:brown",  "marker": "X"},
    "k60":      {"color": "tab:green",  "marker": "D"},
    "k92":      {"color": "tab:cyan",   "marker": "s"},
    "k100":     {"color": "tab:purple", "marker": "v"},
    "k200":     {"color": "black",      "marker": "*"},
}

# ---------------------------------------------------------------------------
# Platform helpers
# ---------------------------------------------------------------------------

def detect_nekstab_root():
    """Walk upward from CASE_DIR to find the nekStab root."""
    env = os.environ.get("NEKSTAB_SOURCE_ROOT")
    if env and Path(env).is_dir():
        return Path(env)
    # Walk up: CASE_DIR = <root>/example/cylinder/stability/direct
    candidate = CASE_DIR.parents[3]
    if (candidate / "src" / "NEKSTAB").exists():
        return candidate
    sys.exit("Cannot find nekStab root. Set NEKSTAB_SOURCE_ROOT.")


def detect_nprocs(nekstab_root):
    """Detect optimal MPI process count from mesh and CPU topology."""
    re2 = CASE_DIR / f"{CASE_NAME}.re2"
    nelgt = None
    if re2.exists():
        with open(re2, "rb") as f:
            header = f.read(80)
            # re2 header: "#v002  NELGT  NDIM  NELGT ..."
            # Skip version prefix (#vNNN), then parse NELGT
            hdr_str = header.decode("ascii", errors="replace")
            tokens = hdr_str.split()
            if len(tokens) >= 2:
                try:
                    nelgt = int(tokens[1])
                except ValueError:
                    pass

    # Detect P-cores on hybrid Intel (e.g., i7-14700: 8P + 12E)
    try:
        out = subprocess.check_output(
            ["lscpu"], text=True, stderr=subprocess.DEVNULL
        )
        # Check for Intel hybrid: if "Model name" has i7/i9 with 12th+ gen
        # use heuristic: P-cores = total_cores - E-cores
        # Simpler: look for "Core(s) per socket" and cap at 8 for safety
        m = re.search(r"Core\(s\) per socket:\s+(\d+)", out)
        physical_cores = int(m.group(1)) if m else (os.cpu_count() or 4)
    except (FileNotFoundError, subprocess.CalledProcessError):
        physical_cores = os.cpu_count() or 4

    # Heuristic: >=200 elements/rank, cap at 8 (P-cores on hybrid Intel)
    p_cores = min(physical_cores, 8)
    if nelgt:
        max_from_mesh = max(1, nelgt // 200)
        nprocs = min(p_cores, max_from_mesh)
    else:
        nprocs = p_cores

    return nprocs


# ---------------------------------------------------------------------------
# Build / Run
# ---------------------------------------------------------------------------

def compile_once(nekstab_root, dry_run=False):
    """Compile the case once (k_dim is runtime via .par)."""
    mks = nekstab_root / "bin" / "mks"
    if not mks.exists():
        sys.exit(f"mks not found at {mks}")

    env = os.environ.copy()
    env["NEKSTAB_SOURCE_ROOT"] = str(nekstab_root)
    env["NEK_SOURCE_ROOT"] = str(nekstab_root / "Nek5000")

    cmd = [str(mks), CASE_NAME]
    if dry_run:
        print(f"[DRY-RUN] cd {CASE_DIR} && {' '.join(cmd)}")
        return True

    # If binary already exists, skip recompilation (k_dim changes at runtime).
    binary = CASE_DIR / "nek5000"
    if binary.exists():
        print(f"Binary exists ({binary}), skipping compilation.")
        return True

    print(f"=== Compiling {CASE_NAME} ===")
    # Pipe "N\n" to stdin: Nek5000's makenek.inc may prompt interactively
    # ("rebuild 3rd party deps?") when .state changes — needs an answer
    # to avoid EOF under set -e when running without a TTY.
    result = subprocess.run(
        cmd, cwd=CASE_DIR, env=env, input="N\n",
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    if result.returncode != 0 or not binary.exists():
        print(result.stdout[-2000:])
        sys.exit("Compilation failed.")

    print("Compilation successful.")
    return True


def modify_par(k_dim, schur_tgt=0, maxmodes=0):
    """Set userParam03/04/07 in .par file."""
    text = PAR_FILE.read_text()
    text = re.sub(r"userParam03\s*=\s*\d+", f"userParam03 = {schur_tgt}", text)
    text = re.sub(r"userParam04\s*=\s*\d+", f"userParam04 = {maxmodes}", text)
    text = re.sub(r"userParam07\s*=\s*\d+", f"userParam07 = {k_dim}", text)
    PAR_FILE.write_text(text)


def clean_outputs():
    """Remove previous run outputs (spectra, field files, logs)."""
    patterns = [
        "Spectre_*.dat", "Spectre_*.info",
        "dRe*.f?????", "dIm*.f?????", "dRe*.nek5000", "dIm*.nek5000",
        "KRY*", "HES*", "*.log.*", "logfile",
    ]
    for pat in patterns:
        for f in CASE_DIR.glob(pat):
            f.unlink(missing_ok=True)


def run_case(label, k_dim, schur_tgt, maxmodes, nprocs, dry_run=False):
    """Run a single eigensolver case and collect results."""
    cmd = ["mpirun", "-np", str(nprocs), "--bind-to", "core",
           str(CASE_DIR / "nek5000")]

    extra = ""
    if schur_tgt > 0:
        extra = f", schur_tgt={schur_tgt}, maxmodes={maxmodes}"

    if dry_run:
        print(f"[DRY-RUN] {label} (k_dim={k_dim}{extra}): {' '.join(cmd)}")
        return None

    out_dir = RESULTS_DIR / label
    out_dir.mkdir(parents=True, exist_ok=True)

    modify_par(k_dim, schur_tgt, maxmodes)
    clean_outputs()

    # Write SESSION.NAME
    session = CASE_DIR / "SESSION.NAME"
    session.write_text(f"{CASE_NAME}\n{CASE_DIR}/\n")

    print(f"\n=== Running {label} (k_dim={k_dim}{extra}, nprocs={nprocs}) ===")
    t0 = time.time()
    result = subprocess.run(
        cmd, cwd=CASE_DIR,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    wall_time = time.time() - t0

    # Save log
    log_file = out_dir / "run.log"
    log_file.write_text(result.stdout)

    if result.returncode != 0:
        print(f"  WARNING: {label} exited with code {result.returncode}")
        # Don't abort — partial results may still be useful

    # Copy output files
    for pat in ["Spectre_*.dat", "Spectre_*.info", "Spectre_*_conv.dat"]:
        for f in CASE_DIR.glob(pat):
            shutil.copy2(f, out_dir / f.name)

    print(f"  Completed in {wall_time:.1f}s")
    return wall_time


# ---------------------------------------------------------------------------
# Results parsing
# ---------------------------------------------------------------------------

def load_spectrum(path):
    """Load spectrum data (handles 2-col and 3-col formats).

    Returns (col1, col2, residuals).  residuals is None for 2-col files.
    """
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] >= 3:
        return data[:, 0], data[:, 1], data[:, 2]
    return data[:, 0], data[:, 1], None


def parse_info(path):
    """Parse Spectre_d.info file for eigensolver metadata."""
    info = {
        "k_dim": None,
        "schur_iterations": None,
        "schur_target": None,
        "nsteps": None,
        "dt": None,
    }
    if not path.exists():
        return info

    text = path.read_text()

    m = re.search(r"k_dim=\s*(\d+)", text)
    if m:
        info["k_dim"] = int(m.group(1))

    m = re.search(r"schur iterations=\s*(\d+)", text)
    if m:
        info["schur_iterations"] = int(m.group(1))

    m = re.search(r"schur_target=\s*(\d+)", text)
    if m:
        info["schur_target"] = int(m.group(1))

    m = re.search(r"nsteps=\s*(\d+)", text)
    if m:
        info["nsteps"] = int(m.group(1))

    m = re.search(r"dt=\s*([\d.Ee+-]+)", text)
    if m:
        info["dt"] = float(m.group(1))

    m = re.search(r"outp=\s*(\d+)", text)
    if m:
        info["outp"] = int(m.group(1))

    m = re.search(r"total matvecs=\s*(\d+)", text)
    if m:
        info["total_matvecs"] = int(m.group(1))

    return info


def estimate_matvecs(info):
    """Get total matvecs from info dict (exact if available, else estimate).

    Exact: read from 'total matvecs=' field (new nekStab format).
    Estimate for old format (upper bound):
      Pure Arnoldi (0 restarts): matvecs = k_dim.
      Krylov-Schur: matvecs ~ k_dim + restarts * (k_dim - schur_target).
      Note: uses schur_target as proxy for nsel (actual retained vectors).
      The real nsel is adaptive and typically larger, so this overestimates.
    """
    # Prefer exact value from new info format
    if info.get("total_matvecs") is not None:
        return info["total_matvecs"]

    k = info.get("k_dim")
    restarts = info.get("schur_iterations")
    nsel = info.get("schur_target")
    if k is None:
        return None
    if restarts is None or restarts == 0:
        return k
    if nsel is None:
        nsel = 2
    # Each restart resumes from nsel vectors, builds up to k
    return k + restarts * (k - nsel)


def collect_results(labels=None):
    """Collect results from all completed cases."""
    results = []
    for label, k_dim, schur_tgt, maxmodes in CASES:
        if labels and label not in labels:
            continue
        out_dir = RESULTS_DIR / label
        ns_file = out_dir / "Spectre_NSd.dat"
        h_file = out_dir / "Spectre_Hd.dat"
        info_file = out_dir / "Spectre_d.info"

        if not ns_file.exists():
            continue

        ns_sigma, ns_omega, ns_res = load_spectrum(ns_file)
        h_re, h_im, h_res = load_spectrum(h_file) if h_file.exists() else (None, None, None)
        info = parse_info(info_file)

        # Count converged eigenvalues
        n_converged = 0
        if ns_res is not None:
            n_converged = int(np.sum(ns_res < EIGENTOL))

        # Parse wall time from log
        wall_time = None
        log_file = out_dir / "run.log"
        if log_file.exists():
            log_text = log_file.read_text()
            m = re.search(r"total elapsed time\s*:\s*([\d.Ee+-]+)", log_text)
            if m:
                wall_time = float(m.group(1))

        total_matvecs = estimate_matvecs(info)

        results.append({
            "label": label,
            "k_dim": k_dim,
            "ns_sigma": ns_sigma,
            "ns_omega": ns_omega,
            "ns_res": ns_res,
            "h_re": h_re,
            "h_im": h_im,
            "h_res": h_res,
            "info": info,
            "restarts": info.get("schur_iterations"),
            "total_matvecs": total_matvecs,
            "n_converged": n_converged,
            "wall_time": wall_time,
        })

    return results


def include_reference_data():
    """Load pre-existing reference data from kdim*/ directories."""
    refs = []
    ref_map = {"kdim50": 50, "kdim92": 92, "kdim200": 200}
    for dirname, k_dim in ref_map.items():
        ref_dir = CASE_DIR / dirname
        ns_file = ref_dir / "Spectre_NSd.dat"
        h_file = ref_dir / "Spectre_Hd.dat"
        info_file = ref_dir / "Spectre_d.info"

        if not ns_file.exists():
            continue

        ns_sigma, ns_omega, ns_res = load_spectrum(ns_file)
        h_re, h_im, h_res = load_spectrum(h_file) if h_file.exists() else (None, None, None)
        info = parse_info(info_file)

        n_converged = 0
        if ns_res is not None:
            n_converged = int(np.sum(ns_res < EIGENTOL))

        refs.append({
            "label": dirname,
            "k_dim": k_dim,
            "ns_sigma": ns_sigma,
            "ns_omega": ns_omega,
            "ns_res": ns_res,
            "h_re": h_re,
            "h_im": h_im,
            "h_res": h_res,
            "info": info,
            "restarts": info.get("schur_iterations"),
            "total_matvecs": estimate_matvecs(info),
            "n_converged": n_converged,
            "wall_time": None,
        })

    return refs


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_converged(results):
    """Check that converged eigenvalues agree across cases and with reference."""
    print("\n=== Convergence Validation ===")

    # Extract leading eigenvalue (largest sigma) from each case
    leading = []
    for r in results:
        if r["ns_res"] is not None:
            mask = r["ns_res"] < EIGENTOL
            if np.any(mask):
                idx = np.argmax(r["ns_sigma"][mask])
                sigma = r["ns_sigma"][mask][idx]
                omega = r["ns_omega"][mask][idx]
                leading.append((r["label"], sigma, omega))

    if not leading:
        print("  No converged eigenvalues found in any case.")
        return False

    ok = True

    # Cross-case agreement
    if len(leading) > 1:
        sigmas = np.array([s for _, s, _ in leading])
        omegas = np.array([o for _, _, o in leading])
        sigma_spread = np.ptp(sigmas)
        omega_spread = np.ptp(omegas)
        print(f"  Leading eigenvalue across {len(leading)} cases:")
        for label, s, o in leading:
            print(f"    {label:>8s}: sigma = {s:+.7f}, omega = {o:.7f}, "
                  f"f = {o/(2*np.pi):.7f}")
        print(f"  Cross-case spread: sigma = {sigma_spread:.2e}, "
              f"omega = {omega_spread:.2e}")
        if sigma_spread > 1e-4 or omega_spread > 1e-4:
            print("  WARNING: Cross-case spread exceeds 1e-4")
            ok = False
        else:
            print("  OK: Cross-case agreement within 1e-4")

    # Marquet comparison (unstable mode: sigma>0 with omega near 0.735)
    ref_sigma = 0.0127451
    ref_omega = 0.735069
    for label, sigma, omega in leading:
        dsig = abs(sigma - ref_sigma)
        domg = abs(omega - ref_omega)
        status = "OK" if dsig < 0.005 and domg < 0.03 else "MISMATCH"
        print(f"  {label} vs Marquet: dsigma = {dsig:.5f}, "
              f"domega = {domg:.5f} [{status}]")
        if status == "MISMATCH":
            ok = False

    return ok


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_results(results, ref_data=None):
    """Generate 2-panel figure: H-plane and NS-plane."""
    params = {
        "text.usetex": False,
        "font.size": 8,
        "legend.fontsize": 7,
        "legend.handlelength": 1.2,
        "axes.linewidth": 0.5,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
    }
    plt.rcParams.update(params)

    fig, (ax_h, ax_ns) = plt.subplots(
        2, 1, figsize=(4.3, 6.8), constrained_layout=True
    )

    # --- Panel (a): H-plane ---
    theta = np.linspace(0, 2 * np.pi, 400)
    ax_h.plot(np.cos(theta), np.sin(theta), lw=0.5, color="tab:red",
              ls="-", zorder=0)
    ax_h.axhline(0, lw=0.2, color="k", ls=":")
    ax_h.axvline(0, lw=0.2, color="k", ls=":")

    for r in results:
        sty = STYLE.get(r["label"], {"color": "gray", "marker": "o"})
        if r["h_re"] is None:
            continue
        res = r["h_res"]
        for k in range(len(r["h_re"])):
            converged = res is not None and res[k] < EIGENTOL
            alpha = 1.0 if converged else 0.25
            fc = sty["color"] if converged else "none"
            ax_h.scatter(
                r["h_re"][k], r["h_im"][k],
                s=12 if converged else 6,
                alpha=alpha, marker=sty["marker"],
                facecolors=fc, edgecolors=sty["color"],
                linewidth=0.3, zorder=2 if converged else 1,
            )

    ax_h.set_xlim(-1.15, 1.15)
    ax_h.set_ylim(-1.15, 1.15)
    ax_h.set_aspect("equal")
    ax_h.set_xlabel(r"$\Re(\mu)$")
    ax_h.set_ylabel(r"$\Im(\mu)$")
    ax_h.set_title("(a) H-plane (Floquet multipliers)", fontsize=9, loc="left")

    # --- Panel (b): NS-plane ---
    ax_ns.axhline(0, lw=0.3, color="k", ls="--")
    ax_ns.axvline(0, lw=0.2, color="k", ls=":")

    # Plot Marquet reference first (background)
    ax_ns.scatter(
        MARQUET_OMEGA / (2 * np.pi), MARQUET_SIGMA,
        s=30, marker="+", color="gray", linewidths=0.8,
        label="Marquet et al.", zorder=3,
    )

    # Plot each case
    for r in results:
        sty = STYLE.get(r["label"], {"color": "gray", "marker": "o"})
        res = r["ns_res"]
        label_set = False
        for k in range(len(r["ns_sigma"])):
            converged = res is not None and res[k] < EIGENTOL
            alpha = 1.0 if converged else 0.25
            fc = sty["color"] if converged else "none"
            lbl = r["label"] if not label_set else None
            label_set = True
            ax_ns.scatter(
                r["ns_omega"][k] / (2 * np.pi), r["ns_sigma"][k],
                s=14 if converged else 6,
                alpha=alpha, marker=sty["marker"],
                facecolors=fc, edgecolors=sty["color"],
                linewidth=0.3, label=lbl, zorder=2 if converged else 1,
            )

    ax_ns.set_ylim(-0.15, 0.03)
    ax_ns.set_xlabel(r"$f = \omega / 2\pi$")
    ax_ns.set_ylabel(r"$\sigma$")
    ax_ns.set_title("(b) NS-plane (growth rate vs frequency)", fontsize=9,
                     loc="left")
    ax_ns.legend(loc="lower left", fontsize=6, ncol=2, framealpha=0.8)

    fig.savefig(OUTPUT_PNG, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"\nFigure saved: {OUTPUT_PNG}")


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

def print_summary(results):
    """Print a summary table of all cases."""
    print("\n" + "=" * 86)
    print(f"{'Label':>10s} {'k_dim':>6s} {'nev':>4s} {'Restarts':>9s} {'Matvecs':>8s} "
          f"{'Conv':>5s} {'Wall (s)':>9s} {'Regime':<14s}")
    print("-" * 86)
    for r in results:
        restarts = r["restarts"]
        if restarts is None:
            regime = "?"
        elif restarts == 0:
            regime = "Pure Arnoldi"
        else:
            regime = "Krylov-Schur"
        nev = r["info"].get("schur_target")
        nev_str = str(nev) if nev is not None else "?"
        rst_str = str(restarts) if restarts is not None else "?"
        mv_str = str(r["total_matvecs"]) if r["total_matvecs"] is not None else "?"
        wt_str = f"{r['wall_time']:.1f}" if r["wall_time"] else "—"
        print(f"{r['label']:>10s} {r['k_dim']:>6d} {nev_str:>4s} {rst_str:>9s} "
              f"{mv_str:>8s} {r['n_converged']:>5d} {wt_str:>9s} {regime:<14s}")
    print("=" * 86)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Eigensolver validation for cylinder at Re=50"
    )
    parser.add_argument("--plot-only", action="store_true",
                        help="Plot from existing results/ (skip compile+run)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show commands without executing")
    parser.add_argument("-n", type=int, default=None,
                        help="Override MPI process count")
    parser.add_argument("--cases", nargs="+", default=None,
                        help="Run subset of cases (e.g. k40 k200)")
    args = parser.parse_args()

    nekstab_root = detect_nekstab_root()
    nprocs = args.n if args.n else detect_nprocs(nekstab_root)

    case_labels = args.cases
    run_cases = [c for c in CASES if case_labels is None or c[0] in case_labels]

    if not args.plot_only:
        # Save original .par to restore later
        par_backup = PAR_FILE.read_text()

        try:
            compile_once(nekstab_root, dry_run=args.dry_run)

            for label, k_dim, schur_tgt, maxmodes in run_cases:
                run_case(label, k_dim, schur_tgt, maxmodes, nprocs,
                         dry_run=args.dry_run)
        finally:
            # Restore original .par
            PAR_FILE.write_text(par_backup)

    # Collect and analyse
    results = collect_results(labels=case_labels)

    # Also load pre-existing reference data from kdim*/ directories
    ref_data = include_reference_data()
    # Merge: use new results where available, fill gaps with reference
    existing_labels = {r["label"] for r in results}
    for ref in ref_data:
        # Map kdimXX -> kXX for style lookup
        mapped = "k" + str(ref["k_dim"])
        if mapped not in existing_labels:
            ref["label"] = mapped
            results.append(ref)

    if not results:
        sys.exit("No results found. Run without --plot-only first.")

    # Sort by k_dim
    results.sort(key=lambda r: r["k_dim"])

    print_summary(results)
    validate_converged(results)
    plot_results(results, ref_data=ref_data)


if __name__ == "__main__":
    main()
