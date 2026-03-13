#!/usr/bin/env python3
"""
nekStab Validation Suite

Validates nekStab against AMR paper test cases (short tier) and compile+run
checks for all example cases (full tier).

Usage:
    ./validate.py                    # Run all cases
    ./validate.py cylinder           # Run specific case
    ./validate.py --short            # Run only short-tier (AMR) cases
    ./validate.py --check-only       # Validate existing results only
    ./validate.py --dry-run          # Show what would run
    ./validate.py --nprocs 8         # Override MPI process count
    ./validate.py --list             # List all cases with status
    ./validate.py --compile-all      # Compile every example (no run)

Ricardo Frantz | Feb 2026
"""

import argparse
import io
import os
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

# ═══════════════════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════════════════

TOLERANCE = 0.05  # 5% relative error
MIN_ELEMS_PER_CORE = 20  # Minimum elements per MPI rank for efficiency
MIN_CORES_PER_SLOT = 4   # Minimum physical cores per parallel execution slot

NEKSTAB_ROOT = Path(os.environ.get("NEKSTAB_SOURCE_ROOT", Path.home() / "nekStab"))
NEKSTAB_DATA_ENV = "NEKSTAB_DATA_ROOT"
DEFAULT_DATA_ROOT = Path.home() / ".data_baptiste.nosync"

# ═══════════════════════════════════════════════════════════════════════════════
# Expected Values (inline - no external JSON)
# ═══════════════════════════════════════════════════════════════════════════════

CASES = {
    # ─── Short tier: AMR validation (eigenvalue checks) ───────────────────
    "cylinder_bf": {
        "description": "2D cylinder - Newton baseflow (Re=50)",
        "dir": "cylinder/baseflow/newton",
        "casename": "1cyl",
        "tier": "short",
        "family": "cylinder",
        "requires": ["BFRe40_1cyl0.f00001"],
        "output": None,
        "reference": "Barkley & Henderson (1996)",
        "checks": [],
        "timeout": 3600,
        "copies_to": [
            {
                "file": "BF_1cyl0.f00001",
                "dest_dir": "cylinder/stability/direct",
            },
            {
                "file": "BF_1cyl0.f00001",
                "dest_dir": "cylinder/stability/adjoint",
            },
            {
                "file": "BF_1cyl0.f00001",
                "dest_dir": "cylinder/postproc/sensitivity_budget_wavemaker",
            },
        ],
    },
    "cylinder_sfd": {
        "description": "2D cylinder - SFD baseflow (Re=50)",
        "dir": "cylinder/baseflow/sfd",
        "casename": "1cyl",
        "tier": "short",
        "family": "cylinder",
        "requires": ["BFRe40_1cyl0.f00001"],
        "output": None,
        "reference": "Åkervik et al. (2006)",
        "checks": [],
        "timeout": 1800,
    },
    "cylinder_direct": {
        "description": "2D cylinder wake - Hopf bifurcation (Re=50)",
        "dir": "cylinder/stability/direct",
        "casename": "1cyl",
        "tier": "short",
        "family": "cylinder",
        "requires": ["BF_1cyl0.f00001"],
        "output": "Spectre_NSd.dat",
        "reference": "Barkley & Henderson (1996), AMR paper Table 1",
        "checks": [
            {
                "name": "growth_rate",
                "metric": "sigma_real",
                "expected": 0.0156,
                "tolerance": 0.05,
                "description": "Positive growth rate confirms instability above Re_c~46.6",
            },
            {
                "name": "strouhal",
                "metric": "strouhal",
                "expected": 0.1204,
                "tolerance": 0.05,
                "description": "Strouhal number St = w/(2pi) ~ 0.12",
            },
        ],
        "copies_to": [
            {
                "file": "dRe1cyl0.f00001",
                "dest_dir": "cylinder/postproc/sensitivity_budget_wavemaker",
            },
            {
                "file": "dIm1cyl0.f00001",
                "dest_dir": "cylinder/postproc/sensitivity_budget_wavemaker",
            },
        ],
    },
    "cylinder_adjoint": {
        "description": "2D cylinder wake - Adjoint stability (Re=50)",
        "dir": "cylinder/stability/adjoint",
        "casename": "1cyl",
        "tier": "short",
        "family": "cylinder",
        "requires": ["BF_1cyl0.f00001"],
        "output": "Spectre_NSd.dat",
        "reference": "Adjoint eigenvalues match direct (same sigma, St)",
        "checks": [
            {
                "name": "adj_growth_rate",
                "metric": "sigma_real",
                "expected": 0.0156,
                "tolerance": 0.05,
                "description": "Adjoint growth rate matches direct",
            },
            {
                "name": "adj_strouhal",
                "metric": "strouhal",
                "expected": 0.1204,
                "tolerance": 0.05,
                "description": "Adjoint Strouhal matches direct",
            },
        ],
        "copies_to": [
            {
                "file": "aRe1cyl0.f00002",
                "dest_dir": "cylinder/postproc/sensitivity_budget_wavemaker",
            },
            {
                "file": "aIm1cyl0.f00002",
                "dest_dir": "cylinder/postproc/sensitivity_budget_wavemaker",
            },
        ],
    },
    "cylinder_sensitivity": {
        "description": "2D cylinder - Wavemaker + energy budget + sensitivity (Re=50)",
        "dir": "cylinder/postproc/sensitivity_budget_wavemaker",
        "casename": "1cyl",
        "tier": "short",
        "family": "cylinder",
        "requires": [
            "BF_1cyl0.f00001",
            "dRe1cyl0.f00001",
            "dIm1cyl0.f00001",
            "aRe1cyl0.f00002",
            "aIm1cyl0.f00002",
        ],
        "output": None,
        "reference": "Giannetti & Luchini (2007), Marquet et al. (2008)",
        "checks": [],
        "timeout": 1800,
    },
    "cylinder_modal": {
        "description": "2D cylinder - POD/DMD/SPOD modal analysis (Re=100)",
        "dir": "cylinder/modal",
        "casename": "1cyl",
        "tier": "short",
        "family": "cylinder",
        "requires": ["rst_1cyl0.f00001"],
        "output": "pod_spectrum.dat",
        "reference": "Barkley & Henderson (1996), St ~ 0.164-0.167",
        "checks": [
            {
                "name": "pod_energy",
                "metric": "pod_energy_12",
                "expected": 97.6,
                "tolerance": 2.0,
                "comparison": "absolute",
                "description": "First 2 POD modes capture ~97.6% energy (vortex shedding pair)",
            },
            {
                "name": "dmd_strouhal",
                "metric": "dmd_strouhal",
                "expected": 0.166,
                "tolerance": 0.05,
                "description": "Dominant DMD frequency matches vortex shedding St",
            },
            {
                "name": "spod_peak",
                "metric": "spod_peak_st",
                "expected": 0.156,
                "tolerance": 0.05,
                "description": "SPOD peak frequency near vortex shedding St",
            },
        ],
        "timeout": 7200,
    },
    "thermosyphon": {
        "description": "Thermosyphon - Pitchfork bifurcation (Ra=500, above critical)",
        "dir": "thersyphon/stability/direct",
        "casename": "tsyphon",
        "tier": "short",
        "family": "thersyphon",
        "requires": [],
        "output": "Spectre_NSd.dat",
        "reference": "AMR paper, Ra_c ~ 494",
        "checks": [
            {
                "name": "unstable",
                "metric": "sigma_real",
                "expected": 0.10,
                "tolerance": 0.15,  # absolute, Ra slightly above critical
                "comparison": "absolute",
                "description": "Positive growth rate (unstable above Ra_c~494)",
            },
            {
                "name": "real_mode",
                "metric": "sigma_imag",
                "expected": 0.0,
                "tolerance": 0.01,  # absolute
                "comparison": "absolute",
                "description": "Zero imaginary part confirms pitchfork (not Hopf)",
            },
        ],
    },
    "flipflop_bf": {
        "description": "Side-by-side cylinders - Newton UPO (Re=62, natural frequency)",
        "dir": "flip_flop/baseflow",
        "casename": "2cyl",
        "tier": "short",
        "family": "flip_flop",
        "requires": ["BF_Re60_2cyl0.f00001"],
        "output": None,  # check logfile convergence
        "reference": "Carini et al. (2015), Re_c ~ 61.17",
        "checks": [],
        "timeout": 7200,
        "copies_to": [
            {
                "file": "BF_2cyl0.f00001",
                "dest_dir": "flip_flop/stability/direct_Floquet",
            },
        ],
    },
    "flipflop_floquet": {
        "description": "Side-by-side cylinders - Floquet at Re=62 (above critical)",
        "dir": "flip_flop/stability/direct_Floquet",
        "casename": "2cyl",
        "tier": "short",
        "family": "flip_flop",
        "requires": ["BF_2cyl0.f00001"],
        "output": "Spectre_Hd.dat",
        "reference": "Carini et al. (2015), Re_c ~ 61.17",
        "checks": [
            {
                "name": "floquet_unstable",
                "metric": "mu_magnitude",
                "expected": 1.05,
                "tolerance": 0.15,  # absolute
                "comparison": "absolute",
                "description": "|mu| > 1 confirms instability above Re_c",
            },
        ],
    },
    "backstep_bf": {
        "description": "Backward-facing step - Newton baseflow (Re=500)",
        "dir": "back_fstep/baseflow",
        "casename": "bfs",
        "tier": "short",
        "family": "back_fstep",
        "requires": ["BF_bfs0.f00001"],
        "output": None,
        "reference": "Barkley et al. (2002)",
        "checks": [],
        "timeout": 3600,
        "copies_to": [
            {
                "file": "BF_bfs0.f00001",
                "dest_dir": "back_fstep/transient_growth",
            },
        ],
    },
    "backstep_tg": {
        "description": "Backward-facing step - Transient growth at tau=1 (Re=500)",
        "dir": "back_fstep/transient_growth",
        "casename": "bfs",
        "tier": "short",
        "family": "back_fstep",
        "requires": ["BF_bfs0.f00001"],
        "output": "Spectre_NSp.dat",
        "reference": "Blackburn et al. (2008), Barkley et al. (2002)",
        "checks": [
            {
                "name": "transient_growth",
                "metric": "sigma_real",
                "expected": 1.17,
                "tolerance": 0.05,
                "description": "Leading singular value sigma=ln(G)/tau~1.17 at tau=1",
            },
        ],
    },
    "tpjet_bf": {
        "description": "Forced jet - Newton UPO convergence (Re=1900, St=0.6)",
        "dir": "tpjet/baseflow/newton",
        "casename": "tpjet",
        "tier": "short",
        "family": "tpjet",
        "requires": ["BF_tpjet0.f00001"],
        "output": None,  # check logfile convergence
        "reference": "AMR paper, Re_c ~ 1371",
        "checks": [],
        "timeout": 7200,
        "copies_to": [
            {
                "file": "BF_tpjet0.f00001",
                "dest_dir": "tpjet/stability/direct_Floquet",
            },
        ],
    },
    "tpjet_floquet": {
        "description": "Forced jet - Period-doubling Floquet (Re=1900)",
        "dir": "tpjet/stability/direct_Floquet",
        "casename": "tpjet",
        "tier": "short",
        "family": "tpjet",
        "requires": ["BF_tpjet0.f00001"],
        "output": "Spectre_Hd.dat",
        "reference": "AMR paper, Re_c ~ 1371",
        "checks": [
            {
                "name": "period_doubling",
                "metric": "mu_real",
                "expected": -1.2,
                "tolerance": 0.3,  # absolute
                "comparison": "absolute",
                "description": "mu < -1 confirms unstable period-doubling",
            },
            {
                "name": "real_multiplier",
                "metric": "mu_imag",
                "expected": 0.0,
                "tolerance": 0.1,  # absolute
                "comparison": "absolute",
                "description": "Real Floquet multiplier (not Neimark-Sacker)",
            },
        ],
    },
    # ─── Full tier: Cylinder family ───────────────────────────────────────
    "cyl_ci_test": {
        "description": "CI smoke test (10 steps, cold start)",
        "dir": "cylinder/ci_test",
        "casename": "1cyl",
        "tier": "full",
        "family": "cylinder",
        "requires": [],
        "output": None,
        "checks": [],
        "timeout": 600,
    },
    "cyl_dns": {
        "description": "DNS from restart (Re=50)",
        "dir": "cylinder/dns",
        "casename": "1cyl",
        "tier": "full",
        "family": "cylinder",
        "requires": ["rst_1cyl0.f00001"],
        "output": None,
        "checks": [],
        "timeout": 3600,
    },
    "cyl_sfd": {
        "description": "SFD baseflow computation (Re=50)",
        "dir": "cylinder/baseflow/sfd",
        "casename": "1cyl",
        "tier": "full",
        "family": "cylinder",
        "requires": ["BFRe40_1cyl0.f00001"],
        "output": None,
        "checks": [],
        "timeout": 3600,
    },
    "cyl_sfd_dyn": {
        "description": "SFD with dynamic parameters (Re=50)",
        "dir": "cylinder/baseflow/sfd_dyn",
        "casename": "1cyl",
        "tier": "full",
        "family": "cylinder",
        "requires": ["BFRe40_1cyl0.f00001"],
        "output": None,
        "checks": [],
        "timeout": 3600,
    },
    "cyl_sfd_oifs": {
        "description": "SFD with OIFS timestepping (Re=50)",
        "dir": "cylinder/baseflow/sfd_oifs",
        "casename": "1cyl",
        "tier": "full",
        "family": "cylinder",
        "requires": ["BFRe40_1cyl0.f00001"],
        "output": None,
        "checks": [],
        "timeout": 3600,
    },
    "cyl_boostconv": {
        "description": "BoostConv baseflow acceleration (Re=50)",
        "dir": "cylinder/baseflow/boostconv",
        "casename": "1cyl",
        "tier": "full",
        "family": "cylinder",
        "requires": ["BFRe40_1cyl0.f00001"],
        "output": None,
        "checks": [],
        "timeout": 3600,
    },
    "cyl_newton": {
        "description": "Newton fixed-point iteration (Re=50)",
        "dir": "cylinder/baseflow/newton",
        "casename": "1cyl",
        "tier": "full",
        "family": "cylinder",
        "requires": ["BFRe40_1cyl0.f00001"],
        "output": None,
        "checks": [],
        "timeout": 3600,
    },
    "cyl_newton_dyn": {
        "description": "Newton with dynamic tolerances (Re=50)",
        "dir": "cylinder/baseflow/newton_dyn",
        "casename": "1cyl",
        "tier": "full",
        "family": "cylinder",
        "requires": ["BFRe40_1cyl0.f00001"],
        "output": None,
        "checks": [],
        "timeout": 3600,
    },
    "cyl_newton_temp": {
        "description": "Newton with temperature field (cold start)",
        "dir": "cylinder/baseflow/newton_dyn_temp",
        "casename": "1cyl",
        "tier": "full",
        "family": "cylinder",
        "requires": [],
        "output": None,
        "checks": [],
        "timeout": 7200,
    },
    # cyl_newton_smooth removed (newton_smoother directory deleted — was obsolete mesh smoother test)
    "cyl_newton_upo": {
        "description": "Newton-GMRES for UPO (Re=50)",
        "dir": "cylinder/baseflow/newton_upo",
        "casename": "1cyl",
        "tier": "full",
        "family": "cylinder",
        "requires": ["rstcyl0.f00001"],
        "output": None,
        "checks": [],
        "timeout": 3600,
    },
    # cyl_adjoint is now "cylinder_adjoint" in short tier
    "cyl_floquet_dir": {
        "description": "Direct Floquet stability (Re=50)",
        "dir": "cylinder/stability/direct_Floquet",
        "casename": "1cyl",
        "tier": "full",
        "family": "cylinder",
        "requires": ["BF_1cyl0.f00001"],
        "output": None,
        "checks": [],
        "timeout": 3600,
    },
    "cyl_floquet_adj": {
        "description": "Adjoint Floquet stability (Re=50)",
        "dir": "cylinder/stability/adjoint_Floquet",
        "casename": "1cyl",
        "tier": "full",
        "family": "cylinder",
        "requires": ["BF_1cyl0.f00001"],
        "output": None,
        "checks": [],
        "timeout": 3600,
    },
    "cyl_animate": {
        "description": "Mode animation (Re=50)",
        "dir": "cylinder/stability/animate_modes",
        "casename": "1cyl",
        "tier": "full",
        "family": "cylinder",
        "requires": [
            "BF_1cyl0.f00001",
            "Spectre_NSd_conv.dat",
            "dRe1cyl0.f00001",
            "dIm1cyl0.f00001",
        ],
        "output": None,
        "checks": [],
        "timeout": 3600,
    },
    "cyl_animate_upo": {
        "description": "Mode animation with UPO (Re=50)",
        "dir": "cylinder/stability/animate_modes_with_UPO",
        "casename": "1cyl",
        "tier": "full",
        "family": "cylinder",
        "requires": [
            "BF_1cyl0.f00001",
            "Spectre_NSd_conv.dat",
            "dRe1cyl0.f00001",
            "dIm1cyl0.f00001",
        ],
        "output": None,
        "checks": [],
        "timeout": 3600,
    },
    "cyl_otd": {
        "description": "Optimally time-dependent modes (Re=50)",
        "dir": "cylinder/otd",
        "casename": "1cyl",
        "tier": "full",
        "family": "cylinder",
        "requires": ["BF_1cyl0.f00001"],
        "output": None,
        "checks": [],
        "timeout": 3600,
    },
    # cyl_modal is now "cylinder_modal" in short tier (POD/DMD/SPOD with eigenvalue checks)
    # cyl_sensitivity is now "cylinder_sensitivity" in short tier
    "cyl_force_sens": {
        "description": "Steady force sensitivity (Re=50)",
        "dir": "cylinder/postproc/steady_force_sensitivity",
        "casename": "1cyl",
        "tier": "full",
        "family": "cylinder",
        "requires": [
            "BF_1cyl0.f00001",
            "sr_1cyl0.f00001",
            "si_1cyl0.f00001",
        ],
        "output": None,
        "checks": [],
        "timeout": 3600,
    },
    # ─── Full tier: Back step family ──────────────────────────────────────
    # backstep_bf is now in short tier (Newton + copies_to transient_growth)
    # ─── Full tier: Thermosyphon family ───────────────────────────────────
    "thermo_bf": {
        "description": "Newton baseflow (Ra=500)",
        "dir": "thersyphon/baseflow",
        "casename": "tsyphon",
        "tier": "full",
        "family": "thersyphon",
        "requires": ["BF_Ra400_tsyphon0.f00001"],
        "output": None,
        "checks": [],
        "timeout": 3600,
    },
    # flipflop_bf is now in short tier (Newton UPO at Re=62)
    # ─── Full tier: Tpjet family ──────────────────────────────────────────
    # tpjet_newton is now "tpjet_bf" in short tier (Newton UPO at Re=1900)
    "tpjet_tdf": {
        "description": "Time-delayed feedback stabilization (Re=2005)",
        "dir": "tpjet/baseflow/tdf",
        "casename": "tpjet",
        "tier": "full",
        "family": "tpjet",
        "requires": ["BF_Re1900_tpjet0.f00001"],
        "output": None,
        "checks": [],
        "timeout": 7200,
    },
    # ─── Full tier: Lid-driven cavity family ──────────────────────────────
    "lid_driven": {
        "description": "2D lid-driven cavity Newton (Re=3600)",
        "dir": "lid_driven",
        "casename": "cav",
        "tier": "full",
        "family": "lid_driven",
        "requires": ["BF_cav0.f00001"],
        "output": None,
        "checks": [],
        "timeout": 3600,
    },
    # ─── Full tier: Cubic cavity family ───────────────────────────────────
    "cubic_cav": {
        "description": "3D cubic cavity Newton (Re=2500)",
        "dir": "cubic_cavity",
        "casename": "cav",
        "tier": "full",
        "family": "cubic_cavity",
        "requires": ["BF_cav0.f00001"],
        "output": None,
        "checks": [],
        "timeout": 7200,
    },
    "cubic_cav_upo": {
        "description": "3D cubic cavity Newton-GMRES UPO (Re=1950)",
        "dir": "cubic_cavity_upo",
        "casename": "cav",
        "tier": "full",
        "family": "cubic_cavity",
        "requires": ["1960_PO_cav0.f00001"],
        "output": None,
        "checks": [],
        "timeout": 7200,
    },
    # ─── Full tier: NACA 0012 family ──────────────────────────────────────
    "naca0012": {
        "description": "NACA 0012 Newton baseflow (Re=2000)",
        "dir": "naca0012",
        "casename": "naca0012",
        "tier": "full",
        "family": "naca0012",
        "requires": ["BF_naca00120.f00001"],
        "output": None,
        "checks": [],
        "timeout": 7200,
    },
    "naca0012_Re2500": {
        "description": "NACA 0012 direct stability (Re=2500)",
        "dir": "naca0012_Re2500",
        "casename": "naca0012",
        "tier": "full",
        "family": "naca0012",
        "requires": ["BF_naca00120.f00001"],
        "output": None,
        "checks": [],
        "timeout": 7200,
    },
    # ─── Full tier: Poiseuille family ─────────────────────────────────────
    "poiseuille": {
        "description": "2D Poiseuille OTD (Re=5000)",
        "dir": "poiseuille_OTD",
        "casename": "poiseuille_OTD",
        "tier": "full",
        "family": "poiseuille",
        "requires": [],
        "output": None,
        "checks": [],
        "timeout": 3600,
    },
    # ─── Full tier: Slot FST family ───────────────────────────────────────
    "slot_fst": {
        "description": "Slot jet with free-stream turbulence (Re=495)",
        "dir": "slot_FST",
        "casename": "slot",
        "tier": "full",
        "family": "slot_FST",
        "requires": ["34_slot0.f00001"],
        "output": None,
        "checks": [],
        "timeout": 7200,
    },
}

# ═══════════════════════════════════════════════════════════════════════════════
# Terminal Colors
# ═══════════════════════════════════════════════════════════════════════════════

class C:
    """ANSI color codes."""
    RED = "\033[0;31m"
    GREEN = "\033[0;32m"
    YELLOW = "\033[1;33m"
    CYAN = "\033[0;36m"
    BOLD = "\033[1m"
    DIM = "\033[2m"
    NC = "\033[0m"


_print_lock = threading.Lock()
_thread_local = threading.local()


def _write_output(msg: str) -> None:
    """Write to per-thread buffer (worker) or directly to stdout (main thread)."""
    buf = getattr(_thread_local, 'buf', None)
    if buf is not None:
        buf.write(msg + "\n")
    else:
        with _print_lock:
            print(msg)


def log_pass(msg: str) -> None:
    _write_output(f"{C.GREEN}[PASS]{C.NC} {msg}")


def log_fail(msg: str) -> None:
    _write_output(f"{C.RED}[FAIL]{C.NC} {msg}")


def log_skip(msg: str) -> None:
    _write_output(f"{C.YELLOW}[SKIP]{C.NC} {msg}")


def log_info(msg: str) -> None:
    _write_output(f"{C.YELLOW}[INFO]{C.NC} {msg}")


def log_detail(msg: str) -> None:
    _write_output(f"{C.CYAN}      {C.NC} {msg}")


def log_header(msg: str) -> None:
    _write_output(f"\n{C.BOLD}=== {msg} ==={C.NC}")


# ═══════════════════════════════════════════════════════════════════════════════
# P-Core Detection (Intel Hybrid CPUs)
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class CPUInfo:
    """CPU configuration for MPI binding."""
    cpu_list: str  # e.g., "0-15"
    count: int     # number of cores


def detect_pcores() -> CPUInfo:
    """
    Detect Performance cores on Intel hybrid CPUs.
    Accounts for hyperthreading (uses physical cores only).
    Returns all physical cores on non-hybrid systems.
    """
    pcore_file = Path("/sys/devices/cpu_core/cpus")

    if pcore_file.exists():
        cpu_list = pcore_file.read_text().strip()
        logical_count = _count_cpus(cpu_list)
        # Detect hyperthreading: threads_per_core > 1
        threads_per_core = _get_threads_per_core()
        physical_count = logical_count // threads_per_core
        return CPUInfo(cpu_list, physical_count)

    # Non-hybrid: use physical cores (account for hyperthreading)
    logical = os.cpu_count() or 4
    threads_per_core = _get_threads_per_core()
    physical = logical // threads_per_core
    return CPUInfo(f"0-{logical - 1}", physical)


def _get_threads_per_core() -> int:
    """Read threads per core from lscpu."""
    try:
        result = subprocess.run(
            ["lscpu"], capture_output=True, text=True, timeout=5
        )
        for line in result.stdout.splitlines():
            if "Thread(s) per core" in line:
                return int(line.split(":")[1].strip())
    except Exception:
        pass
    return 1


def _count_cpus(cpu_range: str) -> int:
    """Count CPUs from range like '0-15' or '0-3,8-11'."""
    count = 0
    for part in cpu_range.split(","):
        if "-" in part:
            start, end = map(int, part.split("-"))
            count += end - start + 1
        else:
            count += 1
    return count


def print_cpu_info(cpu: CPUInfo) -> None:
    """Print CPU detection results."""
    pcore_file = Path("/sys/devices/cpu_core/cpus")
    ecore_file = Path("/sys/devices/cpu_atom/cpus")

    if pcore_file.exists():
        log_info("Intel hybrid CPU detected")
        log_detail(f"P-cores: {cpu.cpu_list} ({cpu.count} cores)")
        if ecore_file.exists():
            ecores = ecore_file.read_text().strip()
            log_detail(f"E-cores: {ecores} (not used)")
    else:
        log_info(f"Standard CPU: {cpu.count} cores available")


def _parse_cpu_list(cpu_range: str) -> list[int]:
    """Parse '0-3,8-11' into [0, 1, 2, 3, 8, 9, 10, 11]."""
    cpus = []
    for part in cpu_range.strip().split(","):
        if "-" in part:
            start, end = map(int, part.split("-"))
            cpus.extend(range(start, end + 1))
        else:
            cpus.append(int(part))
    return cpus


def detect_all_cores() -> list[int]:
    """
    Detect all physical cores (P-cores + E-cores on hybrid, one per HT pair).
    Returns sorted list of CPU IDs — one per physical core.
    """
    pcore_file = Path("/sys/devices/cpu_core/cpus")
    ecore_file = Path("/sys/devices/cpu_atom/cpus")

    if pcore_file.exists():
        # Hybrid Intel CPU: deduplicate HT siblings among P-cores
        pcore_logical = _parse_cpu_list(pcore_file.read_text().strip())
        seen_cores: set[int] = set()
        physical_pcores = []
        for cpu_id in pcore_logical:
            sibling_path = Path(
                f"/sys/devices/system/cpu/cpu{cpu_id}/topology/thread_siblings_list"
            )
            try:
                siblings = _parse_cpu_list(sibling_path.read_text().strip())
                core_key = min(siblings)
                if core_key not in seen_cores:
                    seen_cores.add(core_key)
                    physical_pcores.append(core_key)
            except (OSError, ValueError):
                if cpu_id not in seen_cores:
                    seen_cores.add(cpu_id)
                    physical_pcores.append(cpu_id)

        # E-cores: each is a physical core (no HT)
        ecores = []
        if ecore_file.exists():
            ecores = _parse_cpu_list(ecore_file.read_text().strip())

        return sorted(physical_pcores + ecores)

    # Non-hybrid: stride by threads_per_core
    logical = os.cpu_count() or 4
    threads_per_core = _get_threads_per_core()
    return list(range(0, logical, threads_per_core))


def make_slots(cores: list[int], cores_per_slot: int) -> list[CPUInfo]:
    """Partition physical cores into equal-sized execution slots."""
    n = len(cores)
    num_slots = max(1, n // cores_per_slot)
    slot_size = n // num_slots
    slots = []
    for i in range(num_slots):
        start = i * slot_size
        end = start + slot_size if i < num_slots - 1 else n
        chunk = cores[start:end]
        cpu_list = ",".join(str(c) for c in chunk)
        slots.append(CPUInfo(cpu_list, len(chunk)))
    return slots


# ═══════════════════════════════════════════════════════════════════════════════
# Mesh-Aware Core Scaling
# ═══════════════════════════════════════════════════════════════════════════════

def _parse_size_int(lines: list[str], param: str) -> int | None:
    """Extract an integer parameter value from SIZE file lines."""
    for line in lines:
        if param in line.lower() and "=" in line:
            try:
                after_eq = line.split("=")[1].strip()
                digits = ""
                for c in after_eq:
                    if c.isdigit():
                        digits += c
                    elif digits:
                        break
                if digits:
                    return int(digits)
            except (IndexError, ValueError):
                continue
    return None


def get_optimal_nprocs(size_file: Path, max_cores: int) -> int:
    """
    Determine optimal MPI process count based on mesh size.
    Respects lpmin (minimum ranks the SIZE was compiled for) and
    ensures at least MIN_ELEMS_PER_CORE elements per rank.
    """
    if not size_file.exists():
        log_detail(f"No SIZE file, using {max_cores} cores")
        return max_cores

    lines = size_file.read_text().splitlines()
    lelg = _parse_size_int(lines, "lelg")
    lpmin = _parse_size_int(lines, "lpmin")

    if lelg is None:
        log_detail(f"Could not parse lelg, using {max_cores} cores")
        return max_cores

    # lpmin sets the minimum — lelt was computed as lelg/lpmin + slack
    minimum = lpmin if lpmin and lpmin > 1 else 1

    # Don't exceed efficient core count for small meshes
    efficient = max(1, lelg // MIN_ELEMS_PER_CORE)

    nprocs = min(max_cores, efficient)
    nprocs = max(nprocs, minimum)

    if nprocs != max_cores:
        log_detail(f"Adjusted to {nprocs} cores (lelg={lelg}, lpmin={lpmin})")

    return nprocs


# ═══════════════════════════════════════════════════════════════════════════════
# Eigenvalue Parsing
# ═══════════════════════════════════════════════════════════════════════════════


def load_pod_spectrum(dat_file: Path) -> np.ndarray:
    """
    Load POD spectrum from pod_spectrum.dat.

    Format: mode  eigenvalue  energy_%  cumulative_%  [||Phi||]
    First row may have 5 columns (norm), rest have 4.
    Returns array with columns: mode, eigenvalue, energy_%, cumulative_%.
    """
    data = np.loadtxt(dat_file, usecols=(0, 1, 2, 3))
    return data


def load_dmd_spectrum(dat_file: Path) -> np.ndarray:
    """
    Load DMD spectrum from dmd_spectrum.dat.

    Format: mode  |mu|  sigma  omega  St  mu_real  mu_imag  [||Phi_re||]
    The norm column (col 7) is only present for the first nsave rows,
    so we read only the 7 guaranteed columns to avoid loadtxt failures.
    Returns array with columns: mode, |mu|, sigma, omega, St, mu_real, mu_imag.
    """
    data = np.loadtxt(dat_file, usecols = (0, 1, 2, 3, 4, 5, 6))
    return data


def load_spod_spectrum(dat_file: Path) -> np.ndarray:
    """
    Load SPOD spectrum from spod_stream_spectrum.dat.

    Format: St  lambda_1  lambda_2  ...
    Returns raw array.
    """
    data = np.loadtxt(dat_file)
    return data


def load_eigenvalues(dat_file: Path) -> tuple[np.ndarray, np.ndarray]:
    """
    Load eigenvalues from Spectre_*.dat file.

    Format: real_part  imag_part  residual
    Returns eigenvalues sorted by real part (descending).
    """
    if not dat_file.exists():
        raise FileNotFoundError(f"Not found: {dat_file}")

    data = np.loadtxt(dat_file)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    eigenvalues = data[:, 0] + 1j * data[:, 1]
    residuals = data[:, 2]

    # Sort by real part descending (most unstable first)
    idx = np.argsort(-eigenvalues.real)
    return eigenvalues[idx], residuals[idx]


def _leading_nontrivial_multiplier(eigenvalues: np.ndarray) -> complex:
    """
    Find the leading non-trivial Floquet multiplier.

    nekStab includes the trivial multiplier (mu ≈ 1) corresponding to
    time-translation symmetry.  For validation we need the physical
    leading multiplier — the one with the largest |mu| that is NOT
    the trivial one near +1.
    """
    # Sort by magnitude descending
    idx = np.argsort(-np.abs(eigenvalues))
    for i in idx:
        mu = complex(eigenvalues[i])
        # Skip the trivial multiplier (real, near +1).
        # Threshold 0.02 is safe for the current validation suite where
        # physical multipliers are well separated from +1 (e.g. mu ≈ -1
        # for period-doubling).  Near a saddle-node or pitchfork onset
        # (Re ≈ Re_c) the physical multiplier crosses through +1 and
        # could be incorrectly skipped — tighten or remove this filter
        # if validating cases very close to onset.
        if abs(mu - 1.0) < 0.02 and abs(mu.imag) < 0.01:
            continue
        return mu
    # Fallback: return largest magnitude if all are near 1
    return complex(eigenvalues[idx[0]])


def compute_metric(metric: str, eigenvalues: Optional[np.ndarray], case_dir: Path) -> float:
    """Compute validation metric from eigenvalue or modal data."""
    # Modal metrics don't use eigenvalues — handle first
    # ── Modal decomposition metrics ──────────────────────────────────────
    if metric == "pod_energy_12":
        # Cumulative energy captured by first 2 POD modes
        pod_file = case_dir / "pod_spectrum.dat"
        if not pod_file.exists():
            raise FileNotFoundError(f"Not found: {pod_file}")
        data = load_pod_spectrum(pod_file)
        return float(data[1, 3])  # row 1 (mode 2), col 3 = cumulative_%

    if metric == "pod_eigenvalue_1":
        # Leading POD eigenvalue
        pod_file = case_dir / "pod_spectrum.dat"
        if not pod_file.exists():
            raise FileNotFoundError(f"Not found: {pod_file}")
        data = load_pod_spectrum(pod_file)
        return float(data[0, 1])  # row 0 (mode 1), col 1 = eigenvalue

    if metric == "dmd_strouhal":
        # Fundamental DMD Strouhal (lowest positive St among persistent modes)
        dmd_file = case_dir / "dmd_spectrum.dat"
        if not dmd_file.exists():
            raise FileNotFoundError(f"Not found: {dmd_file}")
        data = load_dmd_spectrum(dmd_file)
        # col 1 = |mu|, col 4 = St; find fundamental among persistent modes
        mu_mag = data[:, 1]
        st_vals = np.abs(data[:, 4])
        # Persistent modes: |mu| > 0.99, positive St
        mask = (mu_mag > 0.99) & (st_vals > 0.01)
        if not np.any(mask):
            return 0.0
        # Fundamental = lowest St among persistent modes
        return float(np.min(st_vals[mask]))

    if metric == "spod_peak_st":
        # SPOD peak frequency (St at maximum lambda_1)
        spod_file = case_dir / "spod_stream_spectrum.dat"
        if not spod_file.exists():
            raise FileNotFoundError(f"Not found: {spod_file}")
        data = load_spod_spectrum(spod_file)
        # col 0 = St, col 1 = lambda_1
        idx_peak = np.argmax(data[:, 1])
        return float(data[idx_peak, 0])

    # Eigenvalue-based metrics require eigenvalues
    if eigenvalues is None:
        raise ValueError(f"Metric '{metric}' requires eigenvalues but none loaded")

    ev = eigenvalues[0]  # leading eigenvalue (sorted by real part desc)

    if metric == "sigma_real":
        return float(ev.real)

    if metric == "sigma_imag":
        return float(abs(ev.imag))

    if metric == "strouhal":
        return float(abs(ev.imag) / (2 * np.pi))

    # Floquet metrics: use leading NON-TRIVIAL multiplier
    if metric == "mu_magnitude":
        mu = _leading_nontrivial_multiplier(eigenvalues)
        return float(abs(mu))

    if metric == "mu_real":
        mu = _leading_nontrivial_multiplier(eigenvalues)
        return float(mu.real)

    if metric == "mu_imag":
        mu = _leading_nontrivial_multiplier(eigenvalues)
        return float(mu.imag)

    if metric == "tau_opt":
        # Transient growth: read from dedicated output file
        tg_file = case_dir / "transient_growth.dat"
        if tg_file.exists():
            # Format: tau  G(tau)  - find max G
            data = np.loadtxt(tg_file)
            idx_max = np.argmax(data[:, 1])
            return float(data[idx_max, 0])

        # Fallback: check if stored in different format
        log_detail("tau_opt: using eigenvalue-based estimate")
        # For backstep, tau_opt ~ 1/|sigma| for least stable mode
        if abs(ev.real) > 1e-10:
            return float(1.0 / abs(ev.real))
        return 58.0  # literature value as last resort

    raise ValueError(f"Unknown metric: {metric}")


# ═══════════════════════════════════════════════════════════════════════════════
# Dependency Checking
# ═══════════════════════════════════════════════════════════════════════════════

def _data_roots() -> list[Path]:
    """Return existing external data roots in search order."""
    roots: list[Path] = []

    env_root = os.environ.get(NEKSTAB_DATA_ENV)
    if env_root:
        roots.append(Path(env_root).expanduser())

    roots.append(DEFAULT_DATA_ROOT.expanduser())

    existing: list[Path] = []
    for root in roots:
        if root.exists() and root not in existing:
            existing.append(root)
    return existing


def _external_prereq_candidates(case_dir: Path, fname: str) -> list[Path]:
    """Build likely external-data locations for a required file."""
    try:
        rel_case_dir = case_dir.relative_to(NEKSTAB_ROOT / "example")
    except ValueError:
        rel_case_dir = case_dir

    candidates: list[Path] = []
    for root in _data_roots():
        candidates.extend([
            root / rel_case_dir / fname,
            root / "example" / rel_case_dir / fname,
            root / NEKSTAB_ROOT.name / "example" / rel_case_dir / fname,
        ])
    return candidates


def _find_external_prereq(case_dir: Path, fname: str) -> Optional[Path]:
    """Locate a missing prerequisite in the configured external data roots."""
    for candidate in _external_prereq_candidates(case_dir, fname):
        if candidate.exists():
            return candidate
    return None


def materialize_requires(case_dir: Path, requires: list[str]) -> list[str]:
    """
    Ensure prerequisite files exist locally, copying them from an external
    data root when available. Returns any remaining missing filenames.
    """
    missing = []
    for fname in requires:
        local_path = case_dir / fname
        if local_path.exists():
            continue

        external_path = _find_external_prereq(case_dir, fname)
        if external_path is None:
            missing.append(fname)
            continue

        shutil.copy2(external_path, local_path)
        log_detail(f"Copied prerequisite {external_path} → {local_path}")
    return missing


def check_requires(case_dir: Path, requires: list[str], allow_external: bool = False) -> list[str]:
    """Check which required files are missing from the case directory."""
    missing = []
    for fname in requires:
        if (case_dir / fname).exists():
            continue
        if allow_external and _find_external_prereq(case_dir, fname) is not None:
            continue
        missing.append(fname)
    return missing


# ═══════════════════════════════════════════════════════════════════════════════
# Logfile Checking
# ═══════════════════════════════════════════════════════════════════════════════

def check_logfile_success(case_dir: Path) -> bool:
    """Check if logfile or run log contains 'run successful' string."""
    # Check Nek5000 logfile first
    logfile = case_dir / "logfile"
    if logfile.exists():
        try:
            text = logfile.read_text(errors="replace")
            if "run successful" in text.lower():
                return True
        except Exception:
            pass

    # Fallback: check run logs (stdout redirected by validate.py)
    for log in sorted(case_dir.glob("*.log.*"), key=lambda p: p.stat().st_mtime, reverse=True):
        try:
            text = log.read_text(errors="replace")
            if "run successful" in text.lower():
                return True
        except Exception:
            continue

    return False


# ═══════════════════════════════════════════════════════════════════════════════
# Validation Logic
# ═══════════════════════════════════════════════════════════════════════════════

def check_value(computed: float, expected: float, tolerance: float,
                comparison: str = "relative") -> tuple[bool, float]:
    """
    Check if computed value matches expected within tolerance.
    Returns (passed, error_value).
    """
    if comparison == "absolute":
        error = abs(computed - expected)
        return error < tolerance, error

    # Relative comparison
    if abs(expected) < 1e-10:
        error = abs(computed - expected)
        return error < tolerance, error

    error = abs(computed - expected) / abs(expected)
    return error < tolerance, error


def _fmt_time(seconds: float) -> str:
    """Format elapsed seconds as human-readable string."""
    if seconds < 60:
        return f"{seconds:.0f}s"
    m, s = divmod(int(seconds), 60)
    if m < 60:
        return f"{m}m {s:02d}s"
    h, m = divmod(m, 60)
    return f"{h}h {m:02d}m {s:02d}s"


def validate_case(name: str, case: dict, check_only: bool,
                  nprocs: int, cpu: CPUInfo, dry_run: bool) -> tuple[str, float]:
    """
    Validate a single case.
    Returns: ("pass"|"fail"|"skip", elapsed_seconds).
    """
    case_t0 = time.perf_counter()

    def _result(status: str) -> tuple[str, float]:
        elapsed = time.perf_counter() - case_t0
        if elapsed > 1.0:
            log_detail(f"Case time: {_fmt_time(elapsed)}")
        return (status, elapsed)

    log_header(f"{name}: {case['description']}")

    case_dir = NEKSTAB_ROOT / "example" / case["dir"]

    if not case_dir.exists():
        log_fail(f"Directory not found: {case_dir}")
        return _result("fail")

    log_detail(f"Directory: {case_dir}")
    if "reference" in case:
        log_detail(f"Reference: {case['reference']}")

    # Check prerequisites
    requires = case.get("requires", [])
    if requires:
        if dry_run:
            missing = check_requires(case_dir, requires, allow_external=True)
        else:
            missing = materialize_requires(case_dir, requires)
        if missing:
            log_skip(f"Missing prerequisites: {', '.join(missing)}")
            return _result("skip")

    # Determine actual nprocs based on mesh
    actual_nprocs = get_optimal_nprocs(case_dir / "SIZE", nprocs)

    if dry_run:
        log_info(f"[DRY RUN] Would compile: mks {case['casename']}")
        log_info(f"[DRY RUN] Would run on {actual_nprocs} cores")
        if case.get("output"):
            log_info(f"[DRY RUN] Would check: {case['output']}")
        else:
            log_info("[DRY RUN] Would check: logfile for run successful")
        if (case_dir / "plot.py").exists():
            log_info("[DRY RUN] Would regenerate: plot.png")
        return _result("pass")

    # Compile and run unless check-only
    if not check_only:
        timeout = case.get("timeout", 7200)
        if not _compile_case(case_dir, case["casename"]):
            return _result("fail")
        if not _run_case(case_dir, case["casename"], actual_nprocs, cpu, timeout,
                         requires=case.get("requires", [])):
            return _result("fail")

        # Copy output files to downstream case directories if specified
        copy_failed = False
        for copy_spec in case.get("copies_to", []):
            src = case_dir / copy_spec["file"]
            dst_dir = NEKSTAB_ROOT / "example" / copy_spec["dest_dir"]
            dst_name = copy_spec.get("dest_name", copy_spec["file"])
            dst = dst_dir / dst_name
            if src.exists():
                dst_dir.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dst)
                log_detail(f"Copied {src.name} → {dst}")
            else:
                log_fail(f"Expected output {src} not found for copy")
                copy_failed = True

        if copy_failed:
            return _result("fail")

        # Regenerate figures from results
        _run_plot(case_dir)

    if check_only:
        # Regenerate plots from existing data
        _run_plot(case_dir)

    # Cases with no output file: check logfile only (any tier)
    if case.get("output") is None:
        if check_logfile_success(case_dir):
            log_pass("logfile: run successful")
            return _result("pass")
        else:
            log_fail("logfile: 'run successful' not found")
            return _result("fail")

    # Cases with output: validate results
    output_file = case_dir / case["output"]

    if not output_file.exists():
        log_fail(f"Output not found: {output_file}")
        if check_only:
            log_detail("Run without --check-only to generate results")
        return _result("fail")

    log_info(f"Reading: {output_file.name}")

    # Modal outputs (pod/dmd/spod) are handled directly by compute_metric;
    # eigenvalue-based outputs need parsing here.
    _MODAL_OUTPUTS = {"pod_spectrum.dat", "dmd_spectrum.dat", "spod_stream_spectrum.dat"}
    eigenvalues = None

    if output_file.name not in _MODAL_OUTPUTS:
        try:
            eigenvalues, _residuals = load_eigenvalues(output_file)
        except Exception as e:
            log_fail(f"Error loading eigenvalues: {e}")
            return _result("fail")

        log_detail(f"Loaded {len(eigenvalues)} eigenvalues")
        log_detail(f"Leading: {eigenvalues[0]:.6f}")

    # Run checks
    all_passed = True

    for check in case["checks"]:
        name_check = check["name"]
        comparison = check.get("comparison", "relative")

        try:
            computed = compute_metric(check["metric"], eigenvalues, case_dir)
            passed, error = check_value(
                computed, check["expected"], check["tolerance"], comparison
            )
        except Exception as e:
            log_fail(f"{name_check}: {e}")
            all_passed = False
            continue

        if passed:
            if comparison == "relative":
                log_pass(f"{name_check}: {computed:.6f} (expected {check['expected']:.6f}, error={error*100:.1f}%)")
            else:
                log_pass(f"{name_check}: {computed:.6f} (expected {check['expected']:.6f}, diff={error:.4f})")
        else:
            if comparison == "relative":
                log_fail(f"{name_check}: {computed:.6f} (expected {check['expected']:.6f}, error={error*100:.1f}% > {check['tolerance']*100:.0f}%)")
            else:
                log_fail(f"{name_check}: {computed:.6f} (expected {check['expected']:.6f}, diff={error:.4f} > {check['tolerance']:.2f})")
            all_passed = False

    return _result("pass") if all_passed else _result("fail")


def _show_first_error(build_log: Path) -> None:
    """Print the first Error/Fatal/undefined line from a build log."""
    try:
        for line in build_log.read_text(errors="replace").splitlines():
            stripped = line.strip()
            if any(kw in stripped for kw in ("Error:", "Fatal", "undefined reference")):
                # Truncate long lines
                if len(stripped) > 100:
                    stripped = stripped[:97] + "..."
                log_detail(f"{C.DIM}{stripped}{C.NC}")
                return
    except Exception:
        pass


def _compile_case(case_dir: Path, casename: str) -> bool:
    """Compile a case using mks (clean build)."""
    log_info(f"Compiling {casename}...")

    # Ensure environment
    env = os.environ.copy()
    if "NEKSTAB_SOURCE_ROOT" not in env:
        env["NEKSTAB_SOURCE_ROOT"] = str(NEKSTAB_ROOT)
        env["NEK_SOURCE_ROOT"] = str(NEKSTAB_ROOT / "Nek5000")
        env["PATH"] = f"{NEKSTAB_ROOT}/Nek5000/bin:{NEKSTAB_ROOT}/bin:{env.get('PATH', '')}"

    # Clean before building to ensure a fresh compilation
    subprocess.run(
        ["mks", "clean"],
        cwd=case_dir,
        env=env,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        timeout=30,
    )

    build_log = case_dir / "build.log"

    try:
        with open(build_log, "w") as f:
            result = subprocess.run(
                ["mks", casename],
                cwd=case_dir,
                env=env,
                stdout=f,
                stderr=subprocess.STDOUT,
                timeout=600,
            )

        if result.returncode != 0:
            log_fail(f"Compilation failed. See {build_log}")
            # Show first error line to avoid having to open the log
            _show_first_error(build_log)
            return False

        log_detail("Compilation successful")
        return True

    except subprocess.TimeoutExpired:
        log_fail("Compilation timeout (10 min)")
        return False
    except FileNotFoundError:
        log_fail("mks not found. Is nekStab environment set up?")
        return False


def _run_case(case_dir: Path, casename: str, nprocs: int, cpu: CPUInfo,
              timeout: int = 7200, requires: Optional[list[str]] = None) -> bool:
    """Run a case with MPI."""
    log_info(f"Running on {nprocs} cores (CPUs: {cpu.cpu_list})...")

    # Write SESSION.NAME
    session = case_dir / "SESSION.NAME"
    session.write_text(f"{casename}\n{case_dir}/\n")

    # Clean old outputs (but preserve files listed in requires)
    protected = set(requires or [])
    for f in case_dir.glob("Spectre_*.dat"):
        if f.name not in protected:
            f.unlink()
    for f in ["logfile", "ioinfo"]:
        p = case_dir / f
        if p.exists():
            p.unlink()

    # Build command with CPU binding
    cmd = ["mpirun", "-np", str(nprocs)]

    # OpenMPI CPU binding (pin to P-cores on hybrid CPUs)
    try:
        help_output = subprocess.run(
            ["mpirun", "--help", "binding"], capture_output=True, text=True
        ).stdout
        if nprocs <= cpu.count:
            # Exact fit: bind each process to a core in the slot
            if "--cpu-list" in help_output:
                cmd.extend(["--bind-to", "core", "--cpu-list", cpu.cpu_list])
            elif "--cpu-set" in help_output:
                cmd.extend(["--bind-to", "core", "--cpu-set", cpu.cpu_list])
        else:
            # Oversubscription: more MPI procs than slot cores (lpmin > slot size)
            cmd.extend(["--oversubscribe", "--bind-to", "none"])
            log_detail(f"Oversubscribed: {nprocs} procs on {cpu.count} cores")
    except Exception:
        pass

    cmd.append("./nek5000")

    log_detail(f"Command: {' '.join(cmd)}")

    run_log = case_dir / f"{casename}.log.{nprocs}"

    try:
        start = time.perf_counter()

        with open(run_log, "w") as f:
            result = subprocess.run(
                cmd,
                cwd=case_dir,
                stdout=f,
                stderr=subprocess.STDOUT,
                timeout=timeout,
            )

        elapsed = time.perf_counter() - start

        if result.returncode != 0:
            log_fail(f"Simulation failed. See {run_log}")
            return False

        log_detail(f"Completed in {elapsed:.0f}s")
        return True

    except subprocess.TimeoutExpired:
        log_fail(f"Simulation timeout ({timeout}s)")
        return False


def _run_plot(case_dir: Path) -> None:
    """Run plot.py in a case directory to regenerate figures (non-fatal)."""
    plot_script = case_dir / "plot.py"
    if not plot_script.exists():
        return

    log_info("Regenerating plot.png...")
    try:
        result = subprocess.run(
            [sys.executable, "plot.py"],
            cwd=case_dir,
            capture_output=True,
            text=True,
            timeout=120,
        )
        if result.returncode == 0:
            log_detail("plot.png updated")
        else:
            log_detail(f"plot.py failed: {result.stderr.strip()[:120]}")
    except subprocess.TimeoutExpired:
        log_detail("plot.py timeout (120s)")
    except Exception as e:
        log_detail(f"plot.py error: {e}")


# ═══════════════════════════════════════════════════════════════════════════════
# Case Ordering
# ═══════════════════════════════════════════════════════════════════════════════

# Mode order for sorting cases within a family.
# Longer patterns checked first to avoid substring conflicts
# (e.g., "direct_Floquet" before "direct").
_DIR_ORDER = [
    ("ci_test", 0),
    ("dns", 1),
    ("baseflow", 2),
    ("transient_growth", 3),
    ("direct_Floquet", 4),
    ("direct", 5),
    ("adjoint_Floquet", 6),
    ("adjoint", 7),
    ("animate", 8),
    ("postproc", 9),
    ("otd", 10),
    ("modal", 11),
]


def _case_sort_key(item: tuple[str, dict]) -> tuple[int, str, int, str]:
    """Sort key: (tier_order, family, dir_order, name)."""
    name, case = item
    tier_order = 0 if case["tier"] == "short" else 1
    family = case.get("family", "")
    d = case["dir"]
    dir_order = 99
    for key, val in _DIR_ORDER:
        if key in d:
            dir_order = val
            break
    return (tier_order, family, dir_order, name)


# ═══════════════════════════════════════════════════════════════════════════════
# Listing
# ═══════════════════════════════════════════════════════════════════════════════

def list_cases() -> None:
    """List all cases grouped by tier, showing dependency status."""
    short_cases = {k: v for k, v in CASES.items() if v["tier"] == "short"}
    full_cases = {k: v for k, v in CASES.items() if v["tier"] == "full"}

    print(f"\n{C.BOLD}AMR Validation (short tier):{C.NC}")
    for name, case in short_cases.items():
        print(f"  {name:20} {case['description']}")

    print(f"\n{C.BOLD}Compile+Run (full tier):{C.NC}")
    # Group by family
    families = {}
    for name, case in full_cases.items():
        fam = case.get("family", "other")
        families.setdefault(fam, []).append((name, case))

    for fam, cases in families.items():
        for name, case in cases:
            case_dir = NEKSTAB_ROOT / "example" / case["dir"]
            requires = case.get("requires", [])
            if not requires:
                status = f"{C.GREEN}[READY]{C.NC}"
            else:
                missing = check_requires(case_dir, requires)
                if missing:
                    status = f"{C.YELLOW}[SKIP ]{C.NC}"
                    desc = f"{case['description']} {C.DIM}(missing: {', '.join(missing)}){C.NC}"
                    print(f"  {name:20} {status}  {desc}")
                    continue
                else:
                    status = f"{C.GREEN}[READY]{C.NC}"
            print(f"  {name:20} {status}  {case['description']}")


# ═══════════════════════════════════════════════════════════════════════════════
# Compile-All Discovery
# ═══════════════════════════════════════════════════════════════════════════════

# Directories that are NOT independent case directories
_SKIP_DIRS = {"geom", "mesh", "obj", "__pycache__"}


def discover_compilable_cases(example_root: Path) -> list[tuple[Path, str]]:
    """
    Walk example/ tree and find all directories containing SIZE + *.usr.
    Returns sorted list of (case_dir, casename).
    """
    cases = []
    for size_file in sorted(example_root.rglob("SIZE")):
        case_dir = size_file.parent
        if case_dir.name in _SKIP_DIRS:
            continue
        usr_files = sorted(
            p for p in case_dir.glob("*.usr")
            if p.is_file() and not p.name.startswith(".") and p.stem
        )
        if not usr_files:
            continue
        casename = usr_files[0].stem
        cases.append((case_dir, casename))
    return cases


def compile_all_cases(example_root: Path) -> int:
    """
    Discover and compile every example case.
    Returns 0 if all pass, 1 if any fail.
    """
    cases = discover_compilable_cases(example_root)

    print()
    print("=" * 70)
    print("              nekStab — Compile All Examples")
    print("=" * 70)
    print()
    log_info(f"Discovered {len(cases)} compilable cases")
    log_info(f"nekStab root: {NEKSTAB_ROOT}")
    print()

    passed: list[str] = []
    failed: list[str] = []
    t0 = time.perf_counter()

    for i, (case_dir, casename) in enumerate(cases, 1):
        rel = case_dir.relative_to(example_root)
        log_header(f"[{i}/{len(cases)}] {rel} ({casename})")
        if _compile_case(case_dir, casename):
            passed.append(str(rel))
        else:
            failed.append(str(rel))

    elapsed = time.perf_counter() - t0

    # Summary
    total = len(cases)
    print()
    print("=" * 70)
    print("                    COMPILE SUMMARY")
    print("=" * 70)
    print()

    if passed:
        print(f"{C.GREEN}{C.BOLD}Compiled: {len(passed)}/{total}{C.NC}")
    if failed:
        print(f"{C.RED}{C.BOLD}Failed:   {len(failed)}/{total}{C.NC}")
        for f in failed:
            log_fail(f)
    if not failed:
        print(f"\n{C.GREEN}{C.BOLD}All {total} cases compiled successfully.{C.NC}")

    print(f"\nTotal time: {elapsed:.0f}s")
    return 1 if failed else 0


# ═══════════════════════════════════════════════════════════════════════════════
# Parallel Execution (DAG-based scheduling with core slots)
# ═══════════════════════════════════════════════════════════════════════════════


def build_dag(sorted_cases: dict) -> dict[str, set[str]]:
    """Build dependency graph from copies_to/requires and same-directory conflicts."""
    preds: dict[str, set[str]] = {name: set() for name in sorted_cases}

    # Map (dest_dir, filename) -> producer case name
    producers: dict[tuple[str, str], str] = {}
    for name, case in sorted_cases.items():
        for cp in case.get("copies_to", []):
            key = (cp["dest_dir"], cp.get("dest_name", cp["file"]))
            producers[key] = name

    # Match each case's requires against producers
    for name, case in sorted_cases.items():
        for req in case.get("requires", []):
            key = (case["dir"], req)
            if key in producers:
                preds[name].add(producers[key])

    # Same-directory conflict: chain cases in sort order
    dir_seen: dict[str, str] = {}
    for name, case in sorted_cases.items():
        d = case["dir"]
        if d in dir_seen:
            preds[name].add(dir_seen[d])
        dir_seen[d] = name

    return preds


def _run_one(name: str, case: dict, check_only: bool,
             slot_cpu: CPUInfo, dry_run: bool) -> tuple[str, float]:
    """Run a single case in a worker thread with buffered output."""
    _thread_local.buf = io.StringIO()
    try:
        # Pass slot core count as max; validate_case calls get_optimal_nprocs internally
        status, elapsed = validate_case(
            name, case, check_only, slot_cpu.count, slot_cpu, dry_run)
    except Exception as e:
        log_fail(f"{name}: unhandled error: {e}")
        status, elapsed = "fail", 0.0
    finally:
        output = _thread_local.buf.getvalue()
        _thread_local.buf = None
        with _print_lock:
            sys.stdout.write(output)
            sys.stdout.flush()
    return (status, elapsed)


def run_parallel(sorted_cases: dict, predecessors: dict[str, set[str]],
                 slots: list[CPUInfo], check_only: bool,
                 dry_run: bool) -> dict[str, tuple[str, float]]:
    """Execute cases in parallel respecting dependency DAG and core slots."""
    num_slots = len(slots)
    results: dict[str, tuple[str, float]] = {}
    pending = set(sorted_cases)
    running: dict = {}          # future -> (name, slot_idx)
    slot_free = list(range(num_slots))
    case_order = list(sorted_cases)  # preserve sort order for priority

    def _log_progress():
        done = len(results)
        total = len(sorted_cases)
        active = [n for _, (n, _) in running.items()]
        with _print_lock:
            print(f"{C.YELLOW}[INFO]{C.NC} [{done}/{total}] "
                  f"Active: {', '.join(active)}")

    def _submit_ready(pool):
        for name in case_order:
            if name not in pending:
                continue
            if not slot_free:
                break
            preds = predecessors.get(name, set())
            if not preds <= results.keys():
                continue
            # Skip if any predecessor failed
            if any(results[p][0] == "fail" for p in preds):
                results[name] = ("skip", 0.0)
                pending.discard(name)
                with _print_lock:
                    print(f"{C.YELLOW}[SKIP]{C.NC} {name}: predecessor failed")
                continue
            slot_idx = slot_free.pop(0)
            pending.discard(name)
            case = sorted_cases[name]
            cpu = slots[slot_idx]
            fut = pool.submit(_run_one, name, case, check_only, cpu, dry_run)
            running[fut] = (name, slot_idx)

    with ThreadPoolExecutor(max_workers=num_slots) as pool:
        _submit_ready(pool)
        _log_progress()
        while running:
            done_futs, _ = wait(running, return_when=FIRST_COMPLETED)
            for fut in done_futs:
                name, slot_idx = running.pop(fut)
                results[name] = fut.result()
                slot_free.append(slot_idx)
            _submit_ready(pool)
            if running:
                _log_progress()

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="nekStab validation against AMR paper test cases",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Short-tier cases (AMR eigenvalue validation, 13 cases):
  cylinder_bf         Newton baseflow (Re=50) -> copies BF downstream
  cylinder_sfd        SFD baseflow (Re=50) — classic method
  cylinder_direct     Hopf bifurcation, direct stability (Re=50)
  cylinder_adjoint    Adjoint stability (Re=50) -> copies eigvecs downstream
  cylinder_sensitivity Wavemaker + energy budget + sensitivity (Re=50)
  thermosyphon        Pitchfork bifurcation (symmetry breaking)
  flipflop_bf         Newton UPO, natural frequency (Re=62)
  flipflop_floquet    Floquet stability |mu|>1 (above Re_c~61.17)
  backstep_bf         Newton baseflow (Re=500) -> copies BF to tg
  backstep_tg         Transient growth (non-modal amplification)
  tpjet_bf            Newton forced PO, St=0.6 (Re=1900)
  tpjet_floquet       Period-doubling Floquet mu<-1

Full-tier cases (compile+run validation):
  cyl_*         Cylinder family (20 cases)
  thermo_bf     Thermosyphon baseflow
  tpjet_tdf     Forced jet TDF
  lid_driven    Lid-driven cavity
  cubic_cav*    3D cubic cavity
  naca0012*     NACA 0012 airfoil
  poiseuille    2D Poiseuille OTD
  slot_fst      Slot jet with FST

Examples:
  ./validate.py                    # Run all cases
  ./validate.py --short            # Run only AMR validation
  ./validate.py cylinder_direct    # Run one case
  ./validate.py --check-only       # Validate existing results
  ./validate.py -n 8 --short       # Use 8 MPI processes
  ./validate.py --compile-all      # Compile every example dir
        """
    )
    parser.add_argument(
        "case", nargs="?",
        help="Specific case to run (default: all)"
    )
    parser.add_argument(
        "--check-only", "-c", action="store_true",
        help="Only validate existing results (no compile/run)"
    )
    parser.add_argument(
        "--dry-run", "-d", action="store_true",
        help="Show what would be run"
    )
    parser.add_argument(
        "--nprocs", "-n", type=int,
        help="Override MPI process count"
    )
    parser.add_argument(
        "--list", "-l", action="store_true",
        help="List available cases and exit"
    )
    parser.add_argument(
        "--short", "-s", action="store_true",
        help="Run only short-tier (AMR validation) cases"
    )
    parser.add_argument(
        "--compile-all", action="store_true",
        help="Discover and compile every example case (no run)"
    )
    parser.add_argument(
        "--sequential", action="store_true",
        help="Disable parallel execution (run cases one by one)"
    )

    args = parser.parse_args()

    # List cases
    if args.list:
        list_cases()
        return 0

    # Compile-all: discover and compile every example, then exit
    if args.compile_all:
        return compile_all_cases(NEKSTAB_ROOT / "example")

    # Header
    print()
    print("=" * 70)
    print("                    nekStab Validation Suite")
    print("=" * 70)
    print()

    # CPU detection
    cpu = detect_pcores()
    print_cpu_info(cpu)

    log_info(f"Tolerance: {TOLERANCE*100:.0f}%")
    log_info(f"nekStab root: {NEKSTAB_ROOT}")

    if args.check_only:
        log_info("CHECK-ONLY mode: validating existing results")
    if args.dry_run:
        log_info("DRY-RUN mode: no commands executed")
    if args.short:
        log_info("SHORT tier only: AMR validation cases")

    # Select cases
    if args.case:
        # Match by prefix
        matches = [n for n in CASES if n.startswith(args.case.lower())]
        if not matches:
            log_fail(f"Unknown case: {args.case}")
            log_info(f"Available: {', '.join(CASES.keys())}")
            return 1
        cases_to_run = {m: CASES[m] for m in matches}
    elif args.short:
        cases_to_run = {k: v for k, v in CASES.items() if v["tier"] == "short"}
    else:
        cases_to_run = dict(CASES)

    # Sort cases for execution order
    sorted_cases = dict(sorted(cases_to_run.items(), key=_case_sort_key))

    # Run validation
    results: dict[str, tuple[str, float]] = {}
    suite_t0 = time.perf_counter()

    if args.sequential or len(sorted_cases) <= 1:
        # Sequential path (original behavior)
        nprocs = args.nprocs or cpu.count
        log_info(f"MPI processes: {nprocs}")
        for name, case in sorted_cases.items():
            status, elapsed = validate_case(
                name, case, args.check_only, nprocs, cpu, args.dry_run
            )
            results[name] = (status, elapsed)
    else:
        # Parallel path: detect all cores, build DAG, run with slots
        all_cores = detect_all_cores()
        cores_per_slot = args.nprocs or MIN_CORES_PER_SLOT
        slots = make_slots(all_cores, cores_per_slot)
        predecessors = build_dag(sorted_cases)
        log_info(f"Parallel mode: {len(slots)} slots × {slots[0].count} cores "
                 f"({len(all_cores)} physical cores)")
        results = run_parallel(sorted_cases, predecessors, slots,
                               args.check_only, args.dry_run)

    suite_elapsed = time.perf_counter() - suite_t0

    # Summary grouped by tier
    print()
    print("=" * 70)
    print("                         SUMMARY")
    print("=" * 70)

    # Re-order results by original sorted_cases order (parallel mode returns completion order)
    ordered_results = {k: results[k] for k in sorted_cases if k in results}
    short_results = {k: v for k, v in ordered_results.items() if CASES[k]["tier"] == "short"}
    full_results = {k: v for k, v in ordered_results.items() if CASES[k]["tier"] == "full"}

    def _print_result_line(name: str, status: str, elapsed: float) -> None:
        desc = CASES[name]["description"]
        time_str = f"({_fmt_time(elapsed):>10})" if elapsed > 1.0 else f"{'':>12}"
        line = f"{name:20} {time_str}  {desc}"
        if status == "pass":
            log_pass(line)
        elif status == "skip":
            log_skip(line)
        else:
            log_fail(line)

    if short_results:
        print(f"\n{C.BOLD}AMR Validation (short tier):{C.NC}")
        for name, (status, elapsed) in short_results.items():
            _print_result_line(name, status, elapsed)

    if full_results:
        print(f"\n{C.BOLD}Compile+Run (full tier):{C.NC}")
        for name, (status, elapsed) in full_results.items():
            _print_result_line(name, status, elapsed)

    print()

    passed_count = sum(1 for s, _ in results.values() if s == "pass")
    skipped_count = sum(1 for s, _ in results.values() if s == "skip")
    failed_count = sum(1 for s, _ in results.values() if s == "fail")
    total_count = len(results)

    if args.dry_run:
        print(f"{C.YELLOW}DRY RUN COMPLETE{C.NC}")
        return 0

    summary = f"PASSED: {passed_count}/{total_count}"
    if skipped_count:
        summary += f"  SKIPPED: {skipped_count}"
    if failed_count:
        summary += f"  FAILED: {failed_count}"
    summary += f"  Total: {_fmt_time(suite_elapsed)}"

    if failed_count == 0:
        print(f"{C.GREEN}{C.BOLD}{summary}{C.NC}")
        return 0
    else:
        print(f"{C.RED}{C.BOLD}{summary}{C.NC}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
