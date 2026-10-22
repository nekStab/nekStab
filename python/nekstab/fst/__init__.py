"""Free-Stream-Turbulence (FST) inflow generation for nekStab.

Generates a synthetic turbulent inflow as a superposition of continuous-spectrum
Orr-Sommerfeld-Squire modes of the Blasius boundary layer (Brandt, Schlatter &
Henningson 2004; Jacobs & Durbin 1998).  This is a pure-Python, case-agnostic
re-implementation of the legacy MATLAB pipeline in
example/slot_FST/preprocessFST/, with the wavenumber-reading bug fixed (the
original collapsed omega = gamma = beta).

Public API
----------
generate_wavenumbers(numk, kkini, kkfin, seed)      -> (n_modes, 3) of (omega, gamma, beta)
solve_oss_mode(omega, gamma, beta, Re, Ny, Ly, rng) -> (y, U, V, W, newton_res)
blasius_profile(H, yvecs)                            -> (U, U', U'')  [base flow]
cheb_collocation(N, X1, Xint)                        -> {D1, D2, D3, D4, y}
generate_fst_inflow(...) / write_inflow_file(...)    -> high-level driver (see generator.py)

CLI
---
    python -m nekstab.fst --re 495 --ny 210 --ly 20 --numk 25 \
        --kmin <kmin> --kmax <kmax> --out <case>/FST_inflow.dat
"""

from .cheb import cheb_collocation
from .blasius import blasius_profile
from .wavenumber import generate_wavenumbers
from .oss import solve_oss_mode
from .generator import generate_fst_inflow, write_inflow_file

__all__ = [
    "cheb_collocation",
    "blasius_profile",
    "generate_wavenumbers",
    "solve_oss_mode",
    "generate_fst_inflow",
    "write_inflow_file",
]
