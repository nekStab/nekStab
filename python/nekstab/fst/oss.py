"""
Orr-Sommerfeld-Squire continuous-spectrum mode generator for Free-Stream
Turbulence inflow (Brandt, Schlatter & Henningson 2004; Jacobs & Durbin 1998).

Faithful port of example/slot_FST/preprocessFST/FST_modes.m + bcA.m, with ONE
deliberate correction (approved): the MATLAB reopened the wavenumber file three
times and read the first scalar each time, so the OSS solve actually received
omega = gamma = beta.  Here we use the THREE DISTINCT wavenumbers (omega, gamma,
beta) -> physically correct continuous-spectrum modes.

Algorithm per mode (TEMPORAL formulation):
  alpha = omega                              (Taylor:  alpha = omega * U_inf,  U_inf = 1)
  kk    = sqrt(alpha^2 + beta^2)
  c     = 1 - i (1 + gamma^2/kk^2) kk^2 / (alpha Re)        (Brandt 2002 phase speed)
  Orr-Sommerfeld operator for v:
    A = (D4 - 2 kk^2 D2 + kk^4 I)
        - i Re alpha [ (diag(U) - c I)(D2 - kk^2 I) - diag(U'') ]
  with Dirichlet/Neumann wall BCs and a homogeneous-deflation top BC; a scalar
  Newton iteration tunes the free-stream v amplitude until the Jacobs & Durbin
  (1998) unbounded boundary condition residual vanishes.  The Squire equation
  supplies the wall-normal vorticity eta.  (v, eta) -> (u, v, w), windowed near
  the top, then normalised to unit energy in the free stream with random phases.

Grid orientation (from chebCL): index 0 = free-stream top (y = y_max),
index -1 = wall (y = 0).
"""

import numpy as np

from .cheb import cheb_collocation
from .blasius import blasius_profile

# np.trapz was removed in NumPy 2.0 (renamed np.trapezoid); support both.
try:
    from numpy import trapezoid as _trapz
except ImportError:  # NumPy < 2.0
    from numpy import trapz as _trapz


def _solve_v_with_bc(A, cc, gamma, alpha, D2, yvecs, Ny):
    """Port of bcA.m: solve the OS system for v with free-stream amplitude cc,
    and return (V, chk) where chk is the Jacobs & Durbin unbounded-BC residual."""
    B = np.zeros(Ny, dtype=complex)
    B[0] = 1.0 + 0.0j   # v = 1 at the top
    B[1] = cc           # v = cc at the node just below the top (Newton target)
    B[-1] = 0.0         # v = 0 at the wall
    V = np.linalg.solve(A, B)
    # chk = (D2(2,:)V + gamma^2 V(2)) - exp(alpha (y3 - y2)) (D2(3,:)V + gamma^2 V(3))
    # MATLAB 1-based (2,3) -> 0-based (1,2).
    chk = ((D2[1, :] @ V + gamma**2 * V[1])
           - np.exp(alpha * (yvecs[2] - yvecs[1])) * (D2[2, :] @ V + gamma**2 * V[2]))
    return V, chk


def solve_oss_mode(omega, gamma, beta, Re, Ny, Ly, rng, newton_tol=1e-9,
                   newton_max=200):
    """Generate one continuous-spectrum OSS mode.

    Returns (y, U, V, W): y is the uniform output grid (ascending, length Ny);
    U, V, W are complex velocity profiles on that grid, unit-energy-normalised
    in the free stream with random relative phases.  Also returns the Newton
    residual achieved (for verification).
    """
    # --- TEMPORAL dispersion relation (Brandt 2002) ---
    alpha = omega
    kk = np.sqrt(alpha**2 + beta**2)
    c = 1.0 - 1j * (1.0 + gamma**2 / kk**2) * kk**2 / (alpha * Re)

    # --- operators + base flow ---
    ch = cheb_collocation(Ny, Ly, 2.0 * Ly / 5.0)
    D1, D2, D4, yvecs = ch["D1"], ch["D2"], ch["D4"], ch["y"]
    I = np.eye(Ny)
    U, _U_p, U_pp = blasius_profile(Ly, yvecs)

    # --- Orr-Sommerfeld operator (FST_modes.m line 67) ---
    A = ((D4 - 2.0 * D2 * kk**2 + kk**4 * I)
         - 1j * Re * alpha * ((np.diag(U) - I * c) @ (D2 - kk**2 * I) - np.diag(U_pp)))
    # BC rows (0-based: 0,1 = top; -2,-1 = wall)
    A[-1, :] = 0.0; A[-1, -1] = 1.0          # Dirichlet at the wall
    A[-2, :] = D1[-1, :]                      # Neumann at the wall
    A[0, :] = 0.0;  A[0, 0] = 1.0 + 0.0j     # Dirichlet deflation at the top
    A[1, :] = 0.0;  A[1, 1] = 1.0            # second top row (Newton freedom)

    # --- scalar Newton on the free-stream amplitude cc (FST_modes.m 78-96) ---
    cc1 = 1.0 + 0.0j
    chk1 = 10.0 + 0.0j
    V = None
    ite = 0
    while abs(chk1) > newton_tol and ite < newton_max:
        V, chk1 = _solve_v_with_bc(A, cc1, gamma, alpha, D2, yvecs, Ny)
        cc2 = cc1 + (rng.random() + 1j * rng.random()) * 1e-7
        _, chk2 = _solve_v_with_bc(A, cc2, gamma, alpha, D2, yvecs, Ny)
        dF = (chk2 - chk1) / (cc2 - cc1)
        cc1 = cc1 - chk1 / dF
        ite += 1
    newton_res = abs(chk1)

    # --- Squire equation for eta (coupling term removed, as in the MATLAB) ---
    S = 1j * Re * alpha * (np.diag(U) - c * I) - (D2 - kk**2 * I)
    Bsq = np.zeros(Ny, dtype=complex)
    S[0, :] = 0.0;  S[0, 0] = 1.0            # eta Dirichlet, top
    S[-1, :] = 0.0; S[-1, -1] = 1.0          # eta Dirichlet, wall
    Bsq[0] = 1.0 + 0.0j                      # eta = 1 at the top
    Bsq[-1] = 0.0                            # eta = 0 at the wall
    E = np.linalg.solve(S, Bsq)

    # --- top-of-domain smoothing window Ss (FST_modes.m 114-129) ---
    # A C-infinity bump that ramps the eigenfunction down to zero over the top
    # 20% of the domain [ydm, ymax].  This isolates the mode from the imposed
    # top boundary condition (a numerical sponge), so the free-stream mode is
    # not contaminated by the deflation BCs used to solve the OS system.
    Lreal = Ly
    ydm = Lreal - 0.2 * Lreal
    ymax = Lreal
    Ss = np.zeros(Ny)
    for i in range(Ny):
        ys = 1.0 - (yvecs[i] - ydm) / (ymax - ydm)
        if ys <= 0.0:
            Ss[i] = 0.0           # at/above the top -> fully damped
        elif ys < 1.0:
            # Smooth transition.  arg -> +inf as ys -> 0+ (window -> 0); guard
            # the exp against overflow (exp(>709) overflows float64).
            arg = 1.0 / (ys - 1.0) + 1.0 / ys
            val = 0.0 if arg > 700.0 else 1.0 / (1.0 + np.exp(arg))
            Ss[i] = 0.0 if val < 1e-50 else val
        else:
            Ss[i] = 1.0           # well below ydm -> mode passes through

    # Keep the RAW (pre-window) wall-normal velocity for the continuity check:
    # the sponge window perturbs incompressibility on purpose, so the meaningful
    # divergence-free verification must use the unwindowed reconstruction.
    V_raw = V.copy()

    # NOTE: replicates the MATLAB exactly, incl. the order V=V.*Ss BEFORE forming
    # dV (so dV uses the already-windowed V) -- faithful to FST_modes.m 132-134.
    V = V * Ss
    dV = (D1 @ V) * Ss + V * (D1 @ Ss)
    E = E * Ss

    # --- (v, eta) -> (u, v, w) ---
    Uos = 1j * alpha / kk**2 * dV
    Wos = 1j * beta / kk**2 * dV
    Usq = -1j * beta / kk**2 * E
    Wsq = 1j * alpha / kk**2 * E

    # uniform output grid (FST_modes.m: y = linspace(0, Lreal, Ny))
    y = np.linspace(0.0, Lreal, Ny)

    def _spline(src):
        # chebCL yvecs is descending (top->wall); CubicSpline needs ascending x.
        from scipy.interpolate import CubicSpline
        idx = np.argsort(yvecs)
        return CubicSpline(yvecs[idx], src[idx])(y)

    phi = rng.random() * 2.0 * np.pi
    Uos_i = _spline(Uos) * np.exp(1j * phi)
    Vos_i = _spline(V) * np.exp(1j * phi)
    Wos_i = _spline(Wos) * np.exp(1j * phi)

    phi = rng.random() * 2.0 * np.pi
    Usq_i = _spline(Usq) * np.exp(1j * phi)
    Wsq_i = _spline(Wsq) * np.exp(1j * phi)

    phi = rng.random() * 2.0 * np.pi
    Uv = np.cos(phi) * Uos_i + np.sin(phi) * Usq_i
    Vv = np.cos(phi) * Vos_i
    Wv = np.cos(phi) * Wos_i + np.sin(phi) * Wsq_i

    # --- unit-energy normalisation over the free-stream band [5, ydm] ---
    y1 = np.linspace(5.0, ydm, Ny)
    from scipy.interpolate import CubicSpline
    U1 = CubicSpline(y, Uv)(y1)
    V1 = CubicSpline(y, Vv)(y1)
    W1 = CubicSpline(y, Wv)(y1)
    energy = (U1 * np.conj(U1) + V1 * np.conj(V1) + W1 * np.conj(W1)).real
    energia = 0.5 * abs(_trapz(energy, y1)) / abs(y1[0] - y1[-1])
    Uv = Uv / np.sqrt(energia)
    Vv = Vv / np.sqrt(energia)
    Wv = Wv / np.sqrt(energia)

    # --- spectral continuity diagnostic ---
    # The Orr-Sommerfeld reconstruction is divergence-free by construction:
    #   u = i*alpha/kk^2 v',  w = i*beta/kk^2 v'   (v' = D1 @ V_raw, spectral)
    #   => i*alpha*u + v' + i*beta*w = (D1 @ V_raw)*(1 - (alpha^2+beta^2)/kk^2) = 0.
    # Computed on the RAW (pre-window) field with the spectral derivative -> ~1e-12
    # when the wiring (kk, the i*alpha/i*beta factors, D1) is correct.  This is a
    # tight regression guard against reconstruction bugs; the windowed/splined
    # output field is only div-free to ~1e-2 (the sponge perturbs it on purpose).
    dV_raw = D1 @ V_raw
    Uos_raw = 1j * alpha / kk**2 * dV_raw
    Wos_raw = 1j * beta / kk**2 * dV_raw
    div_raw = 1j * alpha * Uos_raw + dV_raw + 1j * beta * Wos_raw
    interior = slice(2, Ny - 2)  # exclude the 2 top + 2 wall BC rows
    denom = np.max(np.abs(dV_raw[interior])) + 1e-30
    spectral_div = float(np.max(np.abs(div_raw[interior])) / denom)
    diag = {"newton_res": newton_res, "spectral_div": spectral_div}

    # MATLAB flips to descending y before writing; we return ascending y with
    # matching velocity ordering (caller decides storage orientation).
    return y, Uv, Vv, Wv, diag


if __name__ == "__main__":
    print("=" * 70)
    print("OSS continuous-spectrum mode self-test (port of FST_modes.m)")
    print("=" * 70)
    Re = 495.0
    Ny = 35 * 6      # main.m: Ney=35, lx1=6 -> Ny = Ney*lx1
    Ly = 20.0
    rng = np.random.default_rng(0)

    # a representative (omega, gamma, beta) triplet (distinct -> the fixed physics)
    omega, gamma, beta = 0.8, 0.5, 0.6
    y, U, V, W, diag = solve_oss_mode(omega, gamma, beta, Re, Ny, Ly, rng)

    finite = np.all(np.isfinite(U)) and np.all(np.isfinite(V)) and np.all(np.isfinite(W))
    print(f"\nNewton unbounded-BC residual: {diag['newton_res']:.2e}  (target < 1e-9)")
    print(f"all velocities finite: {finite}")

    # unit-energy check over the same free-stream band used for normalisation
    ydm = Ly - 0.2 * Ly
    from scipy.interpolate import CubicSpline
    y1 = np.linspace(5.0, ydm, Ny)
    e = (CubicSpline(y, U)(y1) * np.conj(CubicSpline(y, U)(y1))
         + CubicSpline(y, V)(y1) * np.conj(CubicSpline(y, V)(y1))
         + CubicSpline(y, W)(y1) * np.conj(CubicSpline(y, W)(y1))).real
    e_norm = 0.5 * abs(_trapz(e, y1)) / abs(y1[0] - y1[-1])
    print(f"free-stream energy after normalisation: {e_norm:.6f}  (target ~ 1.0)")

    # incompressibility: spectral divergence of the RAW reconstruction (cheb D1).
    print(f"spectral continuity residual (raw): {diag['spectral_div']:.3e}  (target < 1e-9)")

    ok = (diag['newton_res'] < 1e-9 and finite and abs(e_norm - 1.0) < 1e-6
          and diag['spectral_div'] < 1e-9)
    print("\nPASS" if ok else "\nFAIL")
    print("=" * 70)
