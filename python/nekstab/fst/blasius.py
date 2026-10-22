"""
Blasius boundary-layer similarity profile, dimensionless by displacement thickness.

Faithful port of:
  example/slot_FST/preprocessFST/blasius_profil.m
  example/slot_FST/preprocessFST/blasiusa.m

IMPORTANT — what the MATLAB actually runs:
The shooting/bisection loops in blasius_profil.m are COMMENTED OUT.  The code
that executes uses the pre-converged constants

    W = 1.480620611513586          (ODE parameter)
    c = f''(0) = 0.571371930837631 (wall shear, displacement-thickness-normalized)

integrates the Blasius system ONCE on [0, H] with f0 = [0, 0, c], then
spline-interpolates the result onto the requested nodes.  We reproduce exactly
that executed path (no shooting).

ODE (blasiusa.m):  f1' = f2,  f2' = f3,  f3' = -W*f1*f3
with f1(0)=0, f2(0)=0, f3(0)=c, where f1=f, f2=f'=U, f3=f''.

Outputs match blasius_profil.m column ordering:
    Blas    = f(:,2) = f'   = U      (streamwise velocity, U_inf = 1)
    Blas_p  = f(:,3) = f''           (velocity gradient; Blas_p(wall) = c)
    Blas_pp = f(:,4) = -W*f*f'' = f''' = U''  (used as U'' in the OSS operator)
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

# Pre-converged constants from blasius_profil.m (the active, non-commented values).
W_BLASIUS = 1.480620611513586
C_WALL = 0.571371930837631  # f''(0)


def _blasiusa(y, f):
    """Blasius ODE RHS (port of blasiusa.m): df = [f2, f3, -W*f1*f3]."""
    return [f[1], f[2], -W_BLASIUS * f[0] * f[2]]


def blasius_profile(H, yvecs, n_int=4000):
    """Return (Blas, Blas_p, Blas_pp) on yvecs, matching blasius_profil.m.

    Parameters
    ----------
    H : float
        Integration height (the MATLAB passes the OSS domain height Ly here).
    yvecs : array_like
        Physical nodes to interpolate the profile onto (e.g. the Chebyshev grid).
    n_int : int
        Number of dense integration samples on [0, H] used to build the spline
        (the MATLAB uses ode45 adaptive points + spline interp1; we use a dense
        grid + CubicSpline, with a tight tolerance so the base profile is well
        resolved for the downstream eigenvalue problem).

    Returns
    -------
    (Blas, Blas_p, Blas_pp) : tuple of ndarray, each shape == np.shape(yvecs)
        Blas    = U      (velocity, -> 1 in the free stream)
        Blas_p  = U'     (= f''; equals C_WALL at the wall)
        Blas_pp = U''    (= f''' = -W*f*f''; enters the Orr-Sommerfeld operator)
    """
    yvecs = np.asarray(yvecs, dtype=float)
    f0 = [0.0, 0.0, C_WALL]
    t_eval = np.linspace(0.0, H, n_int)
    sol = solve_ivp(_blasiusa, [0.0, H], f0, t_eval=t_eval,
                    rtol=1e-12, atol=1e-12, method='RK45')
    y = sol.t
    f1, f2, f3 = sol.y[0], sol.y[1], sol.y[2]
    f4 = -W_BLASIUS * f1 * f3  # = f''' = U''

    blas = CubicSpline(y, f2)(yvecs)     # U
    blas_p = CubicSpline(y, f3)(yvecs)   # f''
    blas_pp = CubicSpline(y, f4)(yvecs)  # U''
    return blas, blas_p, blas_pp


if __name__ == "__main__":
    print("=" * 70)
    print("Blasius profile self-test (faithful port of blasius_profil.m)")
    print("=" * 70)

    # Integrate over a tall enough domain for the BL to reach the free stream.
    H = 20.0
    yv = np.linspace(0.0, H, 200)
    U, Up, Upp = blasius_profile(H, yv)

    # f''(0) MUST equal the hardcoded C_WALL = 0.571371930837631 (this is the
    # MATLAB's actual wall shear, NOT the eta-normalized Howarth value 0.33206).
    f2pp_wall = Up[0]
    print(f"\nf''(0) [should be {C_WALL:.12f}]: {f2pp_wall:.12f}")
    print(f"  abs error vs hardcoded c: {abs(f2pp_wall - C_WALL):.3e}")
    print(f"U_inf = U(y_max) [should -> 1]: {U[-1]:.10f}")
    print(f"U monotonic increasing: {bool(np.all(np.diff(U) >= -1e-12))}")
    print(f"U'' range (Orr-Sommerfeld input): [{Upp.min():.6f}, {Upp.max():.6f}]")

    ok = (abs(f2pp_wall - C_WALL) < 1e-3 and abs(U[-1] - 1.0) < 1e-3
          and bool(np.all(np.diff(U) >= -1e-12)))
    print("\nPASS" if ok else "\nFAIL")
    print("=" * 70)
