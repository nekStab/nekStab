"""
Chebyshev collocation differentiation matrices.

Port of example/slot_FST/preprocessFST/chebCL.m (MATLAB).
Generates spectral differentiation matrices D1, D2, D3, D4 on a
Chebyshev-Gauss-Lobatto grid with optional domain stretching.

Reference:
  chebCL.m (original MATLAB implementation)
  Spectral methods for boundary layer applications.
"""

import numpy as np


def cheb_collocation(N, X1, Xint):
    """
    Chebyshev collocation differentiation matrices on a stretched grid.

    Generates Chebyshev-Gauss-Lobatto collocation nodes on the interval [0, X1]
    with optional stretching controlled by Xint. Computes spectral differentiation
    matrices (1st, 2nd, 3rd, 4th order) via the mapped Chebyshev formulas.

    Parameters
    ----------
    N : int
        Number of collocation points (actual grid has N+1 nodes, MATLAB chebCL decrements N).
    X1 : float
        Right endpoint of physical domain [0, X1].
    Xint : float
        Interior domain stretching parameter. Affects the mapping from standard
        Chebyshev interval [-1, 1] to [0, X1].

    Returns
    -------
    dict
        Dictionary with keys:
        - 'D1': (N+1, N+1) ndarray, first derivative matrix
        - 'D2': (N+1, N+1) ndarray, second derivative matrix (D1 @ D1)
        - 'D3': (N+1, N+1) ndarray, third derivative matrix (D1 @ D2)
        - 'D4': (N+1, N+1) ndarray, fourth derivative matrix (D2 @ D2)
        - 'y': (N+1,) ndarray, physical collocation nodes in [0, X1]

    Notes
    -----
    Faithful port of MATLAB chebCL.m. The algorithm:
    1. Generates Chebyshev-Gauss-Lobatto nodes on [-1, 1]
    2. Maps to [0, X1] via affine transformation with stretching
    3. Computes raw Chebyshev differentiation matrix (spectral endpoint rule)
    4. Scales by Jacobian of coordinate mapping
    5. Composes to form higher-order derivatives

    The spectral endpoint rule uses c-weights c(1)=c(N+1)=2, c(2..N)=1.
    """
    # Decrement N to match MATLAB convention (MATLAB does N=N-1 first)
    N = N - 1

    # Generate Chebyshev-Gauss-Lobatto nodes on [-1, 1]
    # ksi(j) = cos(pi*(j-1)/N) for j = 1, 2, ..., N+1
    ksi = np.array([np.cos(np.pi * j / N) for j in range(N + 1)])

    # Map to physical domain [0, X1] with stretching
    a = Xint * X1 / (X1 - 2 * Xint)
    b = 1.0 + 2.0 * a / X1
    y = a * (1.0 + ksi) / (b - ksi)

    # c-weights for spectral endpoint rule
    c = np.ones(N + 1)
    c[0] = 2.0
    c[N] = 2.0

    # Compute raw Chebyshev differentiation matrix d (in the reference [-1,1]
    # coordinate ksi).  This is the classic collocation-derivative matrix.
    d = np.zeros((N + 1, N + 1))
    for i in range(N + 1):
        for j in range(N + 1):
            if i != j:
                # Off-diagonal: (c(i)/c(j)) * (-1)^(i+j) / (ksi(i) - ksi(j))
                d[i, j] = (c[i] / c[j]) * ((-1) ** (i + j)) / (ksi[i] - ksi[j])
            elif 0 < i < N:
                # Interior diagonal: -ksi(i) / (2*(1 - ksi(i)^2)).  This formula
                # is SINGULAR at the endpoints (ksi = +-1), so we skip i=0 and
                # i=N here and assign their exact spectral values just below.
                d[i, i] = -ksi[i] / (2.0 * (1.0 - ksi[i] ** 2))

    # Endpoint diagonal entries (the closed-form spectral corner values).
    d[0, 0] = (2.0 * N**2 + 1.0) / 6.0
    d[N, N] = -(2.0 * N**2 + 1.0) / 6.0

    # Scale by Jacobian of coordinate mapping
    # S = a * (b + 1) / (y + a)^2 (diagonal matrix)
    S = a * (b + 1.0) / (y + a) ** 2

    # Compute physical derivatives
    D1 = np.diag(S) @ d
    D2 = D1 @ D1
    D3 = D1 @ D2
    D4 = D2 @ D2

    return {
        "D1": D1,
        "D2": D2,
        "D3": D3,
        "D4": D4,
        "y": y,
    }


if __name__ == "__main__":
    # Self-test: verify differentiation on a known function
    print("=" * 70)
    print("Chebyshev Collocation Differentiation Self-Test")
    print("=" * 70)

    # Build grid with moderate N
    N = 48
    X1 = 0.3
    Xint = 0.1

    result = cheb_collocation(N, X1, Xint)
    D1 = result["D1"]
    D2 = result["D2"]
    y = result["y"]

    # Test function: f(y) = sin(k*y) with k = pi
    k = np.pi
    f = np.sin(k * y)

    # Exact derivatives
    df_exact = k * np.cos(k * y)
    d2f_exact = -k**2 * np.sin(k * y)

    # Computed derivatives via spectral matrices
    df_computed = D1 @ f
    d2f_computed = D2 @ f

    # Errors at interior nodes (exclude endpoints which have boundary effects)
    interior = slice(1, -1)
    err_d1 = np.max(np.abs(df_computed[interior] - df_exact[interior]))
    err_d2 = np.max(np.abs(d2f_computed[interior] - d2f_exact[interior]))

    print(f"\nGrid: N+1={N+1} nodes, X1={X1}, Xint={Xint}")
    print(f"Test function: f(y) = sin({k:.4f}*y)")
    print(f"\nFirst derivative (D1):")
    print(f"  Max error (interior nodes): {err_d1:.4e}")
    print(f"Second derivative (D2):")
    print(f"  Max error (interior nodes): {err_d2:.4e}")

    # Check success
    tol = 1e-6
    if err_d1 < tol and err_d2 < tol:
        print(f"\nPASS: Both errors < {tol:.0e}")
    else:
        print(f"\nWARNING: Some errors >= {tol:.0e}")
        if err_d1 >= tol:
            print(f"  D1 error {err_d1:.4e} exceeds tolerance")
        if err_d2 >= tol:
            print(f"  D2 error {err_d2:.4e} exceeds tolerance")

    print("\nMatrix shapes:")
    print(f"  D1: {D1.shape}, D2: {D2.shape}, D3: {result['D3'].shape}, D4: {result['D4'].shape}")
    print(f"  y:  {y.shape}")
    print("=" * 70)
