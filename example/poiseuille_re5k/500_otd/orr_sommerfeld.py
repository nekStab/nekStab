#!/usr/bin/env python3
"""orr_sommerfeld.py: least damped eigenvalues of plane Poiseuille flow, U = 1 - y**2, y in [-1, 1].

Re is based on the centreline velocity and the half-height. The equation is

    (D^2 - a^2)^2 phi = i a Re [(U - c) (D^2 - a^2) phi - U'' phi],    lambda = -i a c,

for the stream function phi(y) exp(i a x), with phi = phi' = 0 at y = +-1. It is solved by
Chebyshev collocation: the four boundary rows replace four rows of the generalised eigenproblem,
and the eigenvalues of huge modulus that this leaves are discarded.

Checks: at Re = 5772.22, a = 1.02056 (the critical point) the growth rate is -3e-9, and the
result at Re = 5000, a = 1 is the same to 1e-8 with N = 120 and N = 200 points.

USAGE: python orr_sommerfeld.py [Re [alpha [N]]]      (defaults 5000, 1, 120)
"""

import sys

import numpy as np
from scipy.linalg import eig


def cheb(n):
    """Chebyshev points x_j = cos(j pi / n) and the first-derivative matrix (Trefethen, 2000, ch. 6)."""
    x = np.cos(np.pi * np.arange(n + 1) / n)
    c = np.hstack([2.0, np.ones(n - 1), 2.0]) * (-1.0) ** np.arange(n + 1)
    big_x = np.tile(x, (n + 1, 1)).T
    dx = big_x - big_x.T
    d = np.outer(c, 1.0 / c) / (dx + np.eye(n + 1))
    d = d - np.diag(d.sum(axis=1))
    return d, x


def spectrum(re, alpha, n=120):
    """Eigenvalues lambda = sigma + i omega of the Orr-Sommerfeld problem, sorted by decreasing sigma."""
    d1, y = cheb(n)
    d2 = d1 @ d1
    d4 = d2 @ d2
    eye = np.eye(n + 1)
    u = 1.0 - y**2
    upp = -2.0 * np.ones(n + 1)
    lap = d2 - alpha**2 * eye
    # A phi = c B phi
    a_mat = (d4 - 2.0 * alpha**2 * d2 + alpha**4 * eye) / (1j * alpha * re) - (
        np.diag(u) @ lap - np.diag(upp)
    )
    b_mat = -lap
    a_mat = a_mat.astype(complex)
    b_mat = b_mat.astype(complex)
    for row, bc in ((0, eye[0]), (1, d1[0]), (n - 1, d1[n]), (n, eye[n])):
        a_mat[row, :] = bc
        b_mat[row, :] = 0.0
    c, _ = eig(a_mat, b_mat)
    lam = -1j * alpha * c[np.isfinite(c)]
    lam = lam[np.abs(lam) < 5.0]
    return lam[np.argsort(-lam.real)]


if __name__ == "__main__":
    re = float(sys.argv[1]) if len(sys.argv) > 1 else 5000.0
    alpha = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 120
    lam = spectrum(re, alpha, n)
    print(f"Re = {re}, alpha = {alpha}, N = {n}")
    for k in range(4):
        print(f"  lambda_{k + 1} = {lam[k].real:+.9f} {lam[k].imag:+.9f} i")
