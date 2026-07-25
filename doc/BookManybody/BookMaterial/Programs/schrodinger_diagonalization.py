"""
Schrodinger's equation solved by diagonalisation.

Companion code to chapter 1 of *Quantum mechanics for Many-particle Systems*.

Instead of integrating the differential equation we discretise it on a grid
and diagonalise the resulting matrix.  With the three-point formula for the
second derivative,

    u''(x_i) = ( u_{i+1} - 2 u_i + u_{i-1} ) / h^2 + O(h^2),

the one-dimensional Schrodinger equation

    -u''(x) + V(x) u(x) = 2 E u(x)

becomes the tridiagonal eigenvalue problem

    ( 2/h^2 + V_i ) u_i - u_{i+1}/h^2 - u_{i-1}/h^2 = 2 E u_i.

The boundary conditions u(R_min) = u(R_max) = 0 are imposed by simply leaving
the endpoints out of the matrix, which therefore has dimension N-1.

For the harmonic oscillator, V(x) = x^2 in the units of the text, the exact
eigenvalues are 2E_k = 2k + 1, k = 0, 1, 2, ...

Author: Morten Hjorth-Jensen
"""

import numpy as np

from householder import SymmetricEigenSolver
from lanczos import Lanczos


class SchrodingerDiagonalization:
    """Discretise -u'' + V(x) u = 2E u and diagonalise the result."""

    def __init__(self, potential, rmin=-10.0, rmax=10.0, nsteps=400):
        self.potential = potential
        self.rmin = rmin
        self.rmax = rmax
        self.nsteps = nsteps
        self.h = (rmax - rmin) / nsteps
        # interior mesh points only; u vanishes at the two endpoints
        self.x = rmin + self.h * np.arange(1, nsteps)
        self.n = len(self.x)

    # ------------------------------------------------------------------
    @property
    def diagonal(self):
        return 2.0 / self.h**2 + self.potential(self.x)

    @property
    def offdiagonal(self):
        return np.full(self.n - 1, -1.0 / self.h**2)

    def matrix(self):
        """The tridiagonal matrix in dense form."""
        return (np.diag(self.diagonal) +
                np.diag(self.offdiagonal, 1) +
                np.diag(self.offdiagonal, -1))

    def matvec(self, v):
        """Apply the Hamiltonian without ever forming the matrix."""
        w = self.diagonal * v
        w[:-1] -= v[1:] / self.h**2
        w[1:] -= v[:-1] / self.h**2
        return w

    # ------------------------------------------------------------------
    def solve(self, k=5):
        """The k lowest eigenvalues and normalised eigenfunctions.

        The matrix is already tridiagonal, so the Householder step is
        unnecessary and we go straight to the QL algorithm.
        """
        values, vectors = SymmetricEigenSolver.tqli(self.diagonal,
                                                    self.offdiagonal)
        order = np.argsort(values)
        values = values[order][:k]
        vectors = vectors[:, order][:, :k]
        vectors = self._normalise(vectors)
        return values, vectors

    def solve_lanczos(self, k=5, m=None):
        """The same eigenvalues from the Lanczos algorithm.

        Only the product H v is used, so this is the route one follows when
        the matrix is far too large to store.
        """
        m = 10 * k if m is None else m
        lan = Lanczos(matvec=self.matvec, n=self.n)
        values, vectors = lan.ritz(m, k=k)
        return values, self._normalise(vectors)

    def _normalise(self, vectors):
        """Normalise with the trapezoidal rule, h * sum |u_i|^2 = 1."""
        norms = np.sqrt(self.h * np.sum(vectors**2, axis=0))
        vectors = vectors / norms
        # fix the overall sign so that the functions start out positive
        for j in range(vectors.shape[1]):
            first = vectors[np.argmax(np.abs(vectors[:, j])), j]
            if first < 0:
                vectors[:, j] *= -1.0
        return vectors


# ---------------------------------------------------------------------------
def harmonic_oscillator(x):
    """V(x) = x^2 in units where k = hbar = m = alpha = 1."""
    return x**2


def _demo():
    print("=" * 70)
    print("One-dimensional harmonic oscillator by diagonalisation")
    print("exact eigenvalues of the rewritten equation: 2E_k = 2k + 1")
    print("=" * 70)
    print(f"{'N':>7s} {'2E_0':>16s} {'2E_1':>16s} {'2E_2':>16s} "
          f"{'max error':>12s}")
    print("-" * 70)

    exact = np.array([1.0, 3.0, 5.0, 7.0, 9.0])
    for nsteps in (100, 200, 400, 800):
        problem = SchrodingerDiagonalization(harmonic_oscillator,
                                             -10.0, 10.0, nsteps)
        values, _ = problem.solve(k=5)
        err = np.max(np.abs(values - exact))
        print(f"{nsteps:7d} {values[0]:16.10f} {values[1]:16.10f} "
              f"{values[2]:16.10f} {err:12.2e}")
    print("-" * 70)
    print(f"{'exact':>7s} {exact[0]:16.10f} {exact[1]:16.10f} "
          f"{exact[2]:16.10f}")
    print()
    print("The error falls off as h^2, as it must for the three-point")
    print("formula for the second derivative.")

    print()
    print("=" * 70)
    print("The same eigenvalues from the Lanczos algorithm")
    print("=" * 70)
    problem = SchrodingerDiagonalization(harmonic_oscillator,
                                         -10.0, 10.0, 400)
    direct, u_direct = problem.solve(k=3)
    ritz, u_lanczos = problem.solve_lanczos(k=3, m=200)
    print(f"{'state':>7s} {'QL (direct)':>18s} {'Lanczos':>18s} "
          f"{'difference':>14s}")
    for j in range(3):
        print(f"{j:7d} {direct[j]:18.10f} {ritz[j]:18.10f} "
              f"{abs(direct[j]-ritz[j]):14.2e}")

    print()
    print("=" * 70)
    print("Check against the analytic wave functions")
    print("=" * 70)
    x = problem.x
    # u_0(x) = pi^{-1/4} exp(-x^2/2),  u_1(x) = (4/pi)^{1/4} x exp(-x^2/2)
    u0 = np.pi**-0.25 * np.exp(-x**2 / 2.0)
    u1 = (4.0 / np.pi)**0.25 * x * np.exp(-x**2 / 2.0)
    print(f"max |u_0 numerical - u_0 exact| = "
          f"{np.max(np.abs(np.abs(u_direct[:, 0]) - np.abs(u0))):.3e}")
    print(f"max |u_1 numerical - u_1 exact| = "
          f"{np.max(np.abs(np.abs(u_direct[:, 1]) - np.abs(u1))):.3e}")


if __name__ == "__main__":
    _demo()
