"""
The Lanczos algorithm for large sparse symmetric eigenvalue problems.

Companion code to chapter 1 of *Quantum mechanics for Many-particle Systems*.

The Lanczos algorithm builds an orthonormal basis of the Krylov space

    K_m(A, q_0) = span{ q_0, A q_0, A^2 q_0, ..., A^{m-1} q_0 }

by the three-term recurrence

    A q_k = beta_{k-1} q_{k-1} + alpha_k q_k + beta_k q_{k+1},

in which the matrix A enters only through the product A q.  In that basis A is
represented by a small tridiagonal matrix T_m whose extremal eigenvalues (the
Ritz values) converge very quickly to the extremal eigenvalues of A.  This is
the method used to obtain ground-state energies in large-scale shell-model and
full configuration-interaction calculations, where the dimension of A can
exceed 10^9 and A is never stored as a dense matrix.

Author: Morten Hjorth-Jensen
"""

import numpy as np

from householder import SymmetricEigenSolver


class Lanczos:
    """Lanczos iteration for a symmetric matrix or a matrix-free operator.

    Parameters
    ----------
    A : array or None
        The matrix.  May be omitted if `matvec` is given.
    matvec : callable
        A function returning the product A v.  This is the interface used for
        the huge sparse Hamiltonians of many-body physics.
    n : int
        Dimension, required when only `matvec` is supplied.
    reorthogonalize : bool
        In exact arithmetic the Lanczos vectors are orthogonal, but in finite
        precision orthogonality is lost as soon as one Ritz value converges,
        producing spurious copies of already converged eigenvalues.  Full
        reorthogonalisation against all previous vectors cures this at the
        cost of storing them.
    """

    def __init__(self, A=None, matvec=None, n=None, reorthogonalize=True):
        if matvec is not None:
            if n is None:
                raise ValueError("n must be given together with matvec")
            self.matvec = matvec
            self.n = n
        else:
            A = np.asarray(A, dtype=float)
            if not np.allclose(A, A.T, atol=1.0e-10):
                raise ValueError("A must be symmetric")
            self.matvec = lambda v: A @ v
            self.n = A.shape[0]
        self.reorthogonalize = reorthogonalize

    # ------------------------------------------------------------------
    def run(self, m, q0=None, tol=1.0e-12):
        """Perform m Lanczos steps.

        Returns (alpha, beta, Q) where alpha holds the diagonal and beta the
        off-diagonal elements of T_m, and the columns of Q are the Lanczos
        vectors.
        """
        n = self.n
        m = min(m, n)
        rng = np.random.default_rng(2026)
        q = rng.normal(size=n) if q0 is None else np.asarray(q0, dtype=float).copy()
        q /= np.linalg.norm(q)

        Q = np.zeros((n, m))
        alpha = np.zeros(m)
        beta = np.zeros(max(m - 1, 0))

        q_prev = np.zeros(n)
        b_prev = 0.0

        for k in range(m):
            Q[:, k] = q
            w = self.matvec(q)
            alpha[k] = q @ w                       # alpha_k = q_k^T A q_k
            w = w - alpha[k] * q - b_prev * q_prev # r_k

            if self.reorthogonalize:
                # two passes of modified Gram-Schmidt are enough in practice
                for _ in range(2):
                    w -= Q[:, :k + 1] @ (Q[:, :k + 1].T @ w)

            if k == m - 1:
                break
            b = np.linalg.norm(w)
            if b < tol:                            # invariant subspace found
                alpha = alpha[:k + 1]
                beta = beta[:k]
                Q = Q[:, :k + 1]
                break
            beta[k] = b
            q_prev, b_prev = q, b
            q = w / b

        self.alpha, self.beta, self.Q = alpha, beta, Q
        return alpha, beta, Q

    # ------------------------------------------------------------------
    def ritz(self, m, q0=None, k=None):
        """Ritz values and vectors after m steps.

        The small tridiagonal matrix T_m is diagonalised with the Householder
        and QL machinery of `householder.py`; its eigenvalues are the Ritz
        values, and the Ritz vectors are Q y.
        """
        alpha, beta, Q = self.run(m, q0)
        T = np.diag(alpha) + np.diag(beta, 1) + np.diag(beta, -1)
        theta, y = SymmetricEigenSolver(T).eigen()
        vectors = Q @ y
        if k is not None:
            theta, vectors = theta[:k], vectors[:, :k]
        return theta, vectors

    def ground_state(self, m, q0=None):
        """The lowest Ritz value and its Ritz vector."""
        theta, vectors = self.ritz(m, q0)
        return theta[0], vectors[:, 0]


# ---------------------------------------------------------------------------
def _harmonic_oscillator_matrix(n=500, rmax=10.0):
    """The tridiagonal matrix of -d^2/dx^2 + x^2 on [-rmax, rmax].

    A physically realistic test case: the low-lying eigenvalues are well
    separated, which is precisely the situation in which Lanczos shines.  In
    the units of the text the eigenvalues are 2E_k = 2k + 1.
    """
    h = 2.0 * rmax / (n + 1)
    x = -rmax + h * np.arange(1, n + 1)
    diag = 2.0 / h**2 + x**2
    off = np.full(n - 1, -1.0 / h**2)
    return np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)


def _demo():
    n = 500
    A = _harmonic_oscillator_matrix(n)
    exact = np.linalg.eigvalsh(A)

    print("=" * 72)
    print(f"Lanczos for the one-dimensional harmonic oscillator, n = {n}")
    print("The exact lowest eigenvalues of this matrix are 1, 3, 5, 7, ...")
    print("=" * 72)
    print(f"{'m':>5s} {'lowest Ritz value':>22s} {'error':>12s} "
          f"{'second':>18s} {'error':>12s}")
    print("-" * 72)
    for m in (5, 10, 20, 30, 40):
        theta, _ = Lanczos(A).ritz(m)
        print(f"{m:5d} {theta[0]:22.12f} {abs(theta[0]-exact[0]):12.2e} "
              f"{theta[1]:18.10f} {abs(theta[1]-exact[1]):12.2e}")
    print("-" * 72)
    print(f"{'exact':>5s} {exact[0]:22.12f} {'':12s} {exact[1]:18.10f}")

    print()
    print("Forty matrix-vector products give the two lowest eigenvalues of a")
    print(f"{n}x{n} matrix to ten digits.  A full diagonalisation would cost")
    print(f"O(n^3) = {n**3:.1e} operations and require storing all of A.")

    print()
    print("=" * 72)
    print("Loss of orthogonality without reorthogonalisation")
    print("=" * 72)
    for flag in (True, False):
        lan = Lanczos(A, reorthogonalize=flag)
        _, _, Q = lan.run(120)
        dev = np.linalg.norm(Q.T @ Q - np.eye(Q.shape[1]))
        theta, _ = lan.ritz(120)
        # spurious ("ghost") eigenvalues appear as near-duplicates
        gaps = np.diff(np.sort(theta))
        ghosts = int(np.sum(gaps < 1.0e-8))
        print(f"reorthogonalize = {str(flag):5s}:  |Q^T Q - I| = {dev:9.2e},"
              f"  duplicated Ritz values after 120 steps: {ghosts}")

    print()
    print("=" * 72)
    print("Matrix-free Lanczos: only the action of A on a vector is needed")
    print("=" * 72)

    h = 2.0 * 10.0 / (n + 1)
    x = -10.0 + h * np.arange(1, n + 1)

    def matvec(v):
        """Apply -d^2/dx^2 + x^2 without ever forming the matrix."""
        w = (2.0 / h**2 + x**2) * v
        w[:-1] -= v[1:] / h**2
        w[1:] -= v[:-1] / h**2
        return w

    theta, vec = Lanczos(matvec=matvec, n=n).ritz(40)
    print(f"lowest Ritz value = {theta[0]:.12f}   (exact {exact[0]:.12f})")
    print(f"|A v - theta v|   = "
          f"{np.linalg.norm(A @ vec[:, 0] - theta[0] * vec[:, 0]):.2e}")


if __name__ == "__main__":
    _demo()
