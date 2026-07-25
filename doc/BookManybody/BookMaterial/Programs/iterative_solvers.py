"""
Iterative solvers for linear systems of equations.

Companion code to chapter 1 of *Quantum mechanics for Many-particle Systems*.

    Jacobi            -- x^{k+1} = D^{-1} (b - (L+U) x^k)
    GaussSeidel       -- uses the already updated components immediately
    SOR               -- Gauss-Seidel with a relaxation parameter omega
    ConjugateGradient -- Krylov method for symmetric positive definite A

Iterative methods never touch the matrix itself, only the product A x.  That
is why they, and not the direct methods, are what one uses for the very large
sparse Hamiltonian matrices of full configuration-interaction calculations.

Author: Morten Hjorth-Jensen
"""

import numpy as np


# ---------------------------------------------------------------------------
class IterativeSolver:
    """Common interface and convergence bookkeeping."""

    def __init__(self, A, tol=1.0e-10, max_iter=10000):
        A = np.asarray(A, dtype=float)
        if A.ndim != 2 or A.shape[0] != A.shape[1]:
            raise ValueError("A must be a square matrix")
        self.A = A
        self.n = A.shape[0]
        self.tol = tol
        self.max_iter = max_iter
        self.iterations = 0
        self.history = []          # ||r_k|| at every iteration

    def solve(self, b, x0=None):
        raise NotImplementedError

    def _initial_guess(self, x0):
        return np.zeros(self.n) if x0 is None else np.asarray(x0, float).copy()

    def _converged(self, r):
        norm = np.linalg.norm(r)
        self.history.append(norm)
        return norm < self.tol


# ---------------------------------------------------------------------------
class Jacobi(IterativeSolver):
    """The Jacobi iteration.

    With A = D + L + U the update reads

        x^{k+1} = D^{-1} ( b - (L + U) x^k ),

    that is, every component of the new iterate is computed from the *old*
    ones only.  This makes the method trivially parallel, but slower than
    Gauss-Seidel.  It converges whenever A is strictly diagonally dominant or
    symmetric positive definite.
    """

    def solve(self, b, x0=None):
        b = np.asarray(b, dtype=float)
        x = self._initial_guess(x0)
        d = np.diag(self.A)
        if np.any(np.abs(d) < 1.0e-14):
            raise ValueError("zero on the diagonal; reorder the equations")
        R = self.A - np.diag(d)                 # R = L + U
        self.history = []
        for k in range(self.max_iter):
            x_new = (b - R @ x) / d
            self.iterations = k + 1
            if self._converged(b - self.A @ x_new):
                return x_new
            x = x_new
        return x


# ---------------------------------------------------------------------------
class GaussSeidel(IterativeSolver):
    """The Gauss-Seidel iteration,

        x_i^{k+1} = ( b_i - sum_{j<i} a_ij x_j^{k+1}
                          - sum_{j>i} a_ij x_j^{k}   ) / a_ii,

    i.e. the components already updated in the current sweep are used at once
    (forward substitution).  This typically halves the number of iterations
    compared with Jacobi, at the price of being inherently sequential.
    """

    def solve(self, b, x0=None):
        b = np.asarray(b, dtype=float)
        x = self._initial_guess(x0)
        A, n = self.A, self.n
        self.history = []
        for k in range(self.max_iter):
            for i in range(n):
                s = A[i, :i] @ x[:i] + A[i, i + 1:] @ x[i + 1:]
                x[i] = (b[i] - s) / A[i, i]
            self.iterations = k + 1
            if self._converged(b - A @ x):
                return x
        return x


# ---------------------------------------------------------------------------
class SOR(IterativeSolver):
    """Successive over-relaxation,

        x_i^{k+1} = (1-omega) x_i^k
                  + (omega/a_ii) ( b_i - sum_{j<i} a_ij x_j^{k+1}
                                       - sum_{j>i} a_ij x_j^{k} ).

    For symmetric positive definite matrices convergence is guaranteed for
    0 < omega < 2; omega = 1 gives back Gauss-Seidel.  A well chosen omega can
    reduce the iteration count by an order of magnitude.
    """

    def __init__(self, A, omega=1.5, tol=1.0e-10, max_iter=10000):
        super().__init__(A, tol, max_iter)
        if not 0.0 < omega < 2.0:
            raise ValueError("omega must lie in (0, 2)")
        self.omega = omega

    def solve(self, b, x0=None):
        b = np.asarray(b, dtype=float)
        x = self._initial_guess(x0)
        A, n, w = self.A, self.n, self.omega
        self.history = []
        for k in range(self.max_iter):
            for i in range(n):
                s = A[i, :i] @ x[:i] + A[i, i + 1:] @ x[i + 1:]
                x[i] = (1.0 - w) * x[i] + w * (b[i] - s) / A[i, i]
            self.iterations = k + 1
            if self._converged(b - A @ x):
                return x
        return x


# ---------------------------------------------------------------------------
class ConjugateGradient(IterativeSolver):
    """The conjugate gradient method for symmetric positive definite A.

    The solution of A x = b is the unique minimum of

        P(x) = (1/2) x^T A x - x^T b,

    and the method minimises P along a sequence of directions p_k that are
    conjugate, p_i^T A p_j = 0 for i != j.  In exact arithmetic it therefore
    terminates after at most n steps, but in practice a good approximation is
    reached in far fewer iterations.

    Only the matrix-vector product A p is required, so a `matvec` callable may
    be supplied instead of a dense matrix -- which is how the method is used
    for large sparse many-body matrices.
    """

    def __init__(self, A, tol=1.0e-10, max_iter=None, matvec=None, n=None):
        if matvec is not None:
            self.matvec = matvec
            self.n = n
            self.A = None
            self.tol = tol
            self.max_iter = n if max_iter is None else max_iter
            self.iterations = 0
            self.history = []
        else:
            super().__init__(A, tol, A.shape[0] if max_iter is None else max_iter)
            self.matvec = lambda v: self.A @ v

    def solve(self, b, x0=None):
        b = np.asarray(b, dtype=float)
        x = self._initial_guess(x0)
        r = b - self.matvec(x)            # residual
        p = r.copy()                      # first search direction
        rs = r @ r
        self.history = [np.sqrt(rs)]

        for k in range(self.max_iter):
            Ap = self.matvec(p)
            alpha = rs / (p @ Ap)
            x += alpha * p                # x_{k+1} = x_k + alpha_k p_k
            r -= alpha * Ap               # r_{k+1} = r_k - alpha_k A p_k
            rs_new = r @ r
            self.iterations = k + 1
            self.history.append(np.sqrt(rs_new))
            if np.sqrt(rs_new) < self.tol:
                break
            p = r + (rs_new / rs) * p     # beta_k = |r_{k+1}|^2 / |r_k|^2
            rs = rs_new
        return x


# ---------------------------------------------------------------------------
def _demo():
    # the discretised second derivative: symmetric, positive definite,
    # strictly diagonally dominant -- every method below converges on it
    n = 50
    h = 1.0 / (n + 1)
    A = (np.diag(np.full(n, 2.0)) +
         np.diag(np.full(n - 1, -1.0), 1) +
         np.diag(np.full(n - 1, -1.0), -1)) / h**2
    x_grid = np.linspace(h, 1.0 - h, n)
    b = np.full(n, 2.0)                    # exact solution x(1-x)
    exact = x_grid * (1.0 - x_grid)

    print("=" * 66)
    print(f"-u'' = 2 on (0,1), u(0)=u(1)=0, n = {n} interior points")
    print("=" * 66)
    print(f"{'method':22s} {'iterations':>12s} {'|x - exact|':>16s}")
    print("-" * 66)

    for solver in (Jacobi(A, tol=1e-8, max_iter=100000),
                   GaussSeidel(A, tol=1e-8, max_iter=100000),
                   SOR(A, omega=1.9, tol=1e-8, max_iter=100000),
                   ConjugateGradient(A, tol=1e-8)):
        x = solver.solve(b)
        name = type(solver).__name__
        if isinstance(solver, SOR):
            name += f" (omega={solver.omega})"
        print(f"{name:22s} {solver.iterations:12d} "
              f"{np.linalg.norm(x - exact):16.3e}")

    print()
    print("The conjugate gradient method needs only the product A p, so a")
    print("matrix-free version costs no storage at all:")
    matvec = lambda v: (np.concatenate(([0.0], v[:-1])) - 2.0 * v +
                        np.concatenate((v[1:], [0.0]))) / (-h**2)
    cg = ConjugateGradient(None, tol=1e-8, matvec=matvec, n=n)
    x = cg.solve(b)
    print(f"matrix-free CG: {cg.iterations} iterations, "
          f"|x - exact| = {np.linalg.norm(x - exact):.3e}")


if __name__ == "__main__":
    _demo()
