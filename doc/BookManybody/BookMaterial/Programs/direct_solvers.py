"""
Direct solvers for linear systems of equations.

Companion code to chapter 1 of *Quantum mechanics for Many-particle Systems*.

The classes implement, in an object-oriented style, the direct (elimination)
methods discussed in the text:

    GaussianElimination  -- Gaussian elimination with partial pivoting
    LUDecomposition      -- Doolittle LU factorisation with partial pivoting
    CholeskyDecomposition-- A = L L^T for symmetric positive definite A
    TridiagonalSolver    -- the O(n) forward/backward substitution algorithm

All indices follow the convention of the book: rows and columns are numbered
from 0 to n-1, exactly as in Python and NumPy.

Author: Morten Hjorth-Jensen
"""

import numpy as np


# ---------------------------------------------------------------------------
class DirectSolver:
    """Common interface for the direct solvers of A x = b."""

    def __init__(self, A):
        A = np.asarray(A, dtype=float)
        if A.ndim != 2 or A.shape[0] != A.shape[1]:
            raise ValueError("A must be a square matrix")
        self.A = A
        self.n = A.shape[0]

    # -- to be provided by the subclasses ----------------------------------
    def solve(self, b):
        raise NotImplementedError

    def determinant(self):
        raise NotImplementedError

    # -- generic services ---------------------------------------------------
    def inverse(self):
        """Column j of A^{-1} is the solution of A x = e_j."""
        columns = [self.solve(e) for e in np.eye(self.n)]
        return np.column_stack(columns)

    def residual(self, x, b):
        """The residual r = b - A x, a cheap a posteriori error estimate."""
        return np.asarray(b, dtype=float) - self.A @ x


# ---------------------------------------------------------------------------
class GaussianElimination(DirectSolver):
    """Gaussian elimination with partial pivoting.

    The matrix is reduced to upper triangular form by row operations, after
    which the solution is obtained by backward substitution.  The cost is
    2n^3/3 + O(n^2) floating point operations.
    """

    def solve(self, b):
        n = self.n
        # work on copies: the elimination destroys its input
        M = self.A.copy()
        y = np.asarray(b, dtype=float).copy()

        # forward elimination
        for k in range(n - 1):
            # partial pivoting: use the largest element in the column as pivot
            p = k + np.argmax(np.abs(M[k:, k]))
            if abs(M[p, k]) < 1.0e-14:
                raise np.linalg.LinAlgError("matrix is singular to machine precision")
            if p != k:
                M[[k, p]] = M[[p, k]]
                y[k], y[p] = y[p], y[k]
            for i in range(k + 1, n):
                factor = M[i, k] / M[k, k]
                M[i, k:] -= factor * M[k, k:]
                y[i] -= factor * y[k]

        return self._back_substitute(M, y)

    @staticmethod
    def _back_substitute(U, y):
        """Solve U x = y for an upper triangular U."""
        n = U.shape[0]
        x = np.zeros(n)
        for m in range(n - 1, -1, -1):
            x[m] = (y[m] - U[m, m + 1:] @ x[m + 1:]) / U[m, m]
        return x

    def determinant(self):
        n = self.n
        M = self.A.copy()
        sign = 1.0
        for k in range(n - 1):
            p = k + np.argmax(np.abs(M[k:, k]))
            if abs(M[p, k]) < 1.0e-14:
                return 0.0
            if p != k:
                M[[k, p]] = M[[p, k]]
                sign = -sign
            for i in range(k + 1, n):
                M[i, k:] -= (M[i, k] / M[k, k]) * M[k, k:]
        return sign * np.prod(np.diag(M))


# ---------------------------------------------------------------------------
class LUDecomposition(DirectSolver):
    """Doolittle LU factorisation with partial pivoting, P A = L U.

    L is unit lower triangular and U upper triangular.  The factorisation is
    computed once in the constructor and can then be reused for any number of
    right-hand sides, which is what makes LU the method of choice when the
    same matrix is used repeatedly -- as in the inverse, where n systems have
    to be solved.
    """

    def __init__(self, A):
        super().__init__(A)
        self._factorize()

    def _factorize(self):
        n = self.n
        LU = self.A.copy()
        perm = np.arange(n)
        sign = 1.0

        for k in range(n):
            p = k + np.argmax(np.abs(LU[k:, k]))
            if abs(LU[p, k]) < 1.0e-14:
                raise np.linalg.LinAlgError("matrix is singular to machine precision")
            if p != k:
                LU[[k, p]] = LU[[p, k]]
                perm[[k, p]] = perm[[p, k]]
                sign = -sign
            # the multipliers l_ik are stored in place, below the diagonal
            LU[k + 1:, k] /= LU[k, k]
            # rank-one update of the trailing submatrix
            LU[k + 1:, k + 1:] -= np.outer(LU[k + 1:, k], LU[k, k + 1:])

        self.LU = LU
        self.perm = perm
        self.sign = sign

    # -- the two triangular factors, for inspection -------------------------
    @property
    def L(self):
        return np.tril(self.LU, -1) + np.eye(self.n)

    @property
    def U(self):
        return np.triu(self.LU)

    def solve(self, b):
        """Solve A x = b in two steps: L y = P b, then U x = y."""
        y = np.asarray(b, dtype=float)[self.perm].copy()
        n = self.n
        # forward substitution, L has unit diagonal
        for i in range(1, n):
            y[i] -= self.LU[i, :i] @ y[:i]
        # backward substitution
        x = np.zeros(n)
        for i in range(n - 1, -1, -1):
            x[i] = (y[i] - self.LU[i, i + 1:] @ x[i + 1:]) / self.LU[i, i]
        return x

    def determinant(self):
        """det(A) = sign * u_00 u_11 ... u_{n-1,n-1}."""
        return self.sign * np.prod(np.diag(self.LU))


# ---------------------------------------------------------------------------
class CholeskyDecomposition(DirectSolver):
    """Cholesky factorisation A = L L^T of a symmetric positive definite A.

    Half the work of a general LU factorisation and, since the pivots are
    guaranteed positive, no pivoting is needed.  This is the algorithm of
    choice for the overlap and Hessian matrices met later in the book.
    """

    def __init__(self, A):
        super().__init__(A)
        if not np.allclose(self.A, self.A.T):
            raise ValueError("Cholesky requires a symmetric matrix")
        self._factorize()

    def _factorize(self):
        n = self.n
        L = np.zeros((n, n))
        for i in range(n):
            s = self.A[i, i] - L[i, :i] @ L[i, :i]
            if s <= 0.0:
                raise np.linalg.LinAlgError("matrix is not positive definite")
            L[i, i] = np.sqrt(s)
            for j in range(i + 1, n):
                L[j, i] = (self.A[j, i] - L[j, :i] @ L[i, :i]) / L[i, i]
        self.L = L

    def solve(self, b):
        n = self.n
        y = np.asarray(b, dtype=float).copy()
        for i in range(n):                      # forward:  L y = b
            y[i] = (y[i] - self.L[i, :i] @ y[:i]) / self.L[i, i]
        x = np.zeros(n)
        for i in range(n - 1, -1, -1):          # backward: L^T x = y
            x[i] = (y[i] - self.L[i + 1:, i] @ x[i + 1:]) / self.L[i, i]
        return x

    def determinant(self):
        return float(np.prod(np.diag(self.L)) ** 2)


# ---------------------------------------------------------------------------
class TridiagonalSolver:
    """Solve a tridiagonal system in O(n) operations.

    The system is

        a_i u_{i-1} + b_i u_i + c_i u_{i+1} = f_i,     i = 0, ..., n-1,

    with a_0 and c_{n-1} not referenced.  Gaussian elimination on the same
    problem would cost 2n^3/3 operations, so for the discretised differential
    equations of the text this specialised algorithm is the only sensible
    choice.
    """

    def __init__(self, a, b, c):
        self.a = np.asarray(a, dtype=float)   # subdiagonal,   a[0] unused
        self.b = np.asarray(b, dtype=float)   # diagonal
        self.c = np.asarray(c, dtype=float)   # superdiagonal, c[n-1] unused
        self.n = len(self.b)

    def solve(self, f):
        n = self.n
        f = np.asarray(f, dtype=float)
        temp = np.zeros(n)
        u = np.zeros(n)

        # forward substitution
        btemp = self.b[0]
        if abs(btemp) < 1.0e-14:
            raise np.linalg.LinAlgError("b[0] vanishes; reorder the equations")
        u[0] = f[0] / btemp
        for i in range(1, n):
            temp[i] = self.c[i - 1] / btemp
            btemp = self.b[i] - self.a[i] * temp[i]
            u[i] = (f[i] - self.a[i] * u[i - 1]) / btemp

        # backward substitution
        for i in range(n - 2, -1, -1):
            u[i] -= temp[i + 1] * u[i + 1]
        return u

    def dense(self):
        """The same matrix in dense form, for testing."""
        A = np.diag(self.b)
        A += np.diag(self.a[1:], -1)
        A += np.diag(self.c[:-1], 1)
        return A


# ---------------------------------------------------------------------------
def _demo():
    rng = np.random.default_rng(2026)
    n = 6

    print("=" * 62)
    print("Gaussian elimination, LU and the inverse")
    print("=" * 62)
    A = rng.normal(size=(n, n)) + n * np.eye(n)
    b = rng.normal(size=n)
    exact = np.linalg.solve(A, b)

    for cls in (GaussianElimination, LUDecomposition):
        solver = cls(A)
        x = solver.solve(b)
        print(f"{cls.__name__:22s} |x-x_exact| = {np.linalg.norm(x-exact):.3e}"
              f"   det = {solver.determinant():+.6e}")
    print(f"{'numpy reference':22s} "
          f"{'':22s}   det = {np.linalg.det(A):+.6e}")

    lu = LUDecomposition(A)
    print(f"|A^-1 A - I|            = "
          f"{np.linalg.norm(lu.inverse() @ A - np.eye(n)):.3e}")

    print()
    print("=" * 62)
    print("Cholesky for a symmetric positive definite matrix")
    print("=" * 62)
    M = rng.normal(size=(n, n))
    S = M @ M.T + n * np.eye(n)
    chol = CholeskyDecomposition(S)
    x = chol.solve(b)
    print(f"|L L^T - A|             = "
          f"{np.linalg.norm(chol.L @ chol.L.T - S):.3e}")
    print(f"|x - x_exact|           = "
          f"{np.linalg.norm(x - np.linalg.solve(S, b)):.3e}")

    print()
    print("=" * 62)
    print("Tridiagonal solver: -u'' = f with u(0) = u(1) = 0")
    print("=" * 62)
    # exact solution u(x) = x(1-x), so that f(x) = 2
    n = 100
    h = 1.0 / (n + 1)
    x = np.linspace(h, 1.0 - h, n)
    a = np.full(n, -1.0 / h**2)
    bdiag = np.full(n, 2.0 / h**2)
    c = np.full(n, -1.0 / h**2)
    f = np.full(n, 2.0)
    u = TridiagonalSolver(a, bdiag, c).solve(f)
    print(f"n = {n}, max |u_numerical - x(1-x)| = "
          f"{np.max(np.abs(u - x * (1.0 - x))):.3e}")


if __name__ == "__main__":
    _demo()
