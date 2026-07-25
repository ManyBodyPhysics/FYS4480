"""
Householder tridiagonalisation and the QL algorithm with implicit shifts.

Companion code to chapter 1 of *Quantum mechanics for Many-particle Systems*.

The eigenvalue problem A x = lambda x for a real symmetric matrix is solved in
the two steps described in the text:

  1. a sequence of Householder similarity transformations reduces A to
     tridiagonal form,   T = S^T A S,   with S orthogonal;
  2. the tridiagonal matrix is diagonalised by the QL algorithm with implicit
     shifts, which is the modern replacement for Francis' original scheme.

Both steps are similarity transformations, so the eigenvalues are those of the
original matrix, and the eigenvectors are recovered by accumulating the
transformations.

Author: Morten Hjorth-Jensen
"""

import numpy as np


class SymmetricEigenSolver:
    """Eigenvalues and eigenvectors of a real symmetric matrix."""

    def __init__(self, A, tol=1.0e-14):
        A = np.asarray(A, dtype=float)
        if A.ndim != 2 or A.shape[0] != A.shape[1]:
            raise ValueError("A must be a square matrix")
        if not np.allclose(A, A.T, atol=1.0e-10):
            raise ValueError("A must be symmetric")
        self.A = A
        self.n = A.shape[0]
        self.tol = tol

    # ------------------------------------------------------------------
    # Step 1:  Householder reduction to tridiagonal form
    # ------------------------------------------------------------------
    def tridiagonalize(self):
        """Return (d, e, S) with S^T A S tridiagonal.

        d holds the n diagonal elements and e the n-1 off-diagonal ones.
        Each step uses the reflector P = I - 2 u u^T acting on the trailing
        submatrix, chosen so that the first column below the diagonal is
        annihilated except for one element.
        """
        n = self.n
        M = self.A.copy()
        S = np.eye(n)

        for k in range(n - 2):
            v = M[k + 1:, k]                     # the column to be reduced
            norm_v = np.linalg.norm(v)
            if norm_v < self.tol:
                continue
            # the sign is chosen to avoid cancellation in v[0] - kappa
            kappa = -np.copysign(norm_v, v[0])
            w = v.copy()
            w[0] -= kappa
            norm_w = np.linalg.norm(w)
            if norm_w < self.tol:
                continue
            u = w / norm_w                       # unit vector, u^T u = 1

            # similarity transformation of the trailing block,
            #   P A P = A - 2 z u^T - 2 u z^T,   z = A u - (u^T A u) u
            block = M[k + 1:, k + 1:]
            p = block @ u
            z = p - (u @ p) * u
            M[k + 1:, k + 1:] = block - 2.0 * np.outer(z, u) - 2.0 * np.outer(u, z)

            # the reduced column and row
            M[k + 1:, k] = 0.0
            M[k + 1, k] = kappa
            M[k, k + 1:] = 0.0
            M[k, k + 1] = kappa

            # accumulate S <- S P
            S[:, k + 1:] -= 2.0 * np.outer(S[:, k + 1:] @ u, u)

        d = np.diag(M).copy()
        e = np.diag(M, -1).copy()
        return d, e, S

    # ------------------------------------------------------------------
    # Step 2:  QL algorithm with implicit shifts for the tridiagonal matrix
    # ------------------------------------------------------------------
    @staticmethod
    def tqli(d, e, z=None, max_sweeps=50):
        """Diagonalise a symmetric tridiagonal matrix.

        Parameters
        ----------
        d : diagonal elements, length n
        e : off-diagonal elements, length n-1
        z : matrix whose columns are accumulated; pass the S of the
            Householder step to obtain the eigenvectors of the original A,
            or the identity to obtain those of the tridiagonal matrix.

        Returns (eigenvalues, eigenvectors) with the eigenvector belonging to
        eigenvalue k in column k.
        """
        d = np.asarray(d, dtype=float).copy()
        n = len(d)
        # pad so that e[i] couples i and i+1, with e[n-1] = 0
        e = np.concatenate([np.asarray(e, dtype=float), [0.0]])
        z = np.eye(n) if z is None else np.asarray(z, dtype=float).copy()
        eps = np.finfo(float).eps

        for l in range(n):
            for _ in range(max_sweeps):
                # look for a negligible off-diagonal element to split off
                m = n - 1
                for i in range(l, n - 1):
                    if abs(e[i]) <= eps * (abs(d[i]) + abs(d[i + 1])):
                        m = i
                        break
                if m == l:                       # d[l] has converged
                    break

                # Wilkinson shift, computed without cancellation
                g = (d[l + 1] - d[l]) / (2.0 * e[l])
                r = np.hypot(g, 1.0)
                g = d[m] - d[l] + e[l] / (g + np.copysign(r, g))

                s = c = 1.0
                p = 0.0
                broke = False
                for i in range(m - 1, l - 1, -1):
                    f = s * e[i]
                    b = c * e[i]
                    r = np.hypot(f, g)
                    e[i + 1] = r
                    if r == 0.0:                 # recover from underflow
                        d[i + 1] -= p
                        e[m] = 0.0
                        broke = True
                        break
                    s = f / r
                    c = g / r
                    g = d[i + 1] - p
                    r = (d[i] - g) * s + 2.0 * c * b
                    p = s * r
                    d[i + 1] = g + p
                    g = c * r - b
                    # accumulate the plane rotation
                    col = z[:, i + 1].copy()
                    z[:, i + 1] = s * z[:, i] + c * col
                    z[:, i] = c * z[:, i] - s * col
                if broke:
                    continue
                d[l] -= p
                e[l] = g
                e[m] = 0.0
            else:
                raise RuntimeError("too many iterations in tqli")

        return d, z

    # ------------------------------------------------------------------
    def eigen(self, sort=True):
        """Full solution: eigenvalues and eigenvectors of A."""
        d, e, S = self.tridiagonalize()
        values, vectors = self.tqli(d, e, S)
        if sort:
            order = np.argsort(values)
            values = values[order]
            vectors = vectors[:, order]
        return values, vectors

    def eigenvalues(self, sort=True):
        return self.eigen(sort)[0]


# ---------------------------------------------------------------------------
class PowerMethod:
    """The power method for the dominant eigenpair.

    Repeated multiplication by A amplifies the component along the eigenvector
    of largest |lambda|, since

        A^k b_0 = c_1 lambda_1^k ( v_1 + sum_{j>1} (c_j/c_1)
                                        (lambda_j/lambda_1)^k v_j ),

    and the bracket tends to v_1.  Convergence is geometric with ratio
    |lambda_2/lambda_1|, so the method is slow when the two largest
    eigenvalues are close.  It is the ancestor of the Krylov methods.
    """

    def __init__(self, A, tol=1.0e-12, max_iter=10000):
        self.A = np.asarray(A, dtype=float)
        self.tol = tol
        self.max_iter = max_iter
        self.iterations = 0

    def solve(self, b0=None):
        n = self.A.shape[0]
        b = np.ones(n) if b0 is None else np.asarray(b0, dtype=float).copy()
        b /= np.linalg.norm(b)
        mu_old = np.inf
        for k in range(self.max_iter):
            w = self.A @ b
            b = w / np.linalg.norm(w)
            mu = b @ (self.A @ b)          # Rayleigh quotient
            self.iterations = k + 1
            if abs(mu - mu_old) < self.tol:
                break
            mu_old = mu
        return mu, b


# ---------------------------------------------------------------------------
def _demo():
    rng = np.random.default_rng(2026)
    n = 8
    M = rng.normal(size=(n, n))
    A = 0.5 * (M + M.T)                    # a random symmetric matrix

    solver = SymmetricEigenSolver(A)

    print("=" * 68)
    print("Step 1: Householder reduction to tridiagonal form")
    print("=" * 68)
    d, e, S = solver.tridiagonalize()
    T = np.diag(d) + np.diag(e, 1) + np.diag(e, -1)
    print(f"|S^T S - I|        = {np.linalg.norm(S.T @ S - np.eye(n)):.3e}")
    print(f"|S T S^T - A|      = {np.linalg.norm(S @ T @ S.T - A):.3e}")
    print("diagonal      d =", np.array2string(d, precision=4))
    print("off-diagonal  e =", np.array2string(e, precision=4))

    print()
    print("=" * 68)
    print("Step 2: QL with implicit shifts")
    print("=" * 68)
    values, vectors = solver.eigen()
    reference = np.linalg.eigvalsh(A)
    print("eigenvalues (Householder + QL):")
    print("  ", np.array2string(values, precision=8))
    print("eigenvalues (numpy.linalg.eigvalsh):")
    print("  ", np.array2string(reference, precision=8))
    print(f"max deviation      = {np.max(np.abs(values - reference)):.3e}")
    print(f"|A V - V diag(l)|  = "
          f"{np.linalg.norm(A @ vectors - vectors @ np.diag(values)):.3e}")
    print(f"|V^T V - I|        = "
          f"{np.linalg.norm(vectors.T @ vectors - np.eye(n)):.3e}")

    print()
    print("=" * 68)
    print("The power method finds only the dominant eigenvalue")
    print("=" * 68)
    power = PowerMethod(A)
    mu, v = power.solve()
    dominant = reference[np.argmax(np.abs(reference))]
    ratio = np.sort(np.abs(reference))[-2] / np.abs(dominant)
    print(f"  lambda (power method) = {mu:.10f}")
    print(f"  lambda (exact)        = {dominant:.10f}")
    print(f"  iterations            = {power.iterations}")
    print(f"  |lambda_2/lambda_1|   = {ratio:.4f}  (convergence ratio)")
    print(f"  |A v - lambda v|      = "
          f"{np.linalg.norm(A @ v - mu * v):.3e}")


if __name__ == "__main__":
    _demo()
