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
def pairing_hamiltonian(levels=12, pairs=6, delta=1.0, g=0.5):
    """The pairing model Hamiltonian in the seniority-zero space.

    We use the standard pairing Hamiltonian of the course,

        H = delta * sum_p (p-1) N_p - g * sum_{pq} P_p^+ P_q,

    with doubly degenerate single-particle levels p = 1, ..., L and N_p the
    number operator.  In the space of states with no broken pairs a basis
    state is a choice of which levels carry a pair, so the dimension is the
    binomial coefficient C(L, n).  The diagonal elements are the unperturbed
    energies minus g (from the p = q terms) and any two configurations that
    differ by moving a single pair are coupled by -g.

    This is a genuine many-body matrix: sparse, with a well separated ground
    state, and it is exactly the kind of problem the Lanczos algorithm was
    designed for.
    """
    from itertools import combinations

    basis = list(combinations(range(1, levels + 1), pairs))
    index = {c: i for i, c in enumerate(basis)}
    dim = len(basis)
    H = np.zeros((dim, dim))

    for c, i in index.items():
        H[i, i] = 2.0 * delta * sum(p - 1 for p in c) - g * pairs
        occupied = set(c)
        for p in c:                                  # move the pair from p
            for q in range(1, levels + 1):           # to an empty level q
                if q in occupied:
                    continue
                new = tuple(sorted(occupied - {p} | {q}))
                H[index[new], i] -= g
    return H, basis


def _demo():
    H, basis = pairing_hamiltonian(levels=12, pairs=6, delta=1.0, g=0.5)
    n = H.shape[0]
    exact = np.linalg.eigvalsh(H)

    print("=" * 72)
    print(f"Lanczos for the pairing model: 12 levels, 6 pairs, "
          f"dimension {n}")
    print("=" * 72)
    print(f"{'m':>5s} {'lowest Ritz value':>22s} {'error':>12s} "
          f"{'second':>18s} {'error':>12s}")
    print("-" * 72)
    for m in (5, 10, 20, 30, 40):
        theta, _ = Lanczos(H).ritz(m)
        print(f"{m:5d} {theta[0]:22.12f} {abs(theta[0]-exact[0]):12.2e} "
              f"{theta[1]:18.10f} {abs(theta[1]-exact[1]):12.2e}")
    print("-" * 72)
    print(f"{'exact':>5s} {exact[0]:22.12f} {'':12s} {exact[1]:18.10f}")

    print()
    print(f"A few tens of matrix-vector products give the ground-state energy")
    print(f"of a {n}x{n} many-body matrix to machine precision, and the")
    print("Lanczos vectors are the only large objects ever stored.")

    print()
    print("=" * 72)
    print("Loss of orthogonality without reorthogonalisation")
    print("=" * 72)
    for flag in (True, False):
        lan = Lanczos(H, reorthogonalize=flag)
        _, _, Q = lan.run(120)
        dev = np.linalg.norm(Q.T @ Q - np.eye(Q.shape[1]))
        theta, _ = lan.ritz(120)
        # the ground state is non-degenerate, so any extra Ritz value sitting
        # on top of it is a spurious "ghost" produced by the loss of
        # orthogonality among the Lanczos vectors
        copies = int(np.sum(np.abs(theta - exact[0]) < 1.0e-8))
        print(f"reorthogonalize = {str(flag):5s}:  |Q^T Q - I| = {dev:9.2e},"
              f"  copies of the ground state after 120 steps: {copies}")

    print()
    print("=" * 72)
    print("Matrix-free Lanczos: only the action of H on a vector is needed")
    print("=" * 72)
    theta, vec = Lanczos(matvec=lambda v: H @ v, n=n).ritz(40)
    print(f"lowest Ritz value = {theta[0]:.12f}   (exact {exact[0]:.12f})")
    print(f"|H v - theta v|   = "
          f"{np.linalg.norm(H @ vec[:, 0] - theta[0] * vec[:, 0]):.2e}")


if __name__ == "__main__":
    _demo()
