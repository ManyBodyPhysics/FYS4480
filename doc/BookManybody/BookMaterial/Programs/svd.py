"""
The singular value decomposition and its many-body applications.

Companion code to chapter 1 of *Quantum mechanics for Many-particle Systems*.

The SVD writes any matrix, square or not, as

    X = U Sigma V^T,

with U and V orthogonal and Sigma diagonal and non-negative.  It always
exists, it is numerically stable, and in many-body physics it turns out to be
the same object as the Schmidt decomposition of a bipartite quantum state.
The classes below cover the four uses discussed in the text:

    SVD                 -- the decomposition itself, rank, pseudoinverse and
                           the optimal low-rank approximation
    CanonicalOrthogonalization
                        -- curing a near-singular overlap matrix
    SchmidtDecomposition-- bipartite entanglement, natural orbitals and
                           occupation numbers
    TwoBodyFactorization-- low-rank factorisation of two-body matrix elements

Author: Morten Hjorth-Jensen
"""

import numpy as np


# ---------------------------------------------------------------------------
class SVD:
    """The singular value decomposition of a general m x n matrix.

    We use the LAPACK driver through NumPy for the factorisation itself: a
    hand-written Golub-Kahan bidiagonalisation would add nothing to the
    discussion, and the whole point of the SVD is that a stable
    implementation already exists.  What the class adds is everything one
    actually does with the decomposition.
    """

    def __init__(self, X, full_matrices=False):
        self.X = np.asarray(X, dtype=float)
        self.m, self.n = self.X.shape
        self.U, self.sigma, self.Vt = np.linalg.svd(self.X,
                                                    full_matrices=full_matrices)

    # ------------------------------------------------------------------
    @property
    def V(self):
        return self.Vt.T

    def reconstruct(self, k=None):
        """The rank-k truncation U_k Sigma_k V_k^T.

        By the Eckart-Young theorem this is the best rank-k approximation of
        X in both the 2-norm and the Frobenius norm, and the error it makes
        is exactly the first discarded singular value,

            || X - X_k ||_2 = sigma_k.
        """
        k = len(self.sigma) if k is None else k
        return (self.U[:, :k] * self.sigma[:k]) @ self.Vt[:k]

    def rank(self, tol=None):
        """Numerical rank: the number of singular values above a threshold."""
        if tol is None:
            tol = max(self.m, self.n) * np.finfo(float).eps * self.sigma[0]
        return int(np.sum(self.sigma > tol))

    def condition_number(self, tol=None):
        """kappa_2(X) = sigma_max / sigma_min, infinite for a singular X."""
        if tol is None:
            tol = max(self.m, self.n) * np.finfo(float).eps
        if self.sigma[-1] <= tol * self.sigma[0]:
            return np.inf
        return self.sigma[0] / self.sigma[-1]

    def pseudoinverse(self, rcond=1.0e-12):
        """The Moore-Penrose pseudoinverse X^+ = V Sigma^+ U^T.

        Singular values below rcond * sigma_max are treated as zero rather
        than inverted, which is what makes the pseudoinverse usable on
        singular and near-singular matrices.
        """
        s_inv = np.zeros_like(self.sigma)
        keep = self.sigma > rcond * self.sigma[0]
        s_inv[keep] = 1.0 / self.sigma[keep]
        return (self.Vt.T * s_inv) @ self.U.T

    # ------------------------------------------------------------------
    def eigen_route(self):
        """Singular values obtained the naive way, from the eigenvalues of X^T X.

        The singular values are the positive square roots of the eigenvalues
        of X^T X, so this *is* correct mathematically.  It is a bad idea
        numerically: forming X^T X squares the condition number, so half of
        the available digits are thrown away before the eigenvalue solver is
        even called.  The comparison in the demonstration below makes the
        loss visible.
        """
        w = np.linalg.eigvalsh(self.X.T @ self.X)
        return np.sqrt(np.maximum(w[::-1], 0.0))


# ---------------------------------------------------------------------------
class CanonicalOrthogonalization:
    """Orthogonalise a non-orthogonal basis with a near-singular overlap.

    Given an overlap matrix S with elements S_ij = <phi_i|phi_j>, the
    generalised eigenvalue problem H c = E S c is reduced to an ordinary one
    by the transformation

        X = V s^{-1/2},        S = V s V^T,

    keeping only those eigenvalues of S above a threshold.  Discarding the
    small ones removes the directions in which the basis is very nearly
    linearly dependent, and with them the ill-conditioning.
    """

    def __init__(self, S, threshold=1.0e-8):
        self.S = np.asarray(S, dtype=float)
        self.threshold = threshold
        w, V = np.linalg.eigh(self.S)
        order = np.argsort(w)[::-1]
        self.eigenvalues = w[order]
        self.eigenvectors = V[:, order]
        keep = self.eigenvalues > threshold * self.eigenvalues[0]
        self.kept = int(np.sum(keep))
        self.X = (self.eigenvectors[:, keep]
                  / np.sqrt(self.eigenvalues[keep]))

    @property
    def condition_number(self):
        return self.eigenvalues[0] / self.eigenvalues[-1]

    def transform(self, H):
        """The Hamiltonian in the orthogonalised basis, X^T H X."""
        return self.X.T @ np.asarray(H, dtype=float) @ self.X

    def solve_generalized(self, H):
        """Solve H c = E S c through the orthogonalised problem."""
        values, vectors = np.linalg.eigh(self.transform(H))
        return values, self.X @ vectors


# ---------------------------------------------------------------------------
class SchmidtDecomposition:
    """The SVD of a bipartitioned state vector.

    A state of a system split into parts A and B,

        |Psi> = sum_{ij} C_ij |i>_A |j>_B,

    has the SVD C = U Sigma V^T, hence

        |Psi> = sum_k sigma_k |u_k>_A |v_k>_B,

    which is the Schmidt decomposition.  The singular values are the Schmidt
    coefficients; the eigenvalues of the two reduced density matrices are
    sigma_k^2, and the entanglement entropy is

        S = - sum_k sigma_k^2 log sigma_k^2.

    For a two-particle wave function the same decomposition is the natural
    expansion: the left singular vectors are the natural orbitals and the
    sigma_k^2 are the occupation numbers.
    """

    def __init__(self, C):
        self.C = np.asarray(C, dtype=float)
        self.U, self.sigma, self.Vt = np.linalg.svd(self.C,
                                                    full_matrices=False)
        # a normalised state has sum_k sigma_k^2 = 1
        self.weights = self.sigma**2

    def schmidt_rank(self, tol=1.0e-10):
        """Number of Schmidt coefficients above a relative threshold."""
        return int(np.sum(self.sigma > tol * self.sigma[0]))

    def entanglement_entropy(self, base=np.e):
        w = self.weights[self.weights > 1.0e-16]
        entropy = -np.sum(w * np.log(w)) / np.log(base)
        return float(entropy) if entropy > 0.0 else 0.0

    def states_for(self, tol):
        """How many Schmidt states are needed to capture all but tol."""
        return 1 + int(np.argmax(np.cumsum(self.weights) > 1.0 - tol))

    def truncation_error(self, k):
        """Weight discarded when only k Schmidt states are kept."""
        return float(np.sum(self.weights[k:]))

    def natural_occupations(self):
        """Occupation numbers, i.e. the squared Schmidt coefficients."""
        return self.weights

    def truncate(self, k):
        """The best rank-k approximation to the state, renormalised."""
        C = (self.U[:, :k] * self.sigma[:k]) @ self.Vt[:k]
        return C / np.linalg.norm(C)


# ---------------------------------------------------------------------------
class TwoBodyFactorization:
    """Low-rank factorisation of a two-body interaction.

    The two-body matrix elements, arranged as a matrix in the composite pair
    index, V_{(pq),(rs)}, are positive semidefinite for a Coulomb-like
    interaction, and their eigenvalues decay rapidly.  Keeping only the
    leading ones gives

        V_{(pq),(rs)} = sum_{gamma=0}^{M-1} L^gamma_{pq} L^gamma_{rs},

    with M much smaller than the number of pairs.  This is the Cholesky
    decomposition, or density fitting, of the two-electron integrals, and the
    first half of the double factorisation used to reduce the number of terms
    in a qubit Hamiltonian.
    """

    def __init__(self, V):
        self.V = np.asarray(V, dtype=float)
        w, U = np.linalg.eigh(0.5 * (self.V + self.V.T))
        order = np.argsort(w)[::-1]
        self.eigenvalues = w[order]
        self.vectors = U[:, order]

    def rank_for(self, tol):
        """Smallest M such that the discarded weight is below tol."""
        w = np.maximum(self.eigenvalues, 0.0)
        tail = np.cumsum(w[::-1])[::-1]
        below = np.where(tail <= tol * np.sum(w))[0]
        return int(below[0]) if len(below) else len(w)

    def cholesky_vectors(self, M):
        """The M leading vectors L^gamma."""
        w = np.maximum(self.eigenvalues[:M], 0.0)
        return (self.vectors[:, :M] * np.sqrt(w)).T

    def reconstruct(self, M):
        L = self.cholesky_vectors(M)
        return L.T @ L


# ---------------------------------------------------------------------------
#  Physical models used in the demonstration
# ---------------------------------------------------------------------------
def two_particles_in_a_trap(n=40, rmax=5.0, strength=2.0, softening=0.5):
    """Two particles in a one-dimensional harmonic trap on a grid.

    The Hamiltonian is  h(x1) + h(x2) + g/sqrt((x1-x2)^2 + a^2)  with
    h = -d^2/dx^2 + x^2 discretised by the three-point formula, exactly as in
    the Schrodinger section of the chapter.  We diagonalise in the product
    basis and return the ground state reshaped as the coefficient matrix
    C_{pq}, whose SVD is the Schmidt decomposition.
    """
    h_step = 2.0 * rmax / (n + 1)
    x = -rmax + h_step * np.arange(1, n + 1)

    # one-body Hamiltonian on the grid
    h1 = np.diag(2.0 / h_step**2 + x**2)
    off = -1.0 / h_step**2
    h1 += np.diag(np.full(n - 1, off), 1) + np.diag(np.full(n - 1, off), -1)

    I = np.eye(n)
    H = np.kron(h1, I) + np.kron(I, h1)
    interaction = strength / np.sqrt((x[:, None] - x[None, :])**2 + softening**2)
    H += np.diag(interaction.ravel())

    values, vectors = np.linalg.eigh(H)
    C = vectors[:, 0].reshape(n, n)
    C = 0.5 * (C + C.T)                    # the ground state is symmetric
    return values[0], C / np.linalg.norm(C), x


def pairing_ground_state(levels=12, pairs=6, delta=1.0, g=0.5):
    """Ground state of the pairing model, as a coefficient matrix.

    The single-particle levels are split into two halves, A and B.  A
    seniority-zero basis state is a choice of occupied levels, which
    factorises into its part in A and its part in B, so the coefficients form
    a matrix C indexed by the subsets of A and of B.  Its SVD is the Schmidt
    decomposition across the cut, and is exactly the object that a
    density-matrix renormalisation group calculation truncates.
    """
    from itertools import combinations

    basis = list(combinations(range(1, levels + 1), pairs))
    index = {c: i for i, c in enumerate(basis)}
    dim = len(basis)
    H = np.zeros((dim, dim))
    for c, i in index.items():
        H[i, i] = 2.0 * delta * sum(p - 1 for p in c) - g * pairs
        occupied = set(c)
        for p in c:
            for q in range(1, levels + 1):
                if q not in occupied:
                    H[index[tuple(sorted(occupied - {p} | {q}))], i] -= g

    values, vectors = np.linalg.eigh(H)
    psi = vectors[:, 0]

    half = levels // 2
    left = list(range(1, half + 1))
    subsets_A = [frozenset(s) for k in range(len(left) + 1)
                 for s in combinations(left, k)]
    right = list(range(half + 1, levels + 1))
    subsets_B = [frozenset(s) for k in range(len(right) + 1)
                 for s in combinations(right, k)]
    iA = {s: i for i, s in enumerate(subsets_A)}
    iB = {s: i for i, s in enumerate(subsets_B)}

    C = np.zeros((len(subsets_A), len(subsets_B)))
    for c, i in index.items():
        a = frozenset(p for p in c if p <= half)
        b = frozenset(p for p in c if p > half)
        C[iA[a], iB[b]] = psi[i]
    return values[0], C


def two_body_matrix(n_orb=12, n_grid=200, rmax=4.0, width=1.4, softening=0.5):
    """Two-body matrix elements in a non-orthogonal Gaussian basis.

    The orbitals are Gaussians centred on a uniform set of points, which is a
    caricature of the basis sets used in quantum chemistry -- and, like them,
    close to linearly dependent.  We return both the overlap matrix and the
    two-body matrix in the composite pair index,

        V_{(pq),(rs)} = int dx dy  phi_p(x) phi_q(x) K(x,y) phi_r(y) phi_s(y),

    with a softened Coulomb kernel K.
    """
    h_step = 2.0 * rmax / (n_grid - 1)
    x = np.linspace(-rmax, rmax, n_grid)
    centres = np.linspace(-rmax * 0.8, rmax * 0.8, n_orb)

    phi = np.exp(-(x[None, :] - centres[:, None])**2 / (2.0 * width**2))
    phi /= np.sqrt(h_step * np.sum(phi**2, axis=1))[:, None]

    S = h_step * (phi @ phi.T)                       # overlap matrix
    rho = (phi[:, None, :] * phi[None, :, :]).reshape(n_orb * n_orb, n_grid)
    K = 1.0 / np.sqrt((x[:, None] - x[None, :])**2 + softening**2)
    V = h_step**2 * (rho @ K @ rho.T)                # pair-index supermatrix
    return S, V


# ---------------------------------------------------------------------------
def _demo():
    print("=" * 74)
    print("1. The SVD exists for every matrix, including defective ones")
    print("=" * 74)
    X = np.array([[1.0, -1.0], [1.0, -1.0]])
    svd = SVD(X)
    print("X =\n", X)
    print("singular values     :", np.round(svd.sigma, 12))
    print(f"numerical rank      : {svd.rank()}   (X is singular: det = "
          f"{np.linalg.det(X):.1f})")
    print(f"condition number    : {svd.condition_number()}")
    print("U =\n", np.round(svd.U, 8))
    print("V^T =\n", np.round(svd.Vt, 8))
    print(f"|U S V^T - X|       = "
          f"{np.linalg.norm(svd.reconstruct() - X):.3e}")

    print()
    print("=" * 74)
    print("2. Never form X^T X: the condition number is squared")
    print("=" * 74)
    rng = np.random.default_rng(2026)
    Q, _ = np.linalg.qr(rng.normal(size=(60, 12)))
    P, _ = np.linalg.qr(rng.normal(size=(12, 12)))
    sigma_true = np.logspace(0, -7, 12)
    X = (Q * sigma_true) @ P.T
    svd = SVD(X)
    direct = svd.sigma
    naive = svd.eigen_route()
    print(f"{'exact sigma':>16s} {'from SVD':>16s} {'from X^T X':>16s}"
          f" {'rel. error':>13s}")
    for k in (0, 5, 9, 10, 11):
        print(f"{sigma_true[k]:16.8e} {direct[k]:16.8e} {naive[k]:16.8e}"
              f" {abs(naive[k]-sigma_true[k])/sigma_true[k]:13.2e}")
    print(f"\nkappa(X) = {svd.condition_number():.2e},  "
          f"kappa(X^T X) = {svd.condition_number()**2:.2e}")

    print()
    print("=" * 74)
    print("3. Near-linearly-dependent basis: canonical orthogonalisation")
    print("=" * 74)
    S, V = two_body_matrix()
    co = CanonicalOrthogonalization(S, threshold=0.0)
    print(f"{S.shape[0]} Gaussians, condition number of the overlap matrix: "
          f"{co.condition_number:.3e}")
    print("eigenvalues of S:")
    print(" ", np.array2string(co.eigenvalues, precision=2))
    print()
    for thr in (0.0, 1.0e-12, 1.0e-8, 1.0e-5):
        c = CanonicalOrthogonalization(S, threshold=thr)
        kept_kappa = (c.eigenvalues[0] / c.eigenvalues[c.kept - 1])
        print(f"threshold {thr:8.0e}:  kept {c.kept:2d} of {S.shape[0]}"
              f"   effective condition number {kept_kappa:9.2e}")
    print()
    print("Discarding the directions in which the Gaussians are almost")
    print("linearly dependent removes the ill-conditioning at the cost of a")
    print("slightly smaller space; this is canonical orthogonalisation.")

    print()
    print("=" * 74)
    print("4. Low-rank factorisation of the two-body matrix elements")
    print("=" * 74)
    fact = TwoBodyFactorization(V)
    npairs = V.shape[0]
    print(f"number of pair indices (n_orb^2)       : {npairs}")
    for tol in (1.0e-4, 1.0e-6, 1.0e-8):
        M = fact.rank_for(tol)
        err = np.linalg.norm(fact.reconstruct(M) - V) / np.linalg.norm(V)
        print(f"rank needed for a relative weight {tol:.0e}: {M:3d}"
              f"   |V - L^T L|/|V| = {err:.2e}")
    print("The eigenvalues of the pair matrix decay geometrically:")
    print(" ", np.array2string(fact.eigenvalues[:8], precision=4))

    print()
    print("=" * 74)
    print("5. Schmidt decomposition of two particles in a trap")
    print("=" * 74)
    print(f"{'g':>6s} {'E_0':>12s} {'entropy':>12s} {'n_0':>10s} "
          f"{'n_1':>10s} {'n_2':>11s} {'k for 1e-8':>12s}")
    print("-" * 74)
    for g in (0.0, 1.0, 2.0, 5.0):
        E0, C, x = two_particles_in_a_trap(n=40, strength=g)
        sd = SchmidtDecomposition(C)
        print(f"{g:6.1f} {E0:12.6f} {sd.entanglement_entropy():12.6f} "
              f"{sd.weights[0]:10.6f} {sd.weights[1]:10.6f} "
              f"{sd.weights[2]:11.3e} {sd.states_for(1e-8):12d}")
    print("-" * 74)
    print("At g = 0 the ground state is a product: one Schmidt state, zero")
    print("entropy, a single occupied natural orbital.  The interaction")
    print("generates entanglement, and the n_k are the occupation numbers")
    print("of the natural orbitals -- the same numbers a natural-orbital")
    print("truncation would use to decide which orbitals to keep.")

    print()
    print("=" * 74)
    print("6. Schmidt spectrum of the pairing model across a cut")
    print("=" * 74)
    E0, C = pairing_ground_state()
    sd = SchmidtDecomposition(C)
    print(f"ground-state energy          : {E0:.10f}")
    print(f"coefficient matrix           : {C.shape[0]} x {C.shape[1]}")
    print(f"Schmidt rank (rel. 1e-10)    : {sd.schmidt_rank()}")
    print(f"entanglement entropy         : {sd.entanglement_entropy():.6f}")
    print("largest Schmidt weights      :",
          np.array2string(sd.weights[:7], precision=6))
    print()
    print(f"{'kept states':>12s} {'discarded weight':>20s}")
    for k in (1, 2, 3, 5, 7, 10, 15):
        print(f"{k:12d} {sd.truncation_error(k):20.3e}")
    print()
    print(f"Keeping {sd.states_for(1e-6):d} of the {C.shape[0]:d} Schmidt "
          f"states already reproduces the state")
    print("to one part in a million.  This decay is exactly what a DMRG")
    print("calculation exploits when it truncates the basis at each cut.")


if __name__ == "__main__":
    _demo()
