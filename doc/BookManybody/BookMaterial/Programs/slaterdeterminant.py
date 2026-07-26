"""
Slater determinants, basis changes and the energy functional.

Companion code to chapter 2 of *Quantum mechanics for Many-particle Systems*.

Everything in chapter 2 is linear algebra applied to antisymmetric wave
functions, and every statement made there can be checked numerically:

    HarmonicOscillatorBasis -- a single-particle basis with one- and two-body
                               matrix elements
    SlaterDeterminant       -- the determinant itself: antisymmetry, the Pauli
                               principle, and the det(C) det(Phi) rule
    EnergyFunctional        -- E[Phi] and E[Psi] in terms of the coefficients
                               C, and the invariance of the energy under
                               unitary rotations within the occupied space
    HartreeFock             -- a minimal self-consistent field loop, included
                               to show that the Hartree-Fock equations really
                               are a repeated eigenvalue problem

The model is a system of spinless fermions in a one-dimensional harmonic trap
interacting through a softened Coulomb repulsion.  Spin is suppressed
throughout: it adds bookkeeping but nothing conceptual.

Author: Morten Hjorth-Jensen
"""

import numpy as np

from householder import SymmetricEigenSolver


# ---------------------------------------------------------------------------
class HarmonicOscillatorBasis:
    """Single-particle basis of harmonic oscillator eigenfunctions.

    The orbitals are the Hermite functions

        phi_n(x) = (2^n n! sqrt(pi))^{-1/2} H_n(x) exp(-x^2/2),

    which satisfy h_0 phi_n = (n + 1/2) phi_n with h_0 = -(1/2) d^2/dx^2 +
    x^2/2.  Because they are eigenfunctions of the one-body Hamiltonian, the
    one-body matrix elements are diagonal, <p|h_0|q> = eps_p delta_pq, which
    is exactly the situation described in the text.
    """

    def __init__(self, n_orb=8, n_grid=601, rmax=8.0, strength=1.0,
                 softening=0.5):
        self.n_orb = n_orb
        self.strength = strength
        self.softening = softening

        self.x, self.h = np.linspace(-rmax, rmax, n_grid, retstep=True)
        self.phi = self._hermite_functions()
        self.epsilon = np.arange(n_orb) + 0.5

    # ------------------------------------------------------------------
    def _hermite_functions(self):
        """The first n_orb Hermite functions, normalised on the grid."""
        x = self.x
        phi = np.zeros((self.n_orb, len(x)))
        gauss = np.exp(-x**2 / 2.0)
        # stable upward recursion  H_{n+1} = 2x H_n - 2n H_{n-1}
        H_prev, H = np.zeros_like(x), np.ones_like(x)
        for n in range(self.n_orb):
            phi[n] = H * gauss
            H_prev, H = H, 2.0 * x * H - 2.0 * n * H_prev
        # normalise numerically rather than analytically: the grid is finite
        phi /= np.sqrt(self.h * np.sum(phi**2, axis=1))[:, None]
        return phi

    # ------------------------------------------------------------------
    @property
    def overlap(self):
        """S_pq = <phi_p|phi_q>, the identity up to quadrature error."""
        return self.h * (self.phi @ self.phi.T)

    @property
    def one_body(self):
        """h_pq = <p|h_0|q> = eps_p delta_pq in this basis."""
        return np.diag(self.epsilon)

    def two_body(self):
        """v[p,q,r,s] = <pq|v|rs>, the direct (non-antisymmetrised) elements.

        With a local interaction v(x,y) the integral factorises into pair
        densities, so the whole four-index array is one matrix product --
        which is the pair-index picture of Section 1.32.
        """
        n, x, h = self.n_orb, self.x, self.h
        K = self.strength / np.sqrt((x[:, None] - x[None, :])**2
                                    + self.softening**2)
        # pair densities rho_{pr}(x) = phi_p(x) phi_r(x)
        rho = (self.phi[:, None, :] * self.phi[None, :, :]).reshape(n * n, -1)
        V = h**2 * (rho @ K @ rho.T)              # V[(pr), (qs)]
        return V.reshape(n, n, n, n).transpose(0, 2, 1, 3)

    @staticmethod
    def antisymmetrize(v):
        """<pq|v|rs>_AS = <pq|v|rs> - <pq|v|sr>."""
        return v - v.transpose(0, 1, 3, 2)


# ---------------------------------------------------------------------------
class SlaterDeterminant:
    """An antisymmetric N-particle wave function built from orbitals.

    The determinant is evaluated directly from its definition,

        Phi(x_1..x_N) = (1/sqrt(N!)) det[ psi_a(x_b) ],

    which is all that is needed to verify antisymmetry and the Pauli
    principle numerically.
    """

    def __init__(self, basis, occupied, coefficients=None):
        self.basis = basis
        self.occupied = list(occupied)
        self.N = len(self.occupied)
        # coefficients[i, lambda]: the new orbital i in terms of the old basis
        if coefficients is None:
            C = np.zeros((self.N, basis.n_orb))
            for i, p in enumerate(self.occupied):
                C[i, p] = 1.0
            self.C = C
        else:
            self.C = np.asarray(coefficients, dtype=float)

    # ------------------------------------------------------------------
    def orbitals_on_grid(self):
        """The occupied orbitals psi_i(x) = sum_lambda C_{i lambda} phi_lambda."""
        return self.C @ self.basis.phi

    def evaluate(self, coordinates):
        """Phi(x_1, ..., x_N) at one set of coordinates, by interpolation."""
        psi = self.orbitals_on_grid()
        M = np.array([[np.interp(x, self.basis.x, psi[i])
                       for x in coordinates] for i in range(self.N)])
        return np.linalg.det(M) / np.sqrt(_factorial(self.N))

    def density_matrix(self):
        """rho_{lambda mu} = sum_i C*_{i lambda} C_{i mu}.

        For a Slater determinant this matrix is a projector: its eigenvalues,
        the occupation numbers, are exactly 1 (N times) and 0.  Compare with
        the natural occupations of a correlated state in Section 1.31.
        """
        return self.C.T @ self.C

    def occupation_numbers(self):
        return np.sort(np.linalg.eigvalsh(self.density_matrix()))[::-1]


def _factorial(n):
    out = 1
    for k in range(2, n + 1):
        out *= k
    return out


# ---------------------------------------------------------------------------
class EnergyFunctional:
    """The Slater-determinant expectation value of the Hamiltonian.

        E[Phi] = sum_i <i|h_0|i> + (1/2) sum_ij <ij|v|ij>_AS ,

    and, written in terms of the expansion coefficients,

        E[Psi] = sum_i sum_{ab} C*_{ia} C_{ib} <a|h_0|b>
               + (1/2) sum_ij sum_{abcd} C*_{ia} C*_{jb} C_{ic} C_{jd}
                       <ab|v|cd>_AS .

    The second form is the one that is varied to obtain the Hartree-Fock
    equations.
    """

    def __init__(self, basis):
        self.basis = basis
        self.h = basis.one_body
        self.v = basis.antisymmetrize(basis.two_body())

    # ------------------------------------------------------------------
    def energy(self, C):
        """E[Psi] for occupied orbitals given by the rows of C."""
        C = np.asarray(C, dtype=float)
        one = np.einsum("ia,ib,ab->", C, C, self.h)
        two = np.einsum("ia,jb,ic,jd,abcd->", C, C, C, C, self.v)
        return one + 0.5 * two

    def energy_from_occupied(self, occupied):
        """E[Phi] when the occupied orbitals are basis states."""
        occ = list(occupied)
        one = sum(self.h[i, i] for i in occ)
        two = sum(self.v[i, j, i, j] for i in occ for j in occ)
        return one + 0.5 * two

    def fock_matrix(self, C):
        """f_{ab} = h_{ab} + sum_i sum_{cd} C*_{ic} C_{id} <ac|v|bd>_AS."""
        rho = C.T @ C                       # the density matrix
        return self.h + np.einsum("cd,acbd->ab", rho, self.v)


# ---------------------------------------------------------------------------
class HartreeFock:
    """A minimal self-consistent field loop.

    The Hartree-Fock equations say that the optimal coefficients are the
    eigenvectors of the Fock matrix -- which itself depends on them.  So one
    diagonalises, rebuilds, and repeats.  The diagonalisation uses the
    Householder plus QL solver of chapter 1; nothing else is needed.

    This is a preview: Hartree-Fock theory proper is developed later in the
    book.
    """

    def __init__(self, functional, n_particles, max_iter=200, tol=1.0e-10):
        self.functional = functional
        self.N = n_particles
        self.max_iter = max_iter
        self.tol = tol

    def solve(self, C0=None):
        n = self.functional.basis.n_orb
        C = np.eye(n)[:self.N] if C0 is None else np.asarray(C0, dtype=float)
        energy_old = np.inf
        self.iterations = 0
        self.history = []

        for k in range(self.max_iter):
            f = self.functional.fock_matrix(C)
            values, vectors = SymmetricEigenSolver(f).eigen()
            C = vectors[:, :self.N].T             # lowest N eigenvectors
            energy = self.functional.energy(C)
            self.history.append(energy)
            self.iterations = k + 1
            if abs(energy - energy_old) < self.tol:
                break
            energy_old = energy
        self.C, self.epsilon = C, values
        return energy, C


# ---------------------------------------------------------------------------
def exact_two_particle(functional):
    """Exact energy of two spinless fermions in the truncated basis.

    The antisymmetric two-particle states |pq> with p < q form a basis of
    dimension n(n-1)/2, and the Hamiltonian in that basis is small enough to
    diagonalise directly.  The difference from the Slater-determinant result
    is the correlation energy.
    """
    n = functional.basis.n_orb
    pairs = [(p, q) for p in range(n) for q in range(p + 1, n)]
    dim = len(pairs)
    H = np.zeros((dim, dim))
    h, v = functional.h, functional.v

    for a, (p, q) in enumerate(pairs):
        for b, (r, s) in enumerate(pairs):
            element = v[p, q, r, s]
            if q == s:
                element += h[p, r]
            if p == r:
                element += h[q, s]
            if q == r:
                element -= h[p, s]
            if p == s:
                element -= h[q, r]
            H[a, b] = element
    values, vectors = np.linalg.eigh(H)
    return values[0], vectors[:, 0], pairs


# ---------------------------------------------------------------------------
def _demo():
    rng = np.random.default_rng(2026)
    basis = HarmonicOscillatorBasis(n_orb=8, strength=1.0)

    print("=" * 74)
    print("1. The single-particle basis")
    print("=" * 74)
    print(f"{basis.n_orb} harmonic oscillator orbitals on a grid of "
          f"{len(basis.x)} points")
    print(f"|S - I|                  = "
          f"{np.linalg.norm(basis.overlap - np.eye(basis.n_orb)):.2e}"
          f"   (orthonormal to quadrature accuracy)")
    print(f"single-particle energies = {basis.epsilon}")

    print()
    print("=" * 74)
    print("2. The Slater determinant is antisymmetric and obeys Pauli")
    print("=" * 74)
    phi = SlaterDeterminant(basis, [0, 1, 2])
    coords = [0.4, -1.1, 2.0]
    swapped = [coords[1], coords[0], coords[2]]
    a, b = phi.evaluate(coords), phi.evaluate(swapped)
    print(f"Phi(x1,x2,x3)            = {a:+.10f}")
    print(f"Phi(x2,x1,x3)            = {b:+.10f}")
    print(f"sum                      = {a + b:+.2e}   "
          f"(antisymmetry under exchange)")

    repeated = SlaterDeterminant(basis, [0, 1, 1])
    print(f"two particles in orbital 1: Phi = "
          f"{repeated.evaluate(coords):+.2e}   (Pauli principle)")

    print()
    print("=" * 74)
    print("3. Occupation numbers of a Slater determinant are 1 and 0")
    print("=" * 74)
    print("occupations:", np.round(SlaterDeterminant(basis,
                                                    [0, 1, 2]).occupation_numbers(), 12))
    print("Compare with the natural occupations of a correlated state, which")
    print("lie strictly between 0 and 1 (Section 1.31).")

    print()
    print("=" * 74)
    print("4. Basis change: det(C) det(Phi), and the energy is invariant")
    print("=" * 74)
    functional = EnergyFunctional(basis)
    N = 3
    C0 = np.eye(basis.n_orb)[:N]
    E0 = functional.energy(C0)

    # a random unitary rotation *within* the occupied space
    M = rng.normal(size=(N, N))
    Q, _ = np.linalg.qr(M)
    C1 = Q @ C0
    E1 = functional.energy(C1)
    print(f"det(Q)                   = {np.linalg.det(Q):+.10f}"
          f"   (|det| = {abs(np.linalg.det(Q)):.10f})")
    print(f"E[Phi]  original basis   = {E0:.12f}")
    print(f"E[Phi]  rotated basis    = {E1:.12f}")
    print(f"difference               = {abs(E1 - E0):.2e}")
    print("The Slater determinant changes only by the phase det(Q), so the")
    print("energy cannot change: this is the invariance used in Hartree-Fock.")

    # a rotation that mixes occupied and virtual orbitals does change it
    M = rng.normal(size=(basis.n_orb, basis.n_orb))
    Qfull, _ = np.linalg.qr(M)
    E2 = functional.energy(Qfull[:, :N].T)
    print(f"E[Phi]  occupied-virtual mixing = {E2:.12f}   -- and this one does"
          f" change")

    print()
    print("=" * 74)
    print("5. The energy functional and the cost of the two-body elements")
    print("=" * 74)
    print(f"E[Phi] from the sum over occupied states = "
          f"{functional.energy_from_occupied(range(N)):.12f}")
    print(f"E[Psi] from the coefficient formula      = {E0:.12f}")
    n = basis.n_orb
    print(f"\ntwo-body matrix elements stored: n^4 = {n**4}")
    print(f"pair-index matrix is {n*n} x {n*n}; its numerical rank is",
          np.linalg.matrix_rank(basis.two_body().transpose(0, 2, 1, 3)
                                .reshape(n*n, n*n), tol=1e-10))
    print("The rank is far below n^2, which is what the low-rank")
    print("factorisation of Section 1.32 exploits.")

    print()
    print("=" * 74)
    print("6. Hartree-Fock as a repeated eigenvalue problem")
    print("=" * 74)
    for N in (2, 3, 4):
        hf = HartreeFock(functional, N)
        E, C = hf.solve()
        E_ref = functional.energy_from_occupied(range(N))
        print(f"N = {N}:  E(lowest N orbitals) = {E_ref:12.8f}   "
              f"E(Hartree-Fock) = {E:12.8f}   "
              f"gain = {E - E_ref:+.2e}   ({hf.iterations} iterations)")

    print()
    print("=" * 74)
    print("7. Correlation: two particles, Hartree-Fock against exact")
    print("=" * 74)
    hf = HartreeFock(functional, 2)
    E_hf, _ = hf.solve()
    E_exact, psi, pairs = exact_two_particle(functional)
    print(f"E(Hartree-Fock)          = {E_hf:.10f}")
    print(f"E(exact, same basis)     = {E_exact:.10f}")
    print(f"correlation energy       = {E_exact - E_hf:+.10f}")
    print(f"largest amplitude        = {np.max(np.abs(psi)):.6f} "
          f"for the pair {pairs[int(np.argmax(np.abs(psi)))]}")
    print("\nThe exact state is not a single Slater determinant; the missing")
    print("energy is precisely what the methods of the following chapters")
    print("are built to recover.")


if __name__ == "__main__":
    _demo()
