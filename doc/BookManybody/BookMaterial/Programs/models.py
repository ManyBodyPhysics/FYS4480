"""
The model Hamiltonians of chapter 4.

Companion code to chapter 4 of *Quantum mechanics for Many-particle Systems*.

Each class builds the Hamiltonian of one schematic model and diagonalises it,
so that every number quoted in the chapter can be reproduced:

    LipkinModel      -- two levels, degeneracy Omega, quasispin SU(2)
    PairingModel     -- L doubly degenerate levels, pair-transfer interaction
    PairingPHModel   -- the same, plus a pair-breaking particle-hole term
    HubbardChain     -- N sites, hopping t and on-site repulsion U
    HeisenbergChain  -- N spins on a ring, exchange coupling J
    CalogeroModel    -- N particles on a line, inverse-square interaction

Most of them are diagonalised in a fixed symmetry sector; the Calogero model
is exactly solvable and we check the analytic ground-state energy against a
numerical solution of the two-body problem.

Two of the classes exist mainly to be compared with others.  PairingPHModel
reduces to PairingModel at f = 0, and reproduces the exact energies of
`rpa.py`, where the same model carries the whole discussion of the
Tamm-Dancoff and random-phase approximations.  HeisenbergChain is the strong-
coupling limit of HubbardChain, and the two are compared directly in the demo.

Author: Morten Hjorth-Jensen
"""

from itertools import combinations

import numpy as np


# ---------------------------------------------------------------------------
class LipkinModel:
    """The Lipkin-Meshkov-Glick model in the quasispin basis.

    Two levels with degeneracy Omega, quantum numbers sigma = +/- 1 and
    p = 1, ..., Omega.  In terms of the quasispin operators

        J_z = (1/2) sum_{p sigma} sigma a^+_{p sigma} a_{p sigma},
        J_+ = sum_p a^+_{p+} a_{p-},          J_- = (J_+)^+,

    the Hamiltonian becomes

        H = eps J_z + (V/2)(J_+^2 + J_-^2) + (W/2)(-N + J_+ J_- + J_- J_+).

    Since [H, J^2] = 0 the Hamiltonian is block diagonal in J, and each block
    is a small matrix in the basis |J, J_z> with J_z = -J, ..., J.  For N
    particles at half filling the maximal block is J = N/2.
    """

    def __init__(self, J=2.0, eps=2.0, V=-1.0/3.0, W=-1.0/4.0, N=None):
        self.J = J
        self.eps = eps
        self.V = V
        self.W = W
        self.N = 2 * J if N is None else N
        self.jz = np.arange(-J, J + 1)           # the 2J+1 projections
        self.dim = len(self.jz)

    # ------------------------------------------------------------------
    @staticmethod
    def _jplus(J, m):
        """<J, m+1| J_+ |J, m> = sqrt(J(J+1) - m(m+1))."""
        value = J * (J + 1) - m * (m + 1)
        return np.sqrt(value) if value > 0 else 0.0

    def matrix(self):
        """The Hamiltonian in the |J, J_z> basis of the chosen J block."""
        J, dim = self.J, self.dim
        H = np.zeros((dim, dim))

        for a, m in enumerate(self.jz):
            H[a, a] += self.eps * m                       # eps J_z

            # (W/2)(-N + J_+ J_- + J_- J_+) is diagonal:
            #   J_+ J_- + J_- J_+ = 2(J(J+1) - J_z^2)
            H[a, a] += 0.5 * self.W * (-self.N
                                       + 2.0 * (J * (J + 1) - m * m))

            # (V/2)(J_+^2 + J_-^2) connects J_z to J_z +/- 2
            if a + 2 < dim:
                element = 0.5 * self.V * self._jplus(J, m) \
                    * self._jplus(J, m + 1)
                H[a + 2, a] += element
                H[a, a + 2] += element
        return H

    def spectrum(self):
        return np.linalg.eigvalsh(self.matrix())

    def ground_state(self):
        values, vectors = np.linalg.eigh(self.matrix())
        return values[0], vectors[:, 0]


# ---------------------------------------------------------------------------
class PairingModel:
    """The pairing model in the space with no broken pairs.

        H = xi sum_{p sigma} (p-1) a^+_{p sigma} a_{p sigma}
          - (g/2) sum_{pq} P^+_p P^-_q ,
        P^+_p = a^+_{p+} a^+_{p-} .

    With total spin S = 0 and no broken pairs, a basis state is a choice of
    which levels carry a pair, so the dimension is binom(L, n).  The
    interaction moves one pair from level q to level p with amplitude -g/2,
    and the p = q terms give a diagonal -g/2 per occupied level.
    """

    def __init__(self, levels=4, pairs=2, g=1.0, xi=1.0):
        self.levels = levels
        self.pairs = pairs
        self.g = g
        self.xi = xi
        self.basis = list(combinations(range(1, levels + 1), pairs))
        self.index = {c: i for i, c in enumerate(self.basis)}
        self.dim = len(self.basis)

    def matrix(self):
        H = np.zeros((self.dim, self.dim))
        for config, i in self.index.items():
            occupied = set(config)
            # unperturbed: 2 xi (p-1) per occupied level (two particles)
            H[i, i] += 2.0 * self.xi * sum(p - 1 for p in config)
            # p = q terms of the pairing interaction
            H[i, i] -= 0.5 * self.g * self.pairs
            # p /= q: move one pair
            for q in config:
                for p in range(1, self.levels + 1):
                    if p in occupied:
                        continue
                    new = tuple(sorted(occupied - {q} | {p}))
                    H[self.index[new], i] -= 0.5 * self.g
        return H

    def spectrum(self):
        return np.linalg.eigvalsh(self.matrix())

    def ground_state_energy(self):
        return float(np.linalg.eigvalsh(self.matrix())[0])

    def reference_energy(self):
        """The energy of the unperturbed ground-state configuration."""
        lowest = tuple(range(1, self.pairs + 1))
        i = self.index[lowest]
        return float(self.matrix()[i, i])


# ---------------------------------------------------------------------------
#  The pairing model with a particle-hole term
# ---------------------------------------------------------------------------
UP, DN = 0, 1


def orbital(level, spin):
    """(level = 1..L, spin) -> spin-orbital index 0..2L-1."""
    return 2 * (level - 1) + spin


def level_of(o):
    return o // 2 + 1


def _annihilate(state, o):
    if not (state >> o) & 1:
        return None, 0
    sign = -1 if bin(state & ((1 << o) - 1)).count("1") & 1 else 1
    return state & ~(1 << o), sign


def _create(state, o):
    if (state >> o) & 1:
        return None, 0
    sign = -1 if bin(state & ((1 << o) - 1)).count("1") & 1 else 1
    return state | (1 << o), sign


def apply_2body(state, i, j, k, l):
    """a^+_i a^+_j a_k a_l |state>, returning (new state, sign)."""
    s, sgn = _annihilate(state, l)
    if s is None:
        return None, 0
    s, t = _annihilate(s, k)
    if s is None:
        return None, 0
    sgn *= t
    s, t = _create(s, j)
    if s is None:
        return None, 0
    sgn *= t
    s, t = _create(s, i)
    if s is None:
        return None, 0
    return s, sgn * t


class PairingPHModel:
    """The pairing model with a pair-breaking particle-hole term,

        H = xi sum_{p sigma} (p-1) a^+_{p sigma} a_{p sigma}
          - (g/2) sum_{pq}  a^+_{p+} a^+_{p-} a_{q-} a_{q+}
          - (f/2) sum_{pqr} (a^+_{p+} a^+_{p-} a_{q-} a_{r+} + h.c.) .

    Because the interaction no longer conserves the seniority, the basis can no
    longer be the set of pair configurations: we work with every determinant of
    N particles in 2L spin-orbitals, restricted to a fixed 2 S_z = N_up - N_dn.
    At f = 0 the model is PairingModel above, and this class reproduces its
    energies exactly; it is also the model of chapter 7, where the same
    Hamiltonian is treated in the Tamm-Dancoff and random-phase approximations
    by `rpa.py`.
    """

    def __init__(self, levels=4, particles=4, g=1.0, f=0.0, xi=1.0, sz2=0):
        self.L = levels
        self.N = particles
        self.g = g
        self.f = f
        self.xi = xi
        self.sz2 = sz2
        self.states = self._basis()
        self.index = {s: i for i, s in enumerate(self.states)}
        self.dim = len(self.states)

    def _basis(self):
        out = []
        for occ in combinations(range(2 * self.L), self.N):
            if self.sz2 is not None:
                n_up = sum(1 for o in occ if o % 2 == UP)
                if n_up - (self.N - n_up) != self.sz2:
                    continue
            bits = 0
            for o in occ:
                bits |= 1 << o
            out.append(bits)
        return sorted(out)

    def matrix(self):
        L, dim = self.L, self.dim
        H = np.zeros((dim, dim))
        for col, state in enumerate(self.states):
            H[col, col] += sum(self.xi * (level_of(o) - 1)
                               for o in range(2 * L) if (state >> o) & 1)
            # pairing: -(g/2) sum_pq a^+_{p+} a^+_{p-} a_{q-} a_{q+}
            for p in range(1, L + 1):
                i, j = orbital(p, UP), orbital(p, DN)
                for q in range(1, L + 1):
                    k, l = orbital(q, DN), orbital(q, UP)
                    s, sgn = apply_2body(state, i, j, k, l)
                    if s is not None and s in self.index:
                        H[self.index[s], col] += -0.5 * self.g * sgn
            # particle-hole: -(f/2) sum_pqr a^+_{p+} a^+_{p-} a_{q-} a_{r+}
            for p in range(1, L + 1):
                i, j = orbital(p, UP), orbital(p, DN)
                for q in range(1, L + 1):
                    k = orbital(q, DN)
                    for r in range(1, L + 1):
                        l = orbital(r, UP)
                        s, sgn = apply_2body(state, i, j, k, l)
                        if s is not None and s in self.index:
                            H[self.index[s], col] += -0.5 * self.f * sgn
            # and its Hermitian conjugate
            for p in range(1, L + 1):
                k, l = orbital(p, DN), orbital(p, UP)
                for q in range(1, L + 1):
                    j = orbital(q, DN)
                    for r in range(1, L + 1):
                        i = orbital(r, UP)
                        s, sgn = apply_2body(state, i, j, k, l)
                        if s is not None and s in self.index:
                            H[self.index[s], col] += -0.5 * self.f * sgn
        return H

    def spectrum(self, k=6):
        return np.sort(np.linalg.eigvalsh(self.matrix()))[:k]

    def ground_state_energy(self):
        return float(np.linalg.eigvalsh(self.matrix())[0])

    def reference(self):
        """The determinant filling the N/2 lowest levels completely."""
        bits = 0
        for p in range(1, self.N // 2 + 1):
            bits |= 1 << orbital(p, UP)
            bits |= 1 << orbital(p, DN)
        return bits

    def reference_energy(self):
        return float(self.matrix()[self.index[self.reference()],
                                   self.index[self.reference()]])

    # -- structure of the excitations ------------------------------------
    @staticmethod
    def seniority(state, levels):
        """The number of particles that are not part of a complete pair."""
        unpaired = 0
        for p in range(1, levels + 1):
            up = (state >> orbital(p, UP)) & 1
            dn = (state >> orbital(p, DN)) & 1
            unpaired += up ^ dn
        return unpaired

    def excitation_rank(self, state):
        """How many spin-orbitals distinguish `state` from the reference.

        Returns the number of particle-hole pairs: 1 for a 1p-1h determinant,
        2 for a 2p-2h determinant, and so on.
        """
        return bin(state ^ self.reference()).count("1") // 2

    def couplings_to_reference(self):
        """<excited| H |reference>, sorted by excitation rank and seniority.

        This is the table that says which terms of the Hamiltonian reach which
        determinants.  With f = 0 only the seniority-zero 2p-2h states are
        connected to the reference; the particle-hole term opens up the 1p-1h
        channel and the 2p-2h states of seniority four.
        """
        H = self.matrix()
        column = H[:, self.index[self.reference()]]
        out = {}
        for state, i in self.index.items():
            if state == self.reference() or abs(column[i]) < 1e-14:
                continue
            key = (self.excitation_rank(state),
                   self.seniority(state, self.L))
            out.setdefault(key, []).append(float(column[i]))
        return {k: out[k] for k in sorted(out)}


# ---------------------------------------------------------------------------
class HubbardChain:
    """The one-dimensional Fermi-Hubbard model,

        H = -t sum_{<ij> sigma} (c^+_{i sigma} c_{j sigma} + h.c.)
            + U sum_i n_{i up} n_{i down} ,

    diagonalised in a fixed (N_up, N_down) sector.  Because the Hamiltonian
    conserves the two particle numbers separately, the sector factorises as
    up (x) down for the hopping part, and the interaction is diagonal in the
    occupation basis.  This is the structure exploited in the accompanying
    notebook `Hubbard.ipynb`, here reduced to its simplest form.
    """

    def __init__(self, sites=4, n_up=2, n_down=2, t=1.0, U=4.0, pbc=True):
        self.n = sites
        self.n_up = n_up
        self.n_down = n_down
        self.t = t
        self.U = U
        self.pbc = pbc
        self.bonds = [(i, (i + 1) % sites) for i in range(sites)]
        if not pbc:
            self.bonds = self.bonds[:-1]
        self.states_up = self._states(n_up)
        self.states_down = self._states(n_down)
        self.dim = len(self.states_up) * len(self.states_down)

    def _states(self, n_particles):
        """All occupation bit patterns with n_particles particles."""
        out = []
        for occ in combinations(range(self.n), n_particles):
            bits = 0
            for site in occ:
                bits |= 1 << site
            out.append(bits)
        return sorted(out)

    def _hopping(self, states):
        """The single-species hopping matrix, with Jordan-Wigner signs."""
        index = {s: i for i, s in enumerate(states)}
        K = np.zeros((len(states), len(states)))
        for s in states:
            for (i, j) in self.bonds:
                for (c, d) in ((i, j), (j, i)):      # c^+_c c_d
                    if not (s >> d) & 1:
                        continue
                    if (s >> c) & 1:
                        continue
                    # sign: parity of the occupied sites strictly between
                    lo, hi = min(c, d), max(c, d)
                    mask = ((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1)
                    sign = (-1) ** bin(s & mask).count("1")
                    new = (s ^ (1 << d)) | (1 << c)
                    K[index[new], index[s]] += -self.t * sign
        return K

    def matrix(self):
        Ku = self._hopping(self.states_up)
        Kd = self._hopping(self.states_down)
        Iu = np.eye(len(self.states_up))
        Id = np.eye(len(self.states_down))
        H = np.kron(Ku, Id) + np.kron(Iu, Kd)

        double = np.array([[bin(su & sd).count("1")
                            for sd in self.states_down]
                           for su in self.states_up]).ravel()
        H += self.U * np.diag(double.astype(float))
        return H

    def spectrum(self):
        return np.linalg.eigvalsh(self.matrix())

    def ground_state_energy(self):
        return float(np.linalg.eigvalsh(self.matrix())[0])

    def double_occupancy(self):
        """<n_up n_down> summed over sites, in the ground state."""
        values, vectors = np.linalg.eigh(self.matrix())
        psi = vectors[:, 0]
        double = np.array([[bin(su & sd).count("1")
                            for sd in self.states_down]
                           for su in self.states_up]).ravel()
        return float(psi @ (double * psi))


# ---------------------------------------------------------------------------
class HeisenbergChain:
    """The spin-1/2 quantum Heisenberg model on a ring,

        H = J sum_{<ij>} S_i . S_j - h sum_i S^z_i
          = J sum_{<ij>} [ S^z_i S^z_j + (1/2)(S^+_i S^-_j + S^-_i S^+_j) ]
            - h sum_i S^z_i .

    A basis state is a bit pattern, bit i set meaning spin up on site i, and
    the total magnetisation sum_i S^z_i is conserved, so the matrix is block
    diagonal in it.  Spin operators on different sites commute, so -- unlike
    every fermionic model in this file -- there are no Jordan-Wigner signs to
    keep track of.

    For J > 0 the model is antiferromagnetic, and it is the strong-coupling
    limit of the half-filled Hubbard model with J = 4 t^2 / U.
    """

    def __init__(self, sites=4, J=1.0, h=0.0, pbc=True, magnetisation=0):
        self.n = sites
        self.J = J
        self.h = h
        self.pbc = pbc
        self.magnetisation = magnetisation      # 2 S^z_total, or None
        bonds = [(i, (i + 1) % sites) for i in range(sites)]
        if not pbc:
            bonds = bonds[:-1]
        # a two-site ring has one bond, not two
        seen, self.bonds = set(), []
        for i, j in bonds:
            if frozenset((i, j)) in seen or i == j:
                continue
            seen.add(frozenset((i, j)))
            self.bonds.append((i, j))
        self.states = self._basis()
        self.index = {s: i for i, s in enumerate(self.states)}
        self.dim = len(self.states)

    def _basis(self):
        out = []
        for state in range(1 << self.n):
            up = bin(state).count("1")
            if self.magnetisation is not None and \
                    2 * up - self.n != self.magnetisation:
                continue
            out.append(state)
        return out

    def matrix(self):
        H = np.zeros((self.dim, self.dim))
        for state, col in self.index.items():
            # -h sum_i S^z_i
            H[col, col] -= self.h * 0.5 * (2 * bin(state).count("1") - self.n)
            for i, j in self.bonds:
                si = (state >> i) & 1
                sj = (state >> j) & 1
                # S^z_i S^z_j
                H[col, col] += self.J * (0.25 if si == sj else -0.25)
                # the flip-flop term connects antiparallel neighbours only
                if si != sj:
                    flipped = state ^ (1 << i) ^ (1 << j)
                    H[self.index[flipped], col] += 0.5 * self.J
        return H

    def sparse_matrix(self):
        """The same Hamiltonian in compressed sparse row form.

        Each row has at most one non-zero entry per bond, so the matrix is
        extremely sparse and Lanczos is the natural solver once the chain gets
        beyond about twelve sites.
        """
        from scipy.sparse import coo_matrix
        rows, cols, data = [], [], []
        for state, col in self.index.items():
            diagonal = -self.h * 0.5 * (2 * bin(state).count("1") - self.n)
            for i, j in self.bonds:
                si = (state >> i) & 1
                sj = (state >> j) & 1
                diagonal += self.J * (0.25 if si == sj else -0.25)
                if si != sj:
                    flipped = state ^ (1 << i) ^ (1 << j)
                    rows.append(self.index[flipped])
                    cols.append(col)
                    data.append(0.5 * self.J)
            rows.append(col)
            cols.append(col)
            data.append(diagonal)
        return coo_matrix((data, (rows, cols)),
                          shape=(self.dim, self.dim)).tocsr()

    def spectrum(self, k=None):
        values = np.linalg.eigvalsh(self.matrix())
        return values if k is None else values[:k]

    def ground_state_energy(self, sparse=None):
        """The lowest eigenvalue, by dense diagonalisation or by Lanczos."""
        if sparse is None:
            sparse = self.dim > 1500
        if not sparse:
            return float(np.linalg.eigvalsh(self.matrix())[0])
        from scipy.sparse.linalg import eigsh
        values = eigsh(self.sparse_matrix(), k=1, which="SA",
                       return_eigenvectors=False, tol=0.0)
        return float(values[0])

    def energy_per_site(self, sparse=None):
        return self.ground_state_energy(sparse) / self.n

    def ground_state(self):
        values, vectors = np.linalg.eigh(self.matrix())
        return vectors[:, 0]

    def correlation(self, i, j):
        """<S_i . S_j> in the ground state."""
        if i == j:
            return 0.75                     # S_i . S_i = s(s+1) = 3/4
        psi = self.ground_state()
        total = 0.0
        for state, col in self.index.items():
            si = (state >> i) & 1
            sj = (state >> j) & 1
            total += psi[col] * psi[col] * (0.25 if si == sj else -0.25)
            if si != sj:
                flipped = state ^ (1 << i) ^ (1 << j)
                total += 0.5 * psi[self.index[flipped]] * psi[col]
        return float(total)

    def staggered_structure_factor(self):
        """(1/N) sum_{ij} (-1)^(i-j) <S_i . S_j>, the antiferromagnetic order."""
        total = 0.0
        for i in range(self.n):
            for j in range(self.n):
                total += (-1) ** (i - j) * self.correlation(i, j)
        return total / self.n


def bethe_energy_per_site():
    """The Bethe-ansatz result for the infinite spin-1/2 antiferromagnetic
    chain with H = J sum_i S_i . S_{i+1},  e_0 = J (1/4 - ln 2)."""
    return 0.25 - np.log(2.0)


# ---------------------------------------------------------------------------
class CalogeroModel:
    """The Calogero model: N particles on a line with inverse-square coupling,

        H = -(1/2) sum_i d^2/dx_i^2 + (1/2) omega^2 sum_i x_i^2
            + sum_{i<j} lambda(lambda-1) / (x_i - x_j)^2 ,

    in units hbar = m = 1.  The model is exactly solvable: the ground state is
    a Jastrow function times a Gaussian,

        Psi_0 = prod_{i<j} |x_i - x_j|^lambda exp(-omega sum_i x_i^2 / 2),

    with energy

        E_0 = omega [ N/2 + lambda N(N-1)/2 ] .

    For N = 2 the relative motion separates and can be solved on a grid, which
    is the check performed below.
    """

    def __init__(self, N=2, lam=2.0, omega=1.0):
        self.N = N
        self.lam = lam
        self.omega = omega

    # ------------------------------------------------------------------
    def exact_energy(self):
        """E_0 = omega [N/2 + lambda N(N-1)/2]."""
        N, lam = self.N, self.lam
        return self.omega * (0.5 * N + 0.5 * lam * N * (N - 1))

    def relative_energy_numerical(self, n_grid=4000, rmax=12.0):
        """Ground-state energy of the relative motion for N = 2, on a grid.

        With x = x_1 - x_2 the relative Hamiltonian is

            h = -d^2/dx^2 + (omega^2/4) x^2 + lambda(lambda-1)/x^2 ,

        whose ground-state energy is omega(lambda + 1/2).  We solve it on the
        half line x > 0, where the wave function |x|^lambda exp(-omega x^2/4)
        vanishes at the origin for lambda > 0.
        """
        lam, w = self.lam, self.omega
        x, h = np.linspace(0.0, rmax, n_grid + 2, retstep=True)
        x = x[1:-1]                                    # interior points
        diag = 2.0 / h**2 + 0.25 * w**2 * x**2 + lam * (lam - 1.0) / x**2
        off = np.full(len(x) - 1, -1.0 / h**2)
        H = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)
        return float(np.linalg.eigvalsh(H)[0])

    def relative_energy_exact(self):
        """omega (lambda + 1/2) for the relative motion of two particles."""
        return self.omega * (self.lam + 0.5)

    def jastrow(self, positions):
        """The unnormalised ground-state wave function."""
        x = np.asarray(positions, dtype=float)
        product = 1.0
        for i in range(len(x)):
            for j in range(i + 1, len(x)):
                product *= abs(x[i] - x[j]) ** self.lam
        return product * np.exp(-0.5 * self.omega * np.sum(x**2))


# ---------------------------------------------------------------------------
def _demo():
    print("=" * 74)
    print("1. The Lipkin model, J = 2 block")
    print("=" * 74)
    for eps, V, W in ((2.0, -1.0/3.0, -1.0/4.0), (2.0, -4.0/3.0, -1.0)):
        model = LipkinModel(J=2.0, eps=eps, V=V, W=W, N=4)
        H = model.matrix()
        values, vectors = np.linalg.eigh(H)
        print(f"\neps = {eps}, V = {V:+.4f}, W = {W:+.4f}")
        print("Hamiltonian matrix:")
        for row in H:
            print("   " + "  ".join(f"{x:8.4f}" for x in row))
        print("eigenvalues:", np.array2string(values, precision=5))
        psi = vectors[:, 0]
        big = [(abs(c), k) for k, c in enumerate(psi) if abs(c) > 1e-3]
        labels = {k: f"|2,{int(m)}>" for k, m in enumerate(model.jz)}
        text = "  ".join(f"{abs(psi[k]):.5f} {labels[k]}" for _, k in
                         sorted(big, reverse=True))
        print(f"ground state E_0 = {values[0]:.5f}")
        print("   " + text)

    print()
    print("=" * 74)
    print("2. The pairing model, four levels and two pairs")
    print("=" * 74)
    print(f"{'g':>7s} {'E_ref':>12s} {'E_0 (exact)':>14s} "
          f"{'correlation':>14s} {'dimension':>11s}")
    for g in (-1.0, -0.5, 0.0, 0.5, 1.0):
        model = PairingModel(levels=4, pairs=2, g=g)
        e_ref = model.reference_energy()
        e_0 = model.ground_state_energy()
        print(f"{g:7.2f} {e_ref:12.6f} {e_0:14.8f} {e_0 - e_ref:14.8f} "
              f"{model.dim:11d}")
    print("\nAt g = 0 the reference determinant is exact; the correlation")
    print("energy grows quadratically in g for small g, as second-order")
    print("perturbation theory predicts.")

    print()
    print("=" * 74)
    print("3. The pairing model with a particle-hole term")
    print("=" * 74)
    print("The particle-hole term breaks pairs, so the seniority-zero space is")
    print("no longer closed and we work with all 36 determinants of the")
    print("S_z = 0 sector.  At f = 0 the two descriptions must agree:")
    print()
    print(f"{'g':>7s} {'PairingModel':>16s} {'PairingPHModel(f=0)':>22s} "
          f"{'difference':>13s}")
    for g in (-1.0, -0.5, 0.0, 0.5, 1.0):
        seniority_zero = PairingModel(levels=4, pairs=2, g=g) \
            .ground_state_energy()
        full = PairingPHModel(levels=4, particles=4, g=g, f=0.0) \
            .ground_state_energy()
        print(f"{g:7.2f} {seniority_zero:16.8f} {full:22.8f} "
              f"{abs(seniority_zero - full):13.1e}")

    print()
    print("Switching f on lowers the ground state and mixes in states that the")
    print("pairing interaction alone can never reach:")
    print()
    print(f"{'g':>6s} {'f':>7s} {'E_ref':>10s} {'E_0':>14s} "
          f"{'E_corr':>13s} {'dim':>6s}")
    for g, f in ((1.0, 0.0), (1.0, 0.05), (1.0, 0.20), (1.0, 0.50)):
        model = PairingPHModel(levels=4, particles=4, g=g, f=f)
        e_ref = model.reference_energy()
        e_0 = model.ground_state_energy()
        print(f"{g:6.2f} {f:7.2f} {e_ref:10.4f} {e_0:14.8f} "
              f"{e_0 - e_ref:13.8f} {model.dim:6d}")

    print()
    print("Which determinants does the Hamiltonian connect to the reference?")
    print("Excitation rank is counted in particle-hole pairs, seniority is the")
    print("number of unpaired particles.")
    for f in (0.0, 0.05):
        model = PairingPHModel(levels=4, particles=4, g=1.0, f=f)
        print(f"\n  g = 1, f = {f}:")
        for (rank, seniority), values in model.couplings_to_reference().items():
            unique = sorted({round(v, 8) for v in values})
            print(f"    {rank}p-{rank}h, seniority {seniority}: "
                  f"{len(values):2d} determinants, "
                  f"matrix elements {unique}")
    print()
    print("With f = 0 only the four seniority-zero 2p-2h determinants -- the")
    print("ones reached by moving a whole pair -- are connected, each by -g/2.")
    print("With f non-zero the pair-breaking terms open two further channels:")
    print("eight 1p-1h determinants and eight 2p-2h determinants with two")
    print("unpaired particles, all with matrix elements +/- f/2.  The 1p-1h")
    print("channel is what makes the Tamm-Dancoff and random-phase")
    print("approximations of chapter 7 non-trivial.")

    print()
    print("=" * 74)
    print("4. The Hubbard chain at half filling")
    print("=" * 74)
    print(f"{'U/t':>6s} {'dimension':>11s} {'E_0/t':>14s} "
          f"{'E_0/(Nt)':>12s} {'double occ.':>13s}")
    for U in (0.0, 2.0, 4.0, 8.0, 16.0):
        model = HubbardChain(sites=4, n_up=2, n_down=2, t=1.0, U=U)
        e0 = model.ground_state_energy()
        print(f"{U:6.1f} {model.dim:11d} {e0:14.8f} {e0/model.n:12.6f} "
              f"{model.double_occupancy():13.6f}")
    print("\nDouble occupancy falls as U grows: the system crosses over from")
    print("a metal to a Mott insulator, and the low-energy physics becomes")
    print("that of the Heisenberg model with J = 4t^2/U.")

    print()
    print("=" * 74)
    print("5. Strong coupling: the Hubbard chain against 4t^2/U")
    print("=" * 74)
    print("At half filling and U -> infinity no hop is possible and E_0 -> 0.")
    print("The leading correction is the effective spin Hamiltonian")
    print("  H_eff = J sum_<ij> ( S_i.S_j - n_i n_j / 4 ),   J = 4 t^2 / U,")
    print("and for the four-site ring the Heisenberg part gives -2J while the")
    print("constant -n_i n_j/4 on four bonds gives -J, so E_0 -> -3J.")
    print()
    print(f"{'U/t':>6s} {'E_0':>18s} {'J = 4t^2/U':>14s} {'E_0 / J':>12s}")
    for U in (8.0, 16.0, 32.0, 64.0, 128.0):
        e0 = HubbardChain(sites=4, n_up=2, n_down=2, t=1.0, U=U) \
            .ground_state_energy()
        J = 4.0 / U
        print(f"{U:6.1f} {e0:18.10f} {J:14.8f} {e0/J:12.5f}")
    print("\nThe ratio converges to -3, which is the prediction. The Hubbard")
    print("model at strong coupling really is the Heisenberg antiferromagnet.")

    print()
    print("=" * 74)
    print("6. The Heisenberg model")
    print("=" * 74)
    print("The same four-site ring, now as a spin model in its own right.")
    print("The Heisenberg part alone gives -2J, and the -n_i n_j/4 term on")
    print("four bonds accounts for the remaining -J of the Hubbard limit:")
    print()
    ring = HeisenbergChain(sites=4, J=1.0)
    print(f"  Heisenberg four-site ring, E_0 = {ring.ground_state_energy():.8f} J")
    print(f"  nearest-neighbour <S_i.S_j>  = {ring.correlation(0, 1):+.8f}")
    print(f"  next-nearest     <S_i.S_j>   = {ring.correlation(0, 2):+.8f}")
    print()
    print(f"{'U/t':>6s} {'E_0(Hubbard)/J':>16s} {'+1 (remove -J)':>16s}")
    for U in (8.0, 16.0, 32.0, 64.0, 128.0):
        e0 = HubbardChain(sites=4, n_up=2, n_down=2, t=1.0, U=U) \
            .ground_state_energy()
        J = 4.0 / U
        print(f"{U:6.1f} {e0/J:16.6f} {e0/J + 1.0:16.6f}")
    print("  -> -2, the Heisenberg ground-state energy of the same ring.")

    print()
    print("Longer rings, against the Bethe-ansatz value J (1/4 - ln 2) for the")
    print("infinite chain:")
    print()
    print(f"{'N':>4s} {'dimension':>11s} {'E_0/J':>14s} {'E_0/(NJ)':>13s} "
          f"{'deviation':>12s}")
    bethe = bethe_energy_per_site()
    for n in (4, 6, 8, 10, 12, 14, 16):
        chain = HeisenbergChain(sites=n, J=1.0)
        e0 = chain.ground_state_energy()
        print(f"{n:4d} {chain.dim:11d} {e0:14.8f} {e0/n:13.8f} "
              f"{e0/n - bethe:+12.6f}")
    print(f"{'inf':>4s} {'--':>11s} {'--':>14s} {bethe:13.8f} "
          f"{0.0:+12.6f}")
    print()
    print("The finite-size energies approach the Bethe value from below, and")
    print("the deviation falls off roughly as 1/N^2 -- slowly, which is why")
    print("the one-dimensional chain needed an exact solution rather than a")
    print("diagonalisation.")

    print()
    print("=" * 74)
    print("7. The Calogero model")
    print("=" * 74)
    print(f"{'N':>4s} {'lambda':>8s} {'E_0 exact':>14s}")
    for N in (2, 3, 5, 10):
        for lam in (1.0, 2.0):
            model = CalogeroModel(N=N, lam=lam)
            print(f"{N:4d} {lam:8.1f} {model.exact_energy():14.6f}")
    print()
    print("Two particles: the relative motion solved on a grid against the")
    print("analytic result omega(lambda + 1/2).")
    print(f"{'lambda':>8s} {'numerical':>14s} {'exact':>14s} {'error':>12s}")
    for lam in (1.0, 1.5, 2.0, 3.0):
        model = CalogeroModel(N=2, lam=lam)
        numerical = model.relative_energy_numerical()
        exact = model.relative_energy_exact()
        print(f"{lam:8.1f} {numerical:14.8f} {exact:14.8f} "
              f"{abs(numerical - exact):12.2e}")
    print()
    print("At lambda = 1 the interaction vanishes and the ground state is the")
    print("Slater determinant of free fermions; the Jastrow factor")
    print("prod |x_i - x_j| is exactly the Vandermonde determinant.")


if __name__ == "__main__":
    _demo()
