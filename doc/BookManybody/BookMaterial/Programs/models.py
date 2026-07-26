"""
The four model Hamiltonians of chapter 4.

Companion code to chapter 4 of *Quantum mechanics for Many-particle Systems*.

Each class builds the Hamiltonian of one schematic model and diagonalises it,
so that every number quoted in the chapter can be reproduced:

    LipkinModel      -- two levels, degeneracy Omega, quasispin SU(2)
    PairingModel     -- L doubly degenerate levels, pair-transfer interaction
    HubbardChain     -- N sites, hopping t and on-site repulsion U
    CalogeroModel    -- N particles on a line, inverse-square interaction

The first three are diagonalised in a fixed symmetry sector; the fourth is
exactly solvable and we check the analytic ground-state energy against a
numerical solution of the two-body problem.

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
    print("3. The Hubbard chain at half filling")
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
    print("4. Strong coupling: the Hubbard chain against 4t^2/U")
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
    print("5. The Calogero model")
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
