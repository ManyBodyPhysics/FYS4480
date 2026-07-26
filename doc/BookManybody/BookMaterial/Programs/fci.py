"""
Full configuration interaction, and what happens when it is truncated.

Companion code to chapter 5 of *Quantum mechanics for Many-particle Systems*.

The Hamiltonian is built in the basis of all Slater determinants that can be
formed from a truncated set of single-particle states, and diagonalised.  Every
determinant is classified by its excitation level relative to a reference, so
that the block structure of the matrix can be displayed and the truncated
schemes -- CIS, CID, CISD and so on -- obtained by restricting the basis.

    SlaterBasis      -- determinants as bit strings, with fermionic phases
    PairingFCI       -- the pairing model in the full space, not just the
                        seniority-zero subspace
    HubbardFCI       -- the Hubbard ring in the momentum basis, where the
                        free-fermion determinant is the natural reference
    hilbert_growth   -- how the dimension explodes

Author: Morten Hjorth-Jensen
"""

from itertools import combinations
from math import comb

import numpy as np


# ---------------------------------------------------------------------------
#  Slater determinants as bit strings
# ---------------------------------------------------------------------------
def popcount(x):
    return bin(x).count("1")


def annihilate(state, p):
    """a_p |state>, returning (sign, new state) or (0, 0) if it vanishes."""
    if not (state >> p) & 1:
        return 0, 0
    sign = -1 if popcount(state & ((1 << p) - 1)) & 1 else 1
    return sign, state ^ (1 << p)


def create(state, p):
    """a_p^dagger |state>, returning (sign, new state) or (0, 0)."""
    if (state >> p) & 1:
        return 0, 0
    sign = -1 if popcount(state & ((1 << p) - 1)) & 1 else 1
    return sign, state | (1 << p)


def apply_string(state, operators):
    """Apply a list of (index, dagger) operators, rightmost first.

    The convention matches the text: the phase is the parity of the occupied
    orbitals below the one acted on, exactly the Jordan-Wigner string.
    """
    sign = 1
    for p, dagger in reversed(operators):
        s, state = (create(state, p) if dagger else annihilate(state, p))
        if s == 0:
            return 0, 0
        sign *= s
    return sign, state


class SlaterBasis:
    """All determinants with a fixed number of particles in n_orbitals."""

    def __init__(self, n_orbitals, n_particles):
        self.n_orbitals = n_orbitals
        self.n_particles = n_particles
        self.states = []
        for occ in combinations(range(n_orbitals), n_particles):
            bits = 0
            for p in occ:
                bits |= 1 << p
            self.states.append(bits)
        self.states.sort()
        self.index = {s: i for i, s in enumerate(self.states)}
        self.dim = len(self.states)

    def excitation_levels(self, reference):
        """Number of particles promoted above the reference, per determinant."""
        return np.array([popcount(s & ~reference) for s in self.states])


# ---------------------------------------------------------------------------
class PairingFCI:
    """The pairing model in the full space of Slater determinants.

        H = xi sum_{p sigma} (p-1) a^+_{p sigma} a_{p sigma}
          - (g/2) sum_{pq} a^+_{p+} a^+_{p-} a_{q-} a_{q+}

    Spin-orbital (p, sigma) is stored in bit 2(p-1) for sigma = +1 and
    2(p-1)+1 for sigma = -1.  Unlike the seniority-zero treatment of
    chapter 4, broken pairs are included here, so that the full CI space and
    its excitation-level structure can be examined.
    """

    def __init__(self, levels=4, n_particles=4, g=1.0, xi=1.0):
        self.levels = levels
        self.g = g
        self.xi = xi
        self.basis = SlaterBasis(2 * levels, n_particles)
        self.reference = (1 << n_particles) - 1        # lowest orbitals filled
        self.levels_of = self.basis.excitation_levels(self.reference)

    # ------------------------------------------------------------------
    def _sp_energy(self, orbital):
        return self.xi * (orbital // 2)                # (p-1) with p = 1, 2, ...

    def matrix(self):
        dim = self.basis.dim
        H = np.zeros((dim, dim))
        for col, state in enumerate(self.basis.states):
            # one-body part: diagonal
            energy = sum(self._sp_energy(p) for p in range(self.basis.n_orbitals)
                         if (state >> p) & 1)
            H[col, col] += energy
            # pairing interaction
            for q in range(self.levels):
                for p in range(self.levels):
                    ops = [(2*p, True), (2*p + 1, True),
                           (2*q + 1, False), (2*q, False)]
                    sign, new = apply_string(state, ops)
                    if sign:
                        H[self.basis.index[new], col] += -0.5 * self.g * sign
        return H

    # ------------------------------------------------------------------
    def energies(self, max_excitation=None):
        """Ground-state energy with the basis truncated at an excitation level.

        max_excitation = 0 gives the reference determinant alone, 1 adds the
        singles (CIS), 2 the doubles (CID/CISD) and so on; None keeps the
        whole space, which is full CI.
        """
        H = self.matrix()
        if max_excitation is None:
            keep = np.ones(self.basis.dim, dtype=bool)
        else:
            keep = self.levels_of <= max_excitation
        sub = H[np.ix_(keep, keep)]
        return float(np.linalg.eigvalsh(sub)[0]), int(keep.sum())

    def reference_energy(self):
        H = self.matrix()
        i = self.basis.index[self.reference]
        return float(H[i, i])

    def block_structure(self):
        """Which excitation blocks of the Hamiltonian are non-zero."""
        H = self.matrix()
        n = int(self.levels_of.max()) + 1
        table = np.zeros((n, n), dtype=bool)
        for a in range(n):
            for b in range(n):
                rows = self.levels_of == a
                cols = self.levels_of == b
                if rows.any() and cols.any():
                    table[a, b] = np.abs(H[np.ix_(rows, cols)]).max() > 1e-12
        return table


# ---------------------------------------------------------------------------
class HubbardFCI:
    """The Hubbard ring in the momentum basis.

        H = sum_{k sigma} eps_k n_{k sigma}
          + (U/N) sum_{k k' q} c^+_{k+q up} c^+_{k'-q down} c_{k' down} c_{k up}

    with eps_k = -2t cos k.  In this basis the free-fermion determinant is the
    natural reference, the one-body part is diagonal, and the excitation-level
    classification means what it says.  Spin-orbital (k, sigma) occupies bit
    k for spin up and N + k for spin down.
    """

    def __init__(self, sites=6, n_up=3, n_down=3, t=1.0, U=4.0):
        self.n = sites
        self.t = t
        self.U = U
        self.eps = -2.0 * t * np.cos(2.0 * np.pi * np.arange(sites) / sites)

        up_basis = SlaterBasis(sites, n_up)
        dn_basis = SlaterBasis(sites, n_down)
        self.states = [(u, d) for u in up_basis.states for d in dn_basis.states]
        self.index = {s: i for i, s in enumerate(self.states)}
        self.dim = len(self.states)

        order = np.argsort(self.eps, kind="stable")
        ref_up = sum(1 << int(k) for k in order[:n_up])
        ref_dn = sum(1 << int(k) for k in order[:n_down])
        self.reference = (ref_up, ref_dn)
        self.levels_of = np.array(
            [popcount(u & ~ref_up) + popcount(d & ~ref_dn)
             for (u, d) in self.states])

    # ------------------------------------------------------------------
    def _combined(self, up, dn):
        """The two spin blocks as one bit string: up occupies the low bits."""
        return up | (dn << self.n)

    def matrix(self):
        n, dim = self.n, self.dim
        H = np.zeros((dim, dim))
        for col, (up, dn) in enumerate(self.states):
            state = self._combined(up, dn)
            H[col, col] += sum(self.eps[k] for k in range(n)
                               if (up >> k) & 1)
            H[col, col] += sum(self.eps[k] for k in range(n)
                               if (dn >> k) & 1)
            for k in range(n):
                for kp in range(n):
                    for q in range(n):
                        ops = [((k + q) % n, True),
                               (n + (kp - q) % n, True),
                               (n + kp, False),
                               (k, False)]
                        sign, new = apply_string(state, ops)
                        if sign:
                            new_up = new & ((1 << n) - 1)
                            new_dn = new >> n
                            row = self.index[(new_up, new_dn)]
                            H[row, col] += self.U / n * sign
        return H

    def energies(self, max_excitation=None):
        H = self.matrix()
        if max_excitation is None:
            keep = np.ones(self.dim, dtype=bool)
        else:
            keep = self.levels_of <= max_excitation
        sub = H[np.ix_(keep, keep)]
        return float(np.linalg.eigvalsh(sub)[0]), int(keep.sum())

    def reference_energy(self):
        H = self.matrix()
        return float(H[self.index[self.reference], self.index[self.reference]])


# ---------------------------------------------------------------------------
class NonInteractingCopies:
    """Two identical, non-interacting pairing systems, A and B.

    Each subsystem has two levels of energies 0 and xi and holds two
    particles, and the pairing interaction acts only *within* a subsystem.
    The exact energy of the pair of systems must then be twice the exact
    energy of one -- this is what size consistency means.

    Truncated CI fails the test.  Describing both subsystems as doubly
    excited at the same time is a *quadruple* excitation of the combined
    system, which CID does not contain, so the combined answer is worse than
    twice the separate one.
    """

    def __init__(self, g=0.5, xi=1.0, subsystems=2):
        self.g = g
        self.xi = xi
        self.subsystems = subsystems
        self.levels = 2 * subsystems              # two levels per subsystem
        n_particles = 2 * subsystems              # one pair per subsystem
        self.basis = SlaterBasis(2 * self.levels, n_particles)
        # the reference has one pair in the lower level of every subsystem
        ref = 0
        for s in range(subsystems):
            low = 2 * s                            # level index of the lower
            ref |= (1 << (2 * low)) | (1 << (2 * low + 1))
        self.reference = ref
        self.levels_of = self.basis.excitation_levels(ref)

    def _subsystem(self, level):
        return level // 2

    def _sp_energy(self, orbital):
        """Level energies repeat: 0, xi, 0, xi, ... one pair per subsystem."""
        level = orbital // 2
        return self.xi * (level % 2)

    def matrix(self):
        dim = self.basis.dim
        H = np.zeros((dim, dim))
        for col, state in enumerate(self.basis.states):
            H[col, col] += sum(self._sp_energy(p)
                               for p in range(self.basis.n_orbitals)
                               if (state >> p) & 1)
            for q in range(self.levels):
                for p in range(self.levels):
                    if self._subsystem(p) != self._subsystem(q):
                        continue                   # no coupling between A and B
                    ops = [(2*p, True), (2*p + 1, True),
                           (2*q + 1, False), (2*q, False)]
                    sign, new = apply_string(state, ops)
                    if sign:
                        H[self.basis.index[new], col] += -0.5 * self.g * sign
        return H

    def energies(self, max_excitation=None):
        H = self.matrix()
        if max_excitation is None:
            keep = np.ones(self.basis.dim, dtype=bool)
        else:
            keep = self.levels_of <= max_excitation
        return float(np.linalg.eigvalsh(H[np.ix_(keep, keep)])[0]), \
            int(keep.sum())


# ---------------------------------------------------------------------------
def hilbert_growth():
    """Dimensions of some many-body spaces."""
    rows = []
    for n_sp, N, label in ((10, 5, "toy model"),
                           (20, 10, "small shell model"),
                           (40, 8, "oxygen-16 neutrons, 4 shells"),
                           (40, 20, "half filling, 40 states"),
                           (80, 40, "half filling, 80 states")):
        rows.append((label, n_sp, N, comb(n_sp, N)))
    return rows


# ---------------------------------------------------------------------------
def demo_full_space():
    print("=" * 74)
    print("1. The pairing model in the full CI space")
    print("=" * 74)
    model = PairingFCI(levels=4, n_particles=4, g=1.0)
    print(f"four levels, four particles: {2*model.levels} spin-orbitals,")
    print(f"dimension of the full space  = {model.basis.dim} "
          f"= C({2*model.levels},4)")
    print(f"seniority-zero subspace      = {comb(4, 2)} "
          f"(chapter 4 worked here)")
    counts = np.bincount(model.levels_of)
    print("\ndeterminants by excitation level relative to the reference:")
    for level, count in enumerate(counts):
        print(f"   {level}p-{level}h: {count:4d}")



def demo_block_structure():
    model = PairingFCI(levels=4, n_particles=4, g=1.0)
    print("=" * 74)
    print("2. The block structure of the Hamiltonian matrix")
    print("=" * 74)
    table = model.block_structure()
    n = table.shape[0]
    header = "        " + "".join(f"{k}p-{k}h  " for k in range(n))
    print(header)
    for a in range(n):
        row = "".join("  x    " if table[a, b] else "  0    " for b in range(n))
        print(f" {a}p-{a}h  {row}")
    print("\nA two-body operator connects determinants differing by at most")
    print("two single-particle states, so the matrix is banded in the")
    print("excitation level: the Condon-Slater rules of chapter 2 read as a")
    print("statement about sparsity.")



def demo_truncations():
    print("=" * 74)
    print("3. Truncating the CI expansion")
    print("=" * 74)
    print(f"{'g':>6s} {'E_ref':>11s} {'CIS':>13s} {'CID':>13s} "
          f"{'CISD':>13s} {'FCI':>13s}")
    for g in (0.25, 0.5, 1.0, 2.0):
        m = PairingFCI(levels=4, n_particles=4, g=g)
        e_ref = m.reference_energy()
        e1, _ = m.energies(1)
        e2, _ = m.energies(2)
        e_fci, _ = m.energies(None)
        print(f"{g:6.2f} {e_ref:11.6f} {e1:13.8f} {e2:13.8f} "
              f"{e2:13.8f} {e_fci:13.8f}")
    m = PairingFCI(levels=4, n_particles=4, g=1.0)
    print("\nbasis sizes: ", end="")
    for k in (0, 1, 2, 3, 4):
        print(f"{k}p-{k}h and below: {m.energies(k)[1]:3d}   ", end="")
    print(f"full: {m.basis.dim}")
    print("\nSingles alone change nothing: the pairing interaction cannot")
    print("break a pair, so <Phi_0|H|Phi_i^a> = 0 and CIS = reference.")
    print("Doubles capture almost all of the correlation energy, and CISD")
    print("equals CID here for the same reason.")



def demo_size_consistency():
    print("=" * 74)
    print("4. Size consistency")
    print("=" * 74)
    print("One subsystem: two levels, two particles.  Then two identical")
    print("copies with no interaction between them.  A size-consistent")
    print("method must give exactly twice the energy.")
    print()
    print(f"{'g':>6s} {'E(A)':>13s} {'2 E(A)':>13s} {'E(AB)':>13s} "
          f"{'error':>12s}   method")
    for g in (0.25, 0.5, 1.0):
        one = NonInteractingCopies(g=g, subsystems=1)
        two = NonInteractingCopies(g=g, subsystems=2)
        for label, level in (("FCI", None), ("CID", 2)):
            e1, _ = one.energies(level)
            e2, _ = two.energies(level)
            print(f"{g:6.2f} {e1:13.8f} {2*e1:13.8f} {e2:13.8f} "
                  f"{e2 - 2*e1:12.2e}   {label}")
    print()
    print("Full CI passes to machine precision.  CID does not: for the")
    print("combined system, both subsystems being doubly excited at once is")
    print("a quadruple excitation, and CID has no room for it.  The error")
    print("grows with the number of subsystems, so truncated CI becomes")
    print("useless for large systems.  This is the failure that")
    print("coupled-cluster theory repairs by exponentiating the cluster")
    print("operator.")



def demo_hubbard():
    print("=" * 74)
    print("5. The Hubbard ring in the momentum basis")
    print("=" * 74)
    print(f"{'U/t':>6s} {'E_ref':>12s} {'CIS':>13s} {'CID':>13s} "
          f"{'CISD':>13s} {'FCI':>13s} {'% recovered':>12s}")
    for U in (1.0, 2.0, 4.0, 8.0):
        h = HubbardFCI(sites=6, n_up=3, n_down=3, t=1.0, U=U)
        e_ref = h.reference_energy()
        e1, _ = h.energies(1)
        e2, _ = h.energies(2)
        e3, _ = h.energies(3)
        e_fci, _ = h.energies(None)
        pct = 100.0 * (e2 - e_ref) / (e_fci - e_ref)
        print(f"{U:6.1f} {e_ref:12.6f} {e1:13.8f} {e2:13.8f} "
              f"{e3:13.8f} {e_fci:13.8f} {pct:11.2f}%")
    h = HubbardFCI(sites=6, n_up=3, n_down=3, U=4.0)
    print(f"\nsix sites at half filling: dimension "
          f"{h.dim} = C(6,3)^2")
    counts = np.bincount(h.levels_of)
    print("determinants by excitation level:",
          "  ".join(f"{k}p-{k}h: {c}" for k, c in enumerate(counts)))
    print("\nAs U grows the reference determinant becomes a worse starting")
    print("point and the doubles recover a smaller fraction of the")
    print("correlation energy -- the signature of strong correlation.")



def demo_exponential_wall():
    print("=" * 74)
    print("6. The exponential wall")
    print("=" * 74)
    print(f"{'system':>32s} {'n_sp':>6s} {'N':>4s} {'dimension':>16s}")
    for label, n_sp, N, dim in hilbert_growth():
        print(f"{label:>32s} {n_sp:6d} {N:4d} {dim:16.3e}")
    print()
    print("At half filling the binomial coefficient becomes exponential:")
    print(f"{'n':>6s} {'C(n, n/2)':>16s} {'2^n / sqrt(n)':>16s}")
    for n in (10, 20, 40, 80, 160):
        print(f"{n:6d} {comb(n, n//2):16.3e} "
              f"{2.0**n/np.sqrt(n):16.3e}")
    from math import comb as _c
    neutrons = _c(40, 8)
    print(f"\nOxygen-16 in the four lowest major shells (0s, 0p, 1s0d, 1p0f)")
    print(f"has 40 single-particle states for each species, so the eight")
    print(f"neutrons give C(40,8) = {neutrons:.3e} determinants and the full")
    print(f"proton-neutron space {neutrons**2:.3e}.  Symmetries reduce this,")
    print("but adding one more shell or relaxing the truncation restores it.")
    print("Direct diagonalisation reaches about 1e5, Lanczos about 1e10.")
    print("Beyond that the expansion itself has to be truncated.")


def _demo():
    for f in (demo_full_space, demo_block_structure, demo_truncations,
              demo_size_consistency, demo_hubbard, demo_exponential_wall):
        f()
        print()


if __name__ == "__main__":
    _demo()
