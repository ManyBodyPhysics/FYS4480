"""
The Jordan-Wigner transformation, verified.

Companion code to chapters 3 and 4 of *Quantum mechanics for Many-particle
Systems*.

A quantum computer acts on qubits, not on fermionic Fock states.  The
occupation-number basis is already binary, so the states map across trivially;
what is not trivial is the operators, because the fermionic minus signs have
to be reproduced by something.  The Jordan-Wigner transformation supplies them
with a string of Pauli Z operators,

    a_p^dagger = (Z_0 Z_1 ... Z_{p-1}) (X_p - i Y_p)/2 ,
    a_p        = (Z_0 Z_1 ... Z_{p-1}) (X_p + i Y_p)/2 ,

and this program checks that the result really does satisfy the fermionic
anticommutation relations, and works out what the many-body Hamiltonians of
the book look like once the transformation has been applied.

Everything here is done with explicit matrices in a small Fock space, so
nothing has to be taken on trust.  The Pauli decomposition is exact: for M
qubits the 4^M Pauli strings are orthogonal under the Hilbert-Schmidt inner
product, so

    H = sum_P c_P P ,       c_P = Tr(P H) / 2^M ,

and the coefficients are read off by a trace.  That is expensive -- 4^M
traces -- so it is only used for the small examples of the text.

Author: Morten Hjorth-Jensen
"""

from itertools import combinations, product

import numpy as np

# ---------------------------------------------------------------------------
#  Pauli matrices and Pauli strings
# ---------------------------------------------------------------------------
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

PAULI = {"I": I2, "X": X, "Y": Y, "Z": Z}


def kron(matrices):
    """Tensor product of a list of matrices, left to right."""
    out = np.array([[1.0 + 0j]])
    for matrix in matrices:
        out = np.kron(out, matrix)
    return out


def pauli_string(label):
    """The matrix of a Pauli string given as a word such as 'XZZY'."""
    return kron([PAULI[c] for c in label])


def single(op, site, n_qubits):
    """One Pauli operator on `site`, identity elsewhere."""
    return kron([op if k == site else I2 for k in range(n_qubits)])


def decompose(matrix, tol=1e-10):
    """Expand a matrix in Pauli strings: {label: coefficient}.

    Uses the orthogonality Tr(P Q) = 2^M delta_PQ of the Pauli basis.  Only
    the coefficients larger than `tol` are kept, and they are real for a
    Hermitian input.
    """
    dim = matrix.shape[0]
    n_qubits = int(round(np.log2(dim)))
    out = {}
    for word in product("IXYZ", repeat=n_qubits):
        label = "".join(word)
        coefficient = np.trace(pauli_string(label) @ matrix) / dim
        if abs(coefficient) > tol:
            out[label] = (float(coefficient.real)
                          if abs(coefficient.imag) < tol
                          else complex(coefficient))
    return out


def show(decomposition, indent="   "):
    """Print a Pauli decomposition, sorted by decreasing weight."""
    for label, coefficient in sorted(decomposition.items(),
                                     key=lambda kv: (-abs(kv[1]), kv[0])):
        if isinstance(coefficient, complex):
            print(f"{indent}{coefficient:+.6f}  {label}")
        else:
            print(f"{indent}{coefficient:+.6f}  {label}")


def weight(label):
    """Number of non-identity factors in a Pauli string."""
    return sum(1 for c in label if c != "I")


# ---------------------------------------------------------------------------
#  The Jordan-Wigner map
# ---------------------------------------------------------------------------
class JordanWigner:
    """Fermion operators on `n_modes` orbitals, as 2^n x 2^n qubit matrices.

    Qubit p carries the occupation of orbital p, with |0> empty and |1>
    occupied, so that n_p = (1 - Z_p)/2.  Orbital 0 is the leftmost factor of
    the tensor product, which fixes the direction of the parity string.
    """

    def __init__(self, n_modes):
        self.n = n_modes
        self.dim = 1 << n_modes

    def _string(self, p):
        """The parity string Z_0 Z_1 ... Z_{p-1}, identity beyond."""
        return kron([Z if k < p else I2 for k in range(self.n)])

    def create(self, p):
        """a_p^dagger = (prod_{j<p} Z_j) (X_p - i Y_p)/2."""
        sigma_plus = single((X - 1j * Y) / 2.0, p, self.n)
        return self._string(p) @ sigma_plus

    def annihilate(self, p):
        """a_p = (prod_{j<p} Z_j) (X_p + i Y_p)/2."""
        return self.create(p).conj().T

    def number(self, p):
        return self.create(p) @ self.annihilate(p)

    def identity(self):
        return np.eye(self.dim, dtype=complex)

    # -- the checks that justify the whole construction -------------------
    def anticommutator_errors(self):
        """max |{a_p, a_q^dagger} - delta_pq|, |{a_p,a_q}|, |{a_p^+,a_q^+}|."""
        worst = 0.0
        eye = self.identity()
        for p in range(self.n):
            for q in range(self.n):
                a_p, a_q = self.annihilate(p), self.annihilate(q)
                c_q = self.create(q)
                target = eye if p == q else 0.0 * eye
                worst = max(worst, np.abs(a_p @ c_q + c_q @ a_p
                                          - target).max())
                worst = max(worst, np.abs(a_p @ a_q + a_q @ a_p).max())
                c_p = self.create(p)
                worst = max(worst, np.abs(c_p @ c_q + c_q @ c_p).max())
        return float(worst)

    # -- building Hamiltonians -------------------------------------------
    def one_body(self, h):
        """sum_pq h_pq a_p^dagger a_q."""
        out = np.zeros((self.dim, self.dim), dtype=complex)
        for p in range(self.n):
            for q in range(self.n):
                if h[p, q] != 0.0:
                    out += h[p, q] * self.create(p) @ self.annihilate(q)
        return out

    def two_body(self, v):
        """(1/4) sum_pqrs v_pqrs a_p^+ a_q^+ a_s a_r, v antisymmetrised."""
        out = np.zeros((self.dim, self.dim), dtype=complex)
        for p, q, r, s in product(range(self.n), repeat=4):
            if v[p, q, r, s] == 0.0:
                continue
            out += 0.25 * v[p, q, r, s] * (self.create(p) @ self.create(q)
                                           @ self.annihilate(s)
                                           @ self.annihilate(r))
        return out

    def term(self, coefficient, creators, annihilators):
        """coefficient * a^+_{c1} a^+_{c2} ... a_{a1} a_{a2} ..."""
        out = coefficient * self.identity()
        for p in creators:
            out = out @ self.create(p)
        for p in annihilators:
            out = out @ self.annihilate(p)
        return out


# ---------------------------------------------------------------------------
#  The pairing Hamiltonian, in the two encodings
# ---------------------------------------------------------------------------
def pairing_full_space(levels, g, xi=1.0):
    """The pairing Hamiltonian on 2L qubits, by Jordan-Wigner.

    Orbital ordering is (1+, 1-, 2+, 2-, ...), so the spin-up and spin-down
    partners of a level are adjacent and the parity strings between them are
    short.  This is the ordering used throughout the book.
    """
    jw = JordanWigner(2 * levels)
    H = np.zeros((jw.dim, jw.dim), dtype=complex)
    for p in range(1, levels + 1):
        for spin in (0, 1):
            H += xi * (p - 1) * jw.number(2 * (p - 1) + spin)
    for p in range(1, levels + 1):
        for q in range(1, levels + 1):
            H += jw.term(-0.5 * g,
                         [2 * (p - 1), 2 * (p - 1) + 1],
                         [2 * (q - 1) + 1, 2 * (q - 1)])
    return H


def pairing_ph_full_space(levels, g, f, xi=1.0):
    """The pairing plus particle-hole Hamiltonian on 2L qubits."""
    jw = JordanWigner(2 * levels)
    H = pairing_full_space(levels, g, xi)
    for p in range(1, levels + 1):
        for q in range(1, levels + 1):
            for r in range(1, levels + 1):
                block = jw.term(-0.5 * f,
                                [2 * (p - 1), 2 * (p - 1) + 1],
                                [2 * (q - 1) + 1, 2 * (r - 1)])
                H += block + block.conj().T
    return H


def pairing_pair_qubits(levels, g, xi=1.0):
    """The pairing Hamiltonian on L qubits, one qubit per pair level.

    In the seniority-zero space a level is either empty or holds a complete
    pair, so one qubit suffices.  With S_p^+ = (X_p - i Y_p)/2 and
    N_p = 1 - Z_p the Hamiltonian becomes

        H = sum_p (2 xi (p-1) - g/2)(1 - Z_p)/... ,

    which the demo below writes out explicitly.  Here we simply build it from
    the spin operators and let `decompose` produce the Pauli form.
    """
    dim = 1 << levels
    H = np.zeros((dim, dim), dtype=complex)
    plus = [single((X - 1j * Y) / 2.0, p, levels) for p in range(levels)]
    minus = [single((X + 1j * Y) / 2.0, p, levels) for p in range(levels)]
    number = [np.eye(dim, dtype=complex) - single(Z, p, levels)
              for p in range(levels)]
    for p in range(levels):
        H += xi * p * number[p]                 # 2 xi (p-1) per pair: N_p = 2
        for q in range(levels):
            H += -0.5 * g * plus[p] @ minus[q]
    return H


def pairing_pair_qubits_closed_form(levels, g, xi=1.0):
    """The closed Pauli form of the pair-qubit pairing Hamiltonian.

        H = sum_p [xi (p-1) - g/4] (1 - Z_p)
            - (g/4) sum_{p<q} (X_p X_q + Y_p Y_q)

    for the normalisation of the book, H_0 = xi sum_{p sigma} (p-1) n_{p
    sigma} and V = -(g/2) sum_{pq} P_p^+ P_q.  Returned as a dictionary of
    Pauli strings, to be compared with `decompose(pairing_pair_qubits(...))`.
    """
    out = {}

    def add(label, value):
        out[label] = out.get(label, 0.0) + value

    for p in range(levels):
        coefficient = xi * p - 0.25 * g          # p here is (level - 1)
        add("I" * levels, coefficient)
        add("".join("Z" if k == p else "I" for k in range(levels)),
            -coefficient)
    for p in range(levels):
        for q in range(p + 1, levels):
            for letter in "XY":
                add("".join(letter if k in (p, q) else "I"
                            for k in range(levels)), -0.25 * g)
    return {k: v for k, v in out.items() if abs(v) > 1e-12}


def closed_form_agrees(levels, g, xi=1.0, tol=1e-10):
    """True if the closed form above reproduces the explicit matrix."""
    predicted = pairing_pair_qubits_closed_form(levels, g, xi)
    computed = decompose(pairing_pair_qubits(levels, g, xi))
    if set(predicted) != set(computed):
        return False
    return all(abs(predicted[k] - computed[k]) < tol for k in predicted)


def two_orbital_model(eps0, eps1, U, t):
    """The two-orbital Hamiltonian of exercise 6 of chapter 3,

        H = eps0 n_0 + eps1 n_1 + U n_0 n_1 + t (a_0^+ a_1 + a_1^+ a_0),

    as a 4 x 4 matrix, together with the Pauli form predicted by hand.
    """
    jw = JordanWigner(2)
    H = (eps0 * jw.number(0) + eps1 * jw.number(1)
         + U * jw.number(0) @ jw.number(1)
         + t * (jw.create(0) @ jw.annihilate(1)
                + jw.create(1) @ jw.annihilate(0)))
    predicted = {
        "II": 0.5 * (eps0 + eps1) + 0.25 * U,
        "ZI": -0.25 * (2 * eps0 + U),
        "IZ": -0.25 * (2 * eps1 + U),
        "ZZ": 0.25 * U,
        "XX": 0.5 * t,
        "YY": 0.5 * t,
    }
    predicted = {k: v for k, v in predicted.items() if abs(v) > 1e-12}
    return H, predicted


def seniority_zero_indices(levels):
    """Rows of the 2^(2L) Fock space that hold only complete pairs."""
    keep = []
    for state in range(1 << (2 * levels)):
        good = True
        for p in range(levels):
            if ((state >> (2 * p)) & 1) != ((state >> (2 * p + 1)) & 1):
                good = False
                break
        if good:
            keep.append(state)
    return keep


def particle_number(n_modes):
    """The diagonal of sum_p n_p in the 2^n Fock space."""
    return np.array([bin(state).count("1") for state in range(1 << n_modes)])


def two_sz(levels):
    """The diagonal of 2 S_z = N_+ - N_- in the 2^(2L) Fock space."""
    out = []
    for state in range(1 << (2 * levels)):
        up = sum((state >> (2 * p)) & 1 for p in range(levels))
        down = sum((state >> (2 * p + 1)) & 1 for p in range(levels))
        out.append(up - down)
    return np.array(out)


def sector_spectrum(H, levels, particles, sz2=0):
    """Eigenvalues of a 2L-qubit Hamiltonian in a fixed (N, S_z) sector."""
    occupancy = particle_number(2 * levels)
    spin = two_sz(levels)
    rows = [k for k in range(1 << (2 * levels))
            if occupancy[k] == particles and (sz2 is None or spin[k] == sz2)]
    return np.sort(np.linalg.eigvalsh(H[np.ix_(rows, rows)]).real), len(rows)


# ---------------------------------------------------------------------------
def _demo():
    print("=" * 74)
    print("1. The transformation is legal: fermionic anticommutators")
    print("=" * 74)
    print("The whole point of the parity string is to make the qubit")
    print("operators anticommute.  Without it they would merely commute on")
    print("different sites, which is the wrong statistics.")
    print()
    for n in (2, 3, 4, 5):
        jw = JordanWigner(n)
        print(f"   {n} modes: max violation of the anticommutators = "
              f"{jw.anticommutator_errors():.1e}")
    print()
    print("For comparison, dropping the parity string breaks them:")
    n = 3
    naive_plus = [single((X - 1j * Y) / 2.0, p, n) for p in range(n)]
    naive_minus = [single((X + 1j * Y) / 2.0, p, n) for p in range(n)]
    worst = max(np.abs(naive_minus[0] @ naive_minus[1]
                       + naive_minus[1] @ naive_minus[0]).max(),
                np.abs(naive_plus[0] @ naive_plus[1]
                       + naive_plus[1] @ naive_plus[0]).max())
    print(f"   max |{{sigma_0, sigma_1}}| without the string = {worst:.3f}"
          "   (should be 0)")

    print()
    print("=" * 74)
    print("2. Number operators and one-body terms")
    print("=" * 74)
    jw = JordanWigner(4)
    print("n_p = a_p^+ a_p:")
    for p in range(4):
        print(f"   p = {p}: ", decompose(jw.number(p)))
    print()
    print("so a diagonal one-body operator is a sum of I and Z only, and the")
    print("qubit Hamiltonian of a non-interacting system is already diagonal.")
    print()
    print("The off-diagonal, hopping-like combination a_p^+ a_q + a_q^+ a_p:")
    for p, q in ((0, 1), (0, 2), (0, 3)):
        matrix = (jw.create(p) @ jw.annihilate(q)
                  + jw.create(q) @ jw.annihilate(p))
        print(f"   (p, q) = ({p}, {q}):")
        show(decompose(matrix))
    print()
    print("The Z's in the middle are the parity string: a hop from q to p has")
    print("to know how many particles sit between them.  Strings grow with")
    print("|p - q|, which is why the orbital ordering matters in practice.")

    print()
    print("=" * 74)
    print("3. Pair creation and the Hermitian pair-transfer operator")
    print("=" * 74)
    for p, q in ((0, 1), (0, 3)):
        matrix = (jw.create(p) @ jw.create(q)
                  + jw.annihilate(q) @ jw.annihilate(p))
        print(f"   a^+_{p} a^+_{q} + a_{q} a_{p}:")
        show(decompose(matrix))
    print()
    print("Note the relative minus sign between the XX and the YY channel --")
    print("the hopping operator of section 2 has a plus.  That single sign is")
    print("the difference between moving a particle and creating a pair.")

    print()
    print("=" * 74)
    print("4. A general two-body term")
    print("=" * 74)
    print("a^+_p a^+_q a_s a_r with all four indices different produces eight")
    print("Pauli strings of weight four, each with coefficient 1/8.")
    print()
    matrix = jw.term(1.0, [0, 1], [3, 2])       # a^+_0 a^+_1 a_3 a_2
    hermitian = matrix + matrix.conj().T
    print("   a^+_0 a^+_1 a_3 a_2 + h.c.:")
    show(decompose(hermitian))
    counts = {}
    for label in decompose(hermitian):
        counts[weight(label)] = counts.get(weight(label), 0) + 1
    print(f"   strings by weight: {dict(sorted(counts.items()))}")
    print()
    print("This is the pairing pair-transfer term in disguise: read the four")
    print("orbitals as (p+, p-, q+, q-) and the operator is P_p^+ P_q + h.c.")
    print()
    print("For a generic two-body Hamiltonian on M modes the string count")
    print("grows as O(M^4), and the longest strings span the whole register")
    print("because of the parity chains:")
    print()
    print(f"{'M':>4s} {'4^M':>8s} {'strings':>9s} {'max weight':>12s}")
    rng = np.random.default_rng(11)
    for n_modes in (2, 3, 4, 5, 6):
        jw_m = JordanWigner(n_modes)
        h = rng.normal(size=(n_modes, n_modes))
        h = h + h.T
        v = rng.normal(size=(n_modes,) * 4)
        v = v - v.transpose(1, 0, 2, 3)
        v = v - v.transpose(0, 1, 3, 2)
        v = v + v.transpose(2, 3, 0, 1)
        decomposition = decompose(jw_m.one_body(h) + jw_m.two_body(v))
        print(f"{n_modes:4d} {4**n_modes:8d} {len(decomposition):9d} "
              f"{max(weight(label) for label in decomposition):12d}")

    print()
    print("=" * 74)
    print("5. The pairing Hamiltonian: full orbital space")
    print("=" * 74)
    levels, g, xi = 2, 1.0, 1.0
    H_full = pairing_full_space(levels, g, xi)
    print(f"L = {levels} levels -> {2*levels} qubits, "
          f"dimension {H_full.shape[0]}")
    decomposition = decompose(H_full)
    print(f"   {len(decomposition)} Pauli strings:")
    show(decomposition)
    counts = {}
    for label in decomposition:
        counts[weight(label)] = counts.get(weight(label), 0) + 1
    print(f"   strings by weight: {dict(sorted(counts.items()))}")

    print()
    print("=" * 74)
    print("6. The pairing Hamiltonian: one qubit per pair level")
    print("=" * 74)
    H_pair = pairing_pair_qubits(levels, g, xi)
    print(f"L = {levels} levels -> {levels} qubits, "
          f"dimension {H_pair.shape[0]}")
    decomposition = decompose(H_pair)
    print(f"   {len(decomposition)} Pauli strings:")
    show(decomposition)
    print()
    print("Compare with the closed form of the text,")
    print("   H = sum_p [xi (p-1) - g/4] (1 - Z_p)")
    print("       - (g/4) sum_{p<q} (X_p X_q + Y_p Y_q),")
    print("evaluated here for xi = 1, g = 1:")
    print(f"   closed form reproduces the decomposition: "
          f"{closed_form_agrees(levels, g, xi)}")

    print()
    print("=" * 74)
    print("7. The two encodings describe the same physics")
    print("=" * 74)
    print("The compact encoding is a restriction, not an approximation: its")
    print("spectrum is exactly the seniority-zero part of the full one, at")
    print("fixed particle number.")
    print()
    for levels in (2, 3):
        H_full = pairing_full_space(levels, g=1.0)
        H_pair = pairing_pair_qubits(levels, g=1.0)
        keep = seniority_zero_indices(levels)
        occupancy = particle_number(2 * levels)
        for pairs in range(1, levels + 1):
            rows = [k for k in keep if occupancy[k] == 2 * pairs]
            block = H_full[np.ix_(rows, rows)]
            full_values = np.sort(np.linalg.eigvalsh(block).real)
            pair_rows = [k for k in range(1 << levels)
                         if bin(k).count("1") == pairs]
            pair_block = H_pair[np.ix_(pair_rows, pair_rows)]
            pair_values = np.sort(np.linalg.eigvalsh(pair_block).real)
            same = np.allclose(full_values, pair_values)
            print(f"   L = {levels}, {pairs} pair(s): "
                  f"{len(rows):3d} states, agree: {same},  "
                  f"E_0 = {full_values[0]:.8f}")

    print()
    print("=" * 74)
    print("8. What the particle-hole term costs")
    print("=" * 74)
    print("The particle-hole term breaks pairs, so the seniority-zero space")
    print("is no longer invariant and the one-qubit-per-level encoding is not")
    print("available: the full 2L-qubit mapping is the only option.")
    print()
    levels = 2
    for f in (0.0, 0.1):
        H = pairing_ph_full_space(levels, g=1.0, f=f)
        decomposition = decompose(H)
        counts = {}
        for label in decomposition:
            counts[weight(label)] = counts.get(weight(label), 0) + 1
        keep = seniority_zero_indices(levels)
        occupancy = particle_number(2 * levels)
        rows = [k for k in keep if occupancy[k] == 2]
        others = [k for k in range(1 << (2 * levels))
                  if occupancy[k] == 2 and k not in rows]
        leakage = np.abs(H[np.ix_(others, rows)]).max()
        print(f"   f = {f}: {len(decomposition):3d} Pauli strings, "
              f"weights {dict(sorted(counts.items()))}")
        print(f"           largest matrix element out of the "
              f"seniority-zero space: {leakage:.4f}")
    print()
    print("At f = 0 the seniority-zero block is decoupled and the compact")
    print("encoding is exact.  As soon as f is switched on, the Hamiltonian")
    print("connects that block to the rest and the reduction fails.")

    print()
    print("=" * 74)
    print("9. How many Pauli strings?")
    print("=" * 74)
    print(f"{'model':>22s} {'qubits':>8s} {'strings':>9s} {'max weight':>12s}")
    rows = []
    for levels in (2, 3):
        rows.append((f"pairing, pair qubits", levels,
                     pairing_pair_qubits(levels, g=1.0)))
        rows.append((f"pairing, full space", 2 * levels,
                     pairing_full_space(levels, g=1.0)))
        rows.append((f"pairing + p-h, full", 2 * levels,
                     pairing_ph_full_space(levels, g=1.0, f=0.1)))
    for name, n_qubits, matrix in rows:
        decomposition = decompose(matrix)
        heaviest = max(weight(label) for label in decomposition)
        print(f"{name:>22s} {n_qubits:8d} {len(decomposition):9d} "
              f"{heaviest:12d}")
    print()
    print("The compact encoding wins twice over: half the qubits, and Pauli")
    print("strings of weight two instead of four.  Both matter for the")
    print("variational quantum eigensolver, where every string is a separate")
    print("measurement and every unit of weight is extra circuit depth.")

    print()
    print("=" * 74)
    print("10. Cross-check: the same numbers as the rest of the book")
    print("=" * 74)
    print("The qubit Hamiltonian is the same operator written differently, so")
    print("it must give the energies of chapter 4 exactly.  L = 4, N = 4,")
    print("g = 1, in the S_z = 0 sector:")
    print()
    print(f"{'f':>6s} {'dimension':>11s} {'E_0 (Jordan-Wigner)':>21s} "
          f"{'chapter 4':>14s}")
    reference = {0.0: 0.63554847, 0.05: 0.45058234,
                 0.20: -0.18348455, 0.50: -1.69173670}
    for f in (0.0, 0.05, 0.20, 0.50):
        H = pairing_ph_full_space(4, g=1.0, f=f)
        values, dimension = sector_spectrum(H, levels=4, particles=4, sz2=0)
        print(f"{f:6.2f} {dimension:11d} {values[0]:21.8f} "
              f"{reference[f]:14.8f}")
    print()
    print("Those are the entries of table 4.3, and the f = 0 row is the g = 1")
    print("row of table 4.1 -- the same number that fci.py, rpa.py, mbpt.py")
    print("and coupledcluster.py all reproduce.  Nothing has been lost in")
    print("translation to qubits.")

    print()
    print("=" * 74)
    print("11. Exercise 6 of chapter 3: a Hamiltonian on two qubits")
    print("=" * 74)
    eps0, eps1, U, t = 0.5, 1.5, 2.0, 0.3
    H, predicted = two_orbital_model(eps0, eps1, U, t)
    computed = decompose(H)
    print(f"eps_0 = {eps0}, eps_1 = {eps1}, U = {U}, t = {t}")
    print("Pauli decomposition:")
    show(computed)
    agree = (set(predicted) == set(computed)
             and all(abs(predicted[k] - computed[k]) < 1e-10
                     for k in predicted))
    print(f"   agrees with the hand calculation: {agree}")
    print()
    print("The matrix in the basis |00>, |01>, |10>, |11>:")
    for row in H.real:
        print("   " + "  ".join(f"{value:7.3f}" for value in row))
    print()
    one_particle = np.array([[eps1, t], [t, eps0]])
    exact = np.sort(np.linalg.eigvalsh(one_particle))
    closed = np.sort([0.5 * (eps0 + eps1)
                      + sign * np.sqrt((0.5 * (eps0 - eps1)) ** 2 + t ** 2)
                      for sign in (-1, 1)])
    print(f"   one-particle block eigenvalues : {exact[0]:.8f}, {exact[1]:.8f}")
    print(f"   closed form E_+- of the text   : {closed[0]:.8f}, "
          f"{closed[1]:.8f}")
    print(f"   agree: {bool(np.allclose(exact, closed))}")


if __name__ == "__main__":
    _demo()
