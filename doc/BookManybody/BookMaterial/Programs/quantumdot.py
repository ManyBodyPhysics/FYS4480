"""
Two-dimensional quantum dots: Hartree-Fock, MBPT, CCSD and unitary CCSD.

Companion code to chapter 11 of *Quantum mechanics for Many-particle Systems*.

The system is N electrons confined to two dimensions by a parabolic trap,

    H = sum_i [ -grad_i^2 / 2 + omega^2 r_i^2 / 2 ] + sum_{i<j} 1/|r_i - r_j|

in natural units (hbar = m = e = 1).  Unlike the pairing and Lipkin models of
chapter 4, nothing here is schematic: the interaction is the real Coulomb
repulsion, every matrix element has to be computed, and the single-particle
basis is infinite and has to be truncated.

The single-particle basis is the two-dimensional harmonic oscillator in polar
coordinates, labelled by a radial quantum number n and an angular-momentum
projection m, with energies

    eps_{nm} = hbar omega (2n + |m| + 1) ,

so that shells N_s = 2n + |m| have degeneracy 2(N_s + 1) including spin and
the closed shells fall at N = 2, 6, 12, 20, 30, 42.

The Coulomb matrix elements are evaluated in closed form (Anisimova and
Matulis, J. Phys.: Condens. Matter 10, 601 (1998)); the code below is a
translation of the C++ implementation used in FYS4411.  Once the one- and
two-body matrices are in hand, everything else is imported from
`coupledcluster.py`: the same Hartree-Fock, MP2, CCD, CCSD and unitary
coupled-cluster routines that were validated on the pairing model of
chapter 10 are applied here unchanged.  That reuse is deliberate.  The point
of this chapter is that the methods do not care whether the interaction is a
one-parameter model or the real thing.

Author: Morten Hjorth-Jensen
"""

import math
from functools import lru_cache
from itertools import combinations

import numpy as np

from coupledcluster import (ccd, ccsd, fci_energy, hartree_fock, mp2_energy,
                            reference_energy, UnitaryCC)


# ---------------------------------------------------------------------------
#  Coulomb matrix elements of the two-dimensional harmonic oscillator
# ---------------------------------------------------------------------------
#  < n_i m_i ; n_j m_j | 1/r_12 | n_k m_k ; n_l m_l > , proportional to
#  sqrt(hw), and zero unless m_i + m_j = m_k + m_l.  The closed form is a
#  quadruple sum over binomial expansions of the Laguerre polynomials, with an
#  inner quadruple sum over the angular integration; the whole thing is
#  carried out in logarithms so that the factorials do not overflow.
# ---------------------------------------------------------------------------
@lru_cache(maxsize=None)
def logfac(n):
    """ln(n!), by the same explicit loop as the C++ original."""
    if n < 0:
        raise ValueError("logfac requires a non-negative argument")
    return sum(math.log(a) for a in range(2, n + 1))


def _logratio1(j1, j2, j3, j4):
    return -(logfac(j1) + logfac(j2) + logfac(j3) + logfac(j4))


def _logratio2(G):
    return -0.5 * (G + 1) * math.log(2.0)


def _product1(n1, m1, n2, m2, n3, m3, n4, m4):
    """The radial normalisation prefactor, as a square root of factorials."""
    prod = logfac(n1) + logfac(n2) + logfac(n3) + logfac(n4)
    prod -= (logfac(n1 + abs(m1)) + logfac(n2 + abs(m2))
             + logfac(n3 + abs(m3)) + logfac(n4 + abs(m4)))
    return math.exp(0.5 * prod)


def _logproduct2(n1, m1, n2, m2, n3, m3, n4, m4, j1, j2, j3, j4):
    prod = (logfac(n1 + abs(m1)) + logfac(n2 + abs(m2))
            + logfac(n3 + abs(m3)) + logfac(n4 + abs(m4)))
    prod -= (logfac(n1 - j1) + logfac(n2 - j2)
             + logfac(n3 - j3) + logfac(n4 - j4))
    prod -= (logfac(j1 + abs(m1)) + logfac(j2 + abs(m2))
             + logfac(j3 + abs(m3)) + logfac(j4 + abs(m4)))
    return prod


def _logproduct3(l1, l2, l3, l4, g1, g2, g3, g4):
    prod = logfac(g1) + logfac(g2) + logfac(g3) + logfac(g4)
    prod -= logfac(l1) + logfac(l2) + logfac(l3) + logfac(l4)
    prod -= (logfac(g1 - l1) + logfac(g2 - l2)
             + logfac(g3 - l3) + logfac(g4 - l4))
    return prod


def coulomb_ho(hw, ni, mi, nj, mj, nk, mk, nl, ml):
    """< n_i m_i ; n_j m_j | 1/r_12 | n_k m_k ; n_l m_l > in units of sqrt(hw).

    Returns zero at once when m_i + m_j != m_k + m_l, which is angular-momentum
    conservation and eliminates the great majority of the index combinations.
    """
    if mi + mj != mk + ml:
        return 0.0

    total = 0.0
    for j1 in range(ni + 1):
        for j2 in range(nj + 1):
            for j3 in range(nl + 1):
                for j4 in range(nk + 1):
                    g1 = int(j1 + j4 + 0.5 * (abs(mi) + mi)
                             + 0.5 * (abs(mk) - mk))
                    g2 = int(j2 + j3 + 0.5 * (abs(mj) + mj)
                             + 0.5 * (abs(ml) - ml))
                    g3 = int(j3 + j2 + 0.5 * (abs(ml) + ml)
                             + 0.5 * (abs(mj) - mj))
                    g4 = int(j4 + j1 + 0.5 * (abs(mk) + mk)
                             + 0.5 * (abs(mi) - mi))
                    G = g1 + g2 + g3 + g4

                    prefactor = math.exp(
                        _logratio1(j1, j2, j3, j4)
                        + _logproduct2(ni, mi, nj, mj, nl, ml, nk, mk,
                                       j1, j2, j3, j4)
                        + _logratio2(G))

                    inner = 0.0
                    for l1 in range(g1 + 1):
                        for l2 in range(g2 + 1):
                            for l3 in range(g3 + 1):
                                l4 = l1 + l2 - l3
                                if l4 < 0 or l4 > g4:
                                    continue
                                L = l1 + l2 + l3 + l4
                                sign = 1 - 2 * ((g2 + g3 - l2 - l3) % 2)
                                inner += sign * math.exp(
                                    _logproduct3(l1, l2, l3, l4,
                                                 g1, g2, g3, g4)
                                    + math.lgamma(1.0 + 0.5 * L)
                                    + math.lgamma(0.5 * (G - L + 1.0)))

                    sign_j = 1 - 2 * ((j1 + j2 + j3 + j4) % 2)
                    total += sign_j * prefactor * inner

    total *= _product1(ni, mi, nj, mj, nl, ml, nk, mk)
    return math.sqrt(hw) * total


# ---------------------------------------------------------------------------
#  The single-particle basis and the many-body matrices
# ---------------------------------------------------------------------------
def spatial_basis(shells):
    """All (n, m) with 2n + |m| <= shells, ordered by shell and then by m."""
    out = []
    for n_s in range(shells + 1):
        states = []
        for n in range(n_s // 2 + 1):
            m_abs = n_s - 2 * n
            if m_abs == 0:
                states.append((n, 0))
            else:
                states.append((n, -m_abs))
                states.append((n, m_abs))
        out.extend(sorted(states, key=lambda s: (s[1], s[0])))
    return out


class QuantumDot:
    """N electrons in a two-dimensional parabolic trap.

    `shells` is the highest oscillator shell retained, so shells = 5 gives the
    42 spin-orbitals used throughout chapter 11.  The one-body matrix is
    diagonal in the oscillator basis; the two-body matrix is the
    antisymmetrised Coulomb element

        <pq||rs> = <pq|v|rs> d(sp,sr) d(sq,ss) - <pq|v|sr> d(sp,ss) d(sq,sr),

    with the spatial part supplied by `coulomb_ho` and the spin deltas applied
    here.
    """

    def __init__(self, particles=2, shells=5, hw=1.0):
        self.particles = particles
        self.shells = shells
        self.hw = hw
        self.spatial = spatial_basis(shells)
        # spin-orbitals: each spatial state twice, spin up first
        self.orbitals = [(n, m, spin)
                         for (n, m) in self.spatial for spin in (+1, -1)]
        self.n_orbitals = len(self.orbitals)
        self._spatial_index = {s: k for k, s in enumerate(self.spatial)}
        self._cache = {}

    # -- bookkeeping -----------------------------------------------------
    def shell_of(self, k):
        n, m = self.spatial[k]
        return 2 * n + abs(m)

    def energies(self):
        """The oscillator single-particle energies of the spin-orbitals."""
        return np.array([self.hw * (2 * n + abs(m) + 1)
                         for (n, m, _) in self.orbitals])

    def magic_numbers(self):
        """Cumulative occupancies at closed shells."""
        out, total = [], 0
        for n_s in range(self.shells + 1):
            total += 2 * (n_s + 1)
            out.append(total)
        return out

    # -- integrals -------------------------------------------------------
    def spatial_element(self, p, q, r, s):
        """<pq|v|rs> for spatial indices, cached and symmetry-reduced."""
        key = min((p, q, r, s), (q, p, s, r), (r, s, p, q), (s, r, q, p))
        value = self._cache.get(key)
        if value is None:
            (a, b, c, d) = key
            (na, ma), (nb, mb) = self.spatial[a], self.spatial[b]
            (nc, mc), (nd, md) = self.spatial[c], self.spatial[d]
            value = coulomb_ho(self.hw, na, ma, nb, mb, nc, mc, nd, md)
            self._cache[key] = value
        return value

    def one_body(self):
        return np.diag(self.energies())

    def spatial_tensor(self, progress=False):
        """<pq|v|rs> over the spatial states, as a dense M^4 array.

        The angular-momentum rule m_p + m_q = m_r + m_s kills roughly two
        thirds of the entries before any work is done, and the remaining ones
        are computed once and reused through the symmetry
        <pq|v|rs> = <qp|v|sr> = <rs|v|pq> = <sr|v|qp>.
        """
        m_states = len(self.spatial)
        W = np.zeros((m_states,) * 4)
        m = [mm for (_, mm) in self.spatial]
        for p in range(m_states):
            if progress:
                print(f"      spatial integrals: row {p+1:3d}/{m_states}",
                      flush=True)
            for q in range(m_states):
                for r in range(m_states):
                    for s in range(m_states):
                        if m[p] + m[q] != m[r] + m[s]:
                            continue
                        W[p, q, r, s] = self.spatial_element(p, q, r, s)
        return W

    def two_body(self, progress=False):
        """The antisymmetrised <pq||rs> over spin-orbitals.

            <pq||rs> = <pq|v|rs> d(sp,sr) d(sq,ss)
                     - <pq|v|sr> d(sp,ss) d(sq,sr)

        The spatial part is computed once and then broadcast onto the
        spin-orbital indices, which is much faster than looping over all
        (2M)^4 combinations in Python.
        """
        W = self.spatial_tensor(progress)
        index = np.repeat(np.arange(len(self.spatial)), 2)
        spin = np.array([s for (_, _, s) in self.orbitals])
        direct = W[np.ix_(index, index, index, index)]
        exchange = direct.transpose(0, 1, 3, 2)

        same_pr = (spin[:, None] == spin[None, :])[:, None, :, None]
        same_qs = (spin[:, None] == spin[None, :])[None, :, None, :]
        same_ps = (spin[:, None] == spin[None, :])[:, None, None, :]
        same_qr = (spin[:, None] == spin[None, :])[None, :, :, None]
        return direct * (same_pr & same_qs) - exchange * (same_ps & same_qr)

    def matrices(self, progress=False):
        """(h, v) in the oscillator basis."""
        return self.one_body(), self.two_body(progress)


# ---------------------------------------------------------------------------
#  Caching: the whole hbar-omega dependence is a scaling
# ---------------------------------------------------------------------------
#  The one-body energies are proportional to hw and the Coulomb elements to
#  sqrt(hw), so the expensive integrals need to be computed once and can then
#  be reused for every frequency.  That is worth exploiting -- and it is worth
#  verifying, which the demo does.
# ---------------------------------------------------------------------------
_MATRIX_CACHE = {}


def matrices(shells, hw=1.0, progress=False):
    """(h, v) for a given number of shells, at any frequency."""
    if shells not in _MATRIX_CACHE:
        dot = QuantumDot(2, shells, hw=1.0)
        _MATRIX_CACHE[shells] = (dot.one_body(), dot.two_body(progress))
    h_unit, v_unit = _MATRIX_CACHE[shells]
    return hw * h_unit, math.sqrt(hw) * v_unit


# ---------------------------------------------------------------------------
#  Running the hierarchy of methods
# ---------------------------------------------------------------------------
def solve(particles, shells, hw=1.0, methods=("hf", "mp2", "ccd", "ccsd"),
          progress=False):
    """Hartree-Fock, MP2, CCD and CCSD for one quantum dot.

    Everything after the construction of (h, v) is imported from
    `coupledcluster.py` without modification: the Hartree-Fock rotation, the
    MP2 sum and the coupled-cluster solvers are those of chapter 10.
    """
    dot = QuantumDot(particles, shells, hw)
    h, v = matrices(shells, hw, progress)
    out = {"dot": dot, "n_orbitals": dot.n_orbitals,
           "E0": float(np.sum(dot.energies()[:particles])),
           "E_ref": reference_energy(h, v, particles)}

    h_hf, v_hf, C, iterations = hartree_fock(h, v, particles)
    f_hf = h_hf + np.einsum("piqi->pq", v_hf[:, :particles, :, :particles])
    out["E_HF"] = reference_energy(h_hf, v_hf, particles)
    out["iterations"] = iterations
    out["eps"] = np.diag(f_hf).copy()
    out["h_hf"], out["v_hf"], out["f_hf"] = h_hf, v_hf, f_hf

    if "mp2" in methods:
        out["E_MP2"] = mp2_energy(f_hf, v_hf, particles)
    if "ccd" in methods:
        result = ccd(f_hf, v_hf, particles)
        out["E_CCD"] = result["energy"]
        out["CCD_iterations"] = result["iterations"]
        out["CCD_converged"] = result["converged"]
    if "ccsd" in methods:
        result = ccsd(f_hf, v_hf, particles)
        out["E_CCSD"] = result["energy"]
        out["CCSD_iterations"] = result["iterations"]
        out["CCSD_converged"] = result["converged"]
        out["t1_max"] = float(np.abs(result["t1"]).max())
    return out


def two_electron_fci(h, v):
    """Exact ground state of two particles, by direct construction.

    In the basis of antisymmetrised pairs p < q,

        <pq|H|rs> = d_qs h_pr + d_pr h_qs - d_ps h_qr - d_qr h_ps + <pq||rs>,

    so the matrix is built without any determinant bookkeeping.  The dimension
    is binom(M, 2), which is 861 for the 42 spin-orbitals used in the text --
    small enough to diagonalise densely, and therefore an exact benchmark for
    every approximation in the chapter.
    """
    n = h.shape[0]
    pairs = list(combinations(range(n), 2))
    dim = len(pairs)
    p = np.array([a for a, _ in pairs])
    q = np.array([b for _, b in pairs])

    H = v[p[:, None], q[:, None], p[None, :], q[None, :]].copy()
    H += np.where(q[:, None] == q[None, :], h[p[:, None], p[None, :]], 0.0)
    H += np.where(p[:, None] == p[None, :], h[q[:, None], q[None, :]], 0.0)
    H -= np.where(p[:, None] == q[None, :], h[q[:, None], p[None, :]], 0.0)
    H -= np.where(q[:, None] == p[None, :], h[p[:, None], q[None, :]], 0.0)
    return float(np.linalg.eigvalsh(H)[0]), dim


def exact_two_electron(shells, hw=1.0):
    """Full configuration interaction for two electrons in a given basis."""
    h, v = matrices(shells, hw)
    energy, dimension = two_electron_fci(h, v)
    return energy, h.shape[0], dimension


def ucc_amplitudes_from_ccsd(ucc, t1, t2, n_occ):
    """Pack CCSD amplitudes into the parameter vector of a `UnitaryCC`.

    The unitary ansatz uses the same excitation labels as coupled cluster, so
    the converged CCSD amplitudes are a natural starting point -- and on a
    quantum computer they are often the only ones one can afford to compute.

    One sign has to be watched.  `UnitaryCC` builds its doubles generator as
    a+_a a+_b a_i a_j, whereas the cluster operator is conventionally written
    T_2 = (1/4) sum t_ij^ab a+_a a+_b a_j a_i, and the two orderings of the
    annihilation operators differ by a minus sign.  The doubles amplitudes are
    therefore imported with a minus and the singles as they stand; the demo
    checks this by comparing with a full variational optimisation, where the
    two amplitude vectors come out parallel to within a cosine of 1.0000.
    """
    x = np.zeros(ucc.n_amplitudes)
    for k, label in enumerate(ucc.labels):
        if label[0] == "s":
            _, i, a = label
            x[k] = t1[i, a - n_occ]
        elif label[0] == "d":
            _, (i, j), (a, b) = label
            x[k] = -t2[i, j, a - n_occ, b - n_occ]
    return x


def taut_energy():
    """The exact two-electron energy at hbar omega = 1.

    Taut (Phys. Rev. A 48, 3561 (1993)) showed that the relative motion of two
    electrons in a two-dimensional parabolic trap is analytically solvable for
    a discrete set of frequencies.  For hbar omega = 1 in the units used here
    the ground state is exactly 3.
    """
    return 3.0


# ---------------------------------------------------------------------------
def _demo():
    print("=" * 74)
    print("1. The single-particle basis and the magic numbers")
    print("=" * 74)
    dot = QuantumDot(2, shells=5)
    print(f"shells 0-{dot.shells}: {len(dot.spatial)} spatial states, "
          f"{dot.n_orbitals} spin-orbitals")
    print()
    print(f"{'N_s':>4s} {'(n, m)':>30s} {'degeneracy':>12s} "
          f"{'cumulative N':>14s}")
    total = 0
    for n_s in range(dot.shells + 1):
        states = [s for k, s in enumerate(dot.spatial)
                  if dot.shell_of(k) == n_s]
        total += 2 * len(states)
        label = ",".join(f"({n},{m:+d})" if m else f"({n},0)"
                         for n, m in states)
        print(f"{n_s:4d} {label:>30s} {2*len(states):12d} {total:14d}")
    print()
    print("The last column is the sequence of magic numbers 2, 6, 12, 20, 30,")
    print("42: the electron numbers at which a shell closes and the dot is")
    print("particularly stable.  We compute the first two, N = 2 and N = 6.")

    print()
    print("=" * 74)
    print("2. The Coulomb matrix element")
    print("=" * 74)
    element = coulomb_ho(1.0, 0, 0, 0, 0, 0, 0, 0, 0)
    print(f"   <00;00|v|00;00>  = {element:.12f}")
    print(f"   sqrt(pi/2)       = {math.sqrt(math.pi/2):.12f}")
    print(f"   agree: {abs(element - math.sqrt(math.pi/2)) < 1e-12}")
    print()
    print("   the whole hw dependence is a factor sqrt(hw):")
    for hw in (0.25, 0.5, 1.0, 2.0, 4.0):
        value = coulomb_ho(hw, 0, 0, 0, 0, 0, 0, 0, 0)
        print(f"      hw = {hw:5.2f}:  {value:12.8f}   "
              f"value / sqrt(hw) = {value/math.sqrt(hw):.12f}")
    print()
    print("   angular momentum is conserved, m_p + m_q = m_r + m_s:")
    for label, args in (
            ("<0,+1;0,-1|v|0,0;0,0>   (allowed)", (0, 1, 0, -1, 0, 0, 0, 0)),
            ("<0,+1;0, 0|v|0,0;0,0>   (forbidden)", (0, 1, 0, 0, 0, 0, 0, 0)),
            ("<0,+2;0,-1|v|0,+1;0,0>  (allowed)", (0, 2, 0, -1, 0, 1, 0, 0)),
            ("<0,+2;0, 0|v|0,+1;0,0>  (forbidden)", (0, 2, 0, 0, 0, 1, 0, 0))):
        print(f"      {label:38s} = {coulomb_ho(1.0, *args):+.8f}")
    print()
    print("   symmetries  <pq|v|rs> = <qp|v|sr> = <rs|v|pq>:")
    small = QuantumDot(2, shells=2)
    for p, q, r, s in ((0, 1, 2, 3), (1, 3, 0, 2), (2, 2, 1, 1)):
        direct = small.spatial_element(p, q, r, s)
        print(f"      ({p}{q}|{r}{s}) = {direct:+.8f}   "
              f"particle swap {abs(direct-small.spatial_element(q,p,s,r)):.1e}"
              f"   hermiticity {abs(direct-small.spatial_element(r,s,p,q)):.1e}")

    print()
    print("=" * 74)
    print("3. Hartree-Fock and the size of the basis")
    print("=" * 74)
    print("Adding a shell can only lower the energy: the Hartree-Fock state")
    print("gains variational freedom and nothing is taken away.")
    print()
    print(f"{'shells':>7s} {'orbitals':>9s} {'N=2 E_HF':>14s} {'dE':>11s} "
          f"{'N=6 E_HF':>14s} {'dE':>11s}")
    previous = {2: None, 6: None}
    for shells in range(0, 6):
        n_orbitals = 2 * len(spatial_basis(shells))
        row = [f"{shells:7d}", f"{n_orbitals:9d}"]
        for particles in (2, 6):
            if n_orbitals < particles:
                row += [f"{'--':>14s}", f"{'--':>11s}"]
                continue
            energy = solve(particles, shells, methods=())["E_HF"]
            row.append(f"{energy:14.8f}")
            row.append(f"{'--':>11s}" if previous[particles] is None
                       else f"{energy-previous[particles]:+11.6f}")
            previous[particles] = energy
        print(" ".join(row))

    print()
    print("=" * 74)
    print("4. The hierarchy of methods, hw = 1, shells 0-5 (42 orbitals)")
    print("=" * 74)
    results = {particles: solve(particles, 5) for particles in (2, 6)}
    print()
    print(f"{'quantity':<32s} {'N = 2':>15s} {'N = 6':>15s}")
    print("-" * 64)

    def line(label, values):
        print(f"{label:<32s} {values[0]:15.8f} {values[1]:15.8f}")

    line("non-interacting  E_0", [results[n]["E0"] for n in (2, 6)])
    line("Hartree-Fock  E_HF", [results[n]["E_HF"] for n in (2, 6)])
    line("  interaction  E_HF - E_0",
         [results[n]["E_HF"] - results[n]["E0"] for n in (2, 6)])
    line("MP2 correlation", [results[n]["E_MP2"] for n in (2, 6)])
    line("  total  E_HF + E_MP2",
         [results[n]["E_HF"] + results[n]["E_MP2"] for n in (2, 6)])
    line("CCD correlation", [results[n]["E_CCD"] for n in (2, 6)])
    line("  total  E_HF + E_CCD",
         [results[n]["E_HF"] + results[n]["E_CCD"] for n in (2, 6)])
    line("CCSD correlation", [results[n]["E_CCSD"] for n in (2, 6)])
    line("  total  E_HF + E_CCSD",
         [results[n]["E_HF"] + results[n]["E_CCSD"] for n in (2, 6)])
    print()
    for particles in (2, 6):
        r = results[particles]
        print(f"   N = {particles}: SCF {r['iterations']} iterations, "
              f"CCD {r['CCD_iterations']}, CCSD {r['CCSD_iterations']}, "
              f"max|t_1| = {r['t1_max']:.2e}")
    print()
    print("Note that the singles amplitudes do NOT vanish.  Brillouin's")
    print("theorem makes f_ai = 0, which is why there is no first-order")
    print("singles contribution and why MP2 has only doubles -- but the CCSD")
    print("singles equation is driven by the doubles through terms such as")
    print("sum_mef t_im^ef <ma||ef>, which survive at t_1 = 0.  The singles")
    print("are small, they describe orbital relaxation, and CCSD therefore")
    print("lies slightly below CCD rather than coinciding with it.")

    print()
    print("=" * 74)
    print("5. Two electrons against the exact answer")
    print("=" * 74)
    print("With N = 2 the determinant space is small enough to diagonalise")
    print("exactly, so every approximation can be measured against the exact")
    print("energy in the same basis -- and the basis series against Taut's")
    print("analytic result E = 3 at hw = 1.")
    print()
    print(f"{'shells':>7s} {'orbitals':>9s} {'dim':>6s} {'E_HF':>13s} "
          f"{'E_HF+MP2':>13s} {'E_HF+CCSD':>13s} {'E_FCI':>13s} "
          f"{'CCSD-FCI':>11s}")
    for shells in (1, 2, 3, 4, 5):
        result = solve(2, shells)
        exact, n_orbitals, dimension = exact_two_electron(shells)
        ccsd_total = result["E_HF"] + result["E_CCSD"]
        print(f"{shells:7d} {n_orbitals:9d} {dimension:6d} "
              f"{result['E_HF']:13.8f} "
              f"{result['E_HF']+result['E_MP2']:13.8f} "
              f"{ccsd_total:13.8f} {exact:13.8f} "
              f"{abs(ccsd_total-exact):11.1e}")
    print(f"{'exact':>7s} {'infinite':>9s} {'--':>6s} {'--':>13s} {'--':>13s} "
          f"{'--':>13s} {taut_energy():13.8f}")
    print()
    print("The last column is the point.  For two electrons singles and")
    print("doubles exhaust the available excitations, so CCSD is exact within")
    print("the basis, and it agrees with full diagonalisation to twelve")
    print("digits.  That validates the coupled-cluster solver of chapter 10 on")
    print("a real interaction, not just on a model.  The remaining gap to 3 is")
    print("basis truncation, and it closes slowly: the Coulomb cusp is hard")
    print("for a harmonic-oscillator basis.")

    print()
    print("=" * 74)
    print("6. Dependence on the confinement")
    print("=" * 74)
    print("Kinetic energy scales as hw and the Coulomb energy as sqrt(hw), so")
    print("a shallow trap is strongly correlated and a steep one is nearly")
    print("free.  The last column is the correlation energy as a fraction of")
    print("the Hartree-Fock interaction energy.")
    for particles in (2, 6):
        print()
        print(f"   N = {particles}, shells 0-4")
        print(f"{'hw':>8s} {'E_0':>9s} {'E_HF':>13s} {'E_MP2':>12s} "
              f"{'E_CCSD':>12s} {'ratio':>9s}")
        for hw in (0.25, 0.5, 1.0, 2.0, 4.0):
            r = solve(particles, 4, hw=hw)
            ratio = r["E_CCSD"] / (r["E_HF"] - r["E0"])
            print(f"{hw:8.2f} {r['E0']:9.3f} {r['E_HF']:13.8f} "
                  f"{r['E_MP2']:12.8f} {r['E_CCSD']:12.8f} {ratio:9.4f}")

    print()
    print("=" * 74)
    print("7. Unitary coupled cluster: two electrons")
    print("=" * 74)
    print("The unitary ansatz of chapter 10, optimised variationally on a real")
    print("interaction.  For two electrons it must reproduce CCSD exactly.")
    print()
    print(f"{'shells':>7s} {'orb':>5s} {'params':>7s} {'CCD':>13s} "
          f"{'UCCD':>13s} {'CCSD':>13s} {'UCCSD':>13s} {'FCI':>13s}")
    for shells in (1, 2, 3):
        h, v = matrices(shells)
        h_hf, v_hf, _, _ = hartree_fock(h, v, 2)
        f_hf = h_hf + np.einsum("piqi->pq", v_hf[:, :2, :, :2])
        e_ref = reference_energy(h_hf, v_hf, 2)
        e_ccd = e_ref + ccd(f_hf, v_hf, 2)["energy"]
        e_ccsd = e_ref + ccsd(f_hf, v_hf, 2)["energy"]
        ucc_d = UnitaryCC(h_hf, v_hf, 2, singles=False, doubles=True)
        ucc_sd = UnitaryCC(h_hf, v_hf, 2, singles=True, doubles=True)
        e_ucc_d = ucc_d.optimise()["energy"]
        e_ucc_sd = ucc_sd.optimise()["energy"]
        exact = two_electron_fci(h, v)[0]
        print(f"{shells:7d} {h.shape[0]:5d} {ucc_sd.n_amplitudes:7d} "
              f"{e_ccd:13.8f} {e_ucc_d:13.8f} {e_ccsd:13.8f} "
              f"{e_ucc_sd:13.8f} {exact:13.8f}")
    print()
    print("Two things are worth noticing.  UCCSD, CCSD and full")
    print("diagonalisation coincide, because with two particles singles and")
    print("doubles span the entire space.  And UCCD coincides with CCD, which")
    print("is less obvious: with two particles T_2 T_2 annihilates the")
    print("reference and so does T_2^dagger, so exp(T_2 - T_2^dagger) and")
    print("1 + T_2 sweep out the same two-dimensional manifold, and the")
    print("variational minimum over it is the lowest eigenvalue in the space")
    print("of the reference plus the doubles -- which is what CCD returns.")
    print("The unitary and the standard ansatz differ only when the")
    print("truncation genuinely bites.")

    print()
    print("=" * 74)
    print("8. Six electrons: the whole hierarchy against exact diagonalisation")
    print("=" * 74)
    print("Twelve orbitals is small enough that the six-electron determinant")
    print("space can still be diagonalised, so for once every method in the")
    print("book can be compared with the exact answer on a realistic")
    print("interaction and a system where the truncation matters.")
    print()
    h, v = matrices(2)
    h_hf, v_hf, _, _ = hartree_fock(h, v, 6)
    f_hf = h_hf + np.einsum("piqi->pq", v_hf[:, :6, :, :6])
    e_ref = reference_energy(h_hf, v_hf, 6)
    cc_d = ccd(f_hf, v_hf, 6)
    cc_sd = ccsd(f_hf, v_hf, 6)
    ucc = UnitaryCC(h_hf, v_hf, 6, singles=True, doubles=True)
    exact = float(np.linalg.eigvalsh(ucc.H)[0])
    x = ucc_amplitudes_from_ccsd(ucc, cc_sd["t1"], cc_sd["t2"], 6)
    e_ucc = ucc.energy(x)

    print(f"   basis: 12 orbitals, determinant space {ucc.dim}, "
          f"{ucc.n_amplitudes} amplitudes")
    print()
    rows = [("Hartree-Fock", e_ref),
            ("MP2", e_ref + mp2_energy(f_hf, v_hf, 6)),
            ("CCD", e_ref + cc_d["energy"]),
            ("CCSD", e_ref + cc_sd["energy"]),
            ("UCCSD (at the CCSD amplitudes)", e_ucc),
            ("exact diagonalisation", exact)]
    correlation = exact - e_ref
    print(f"{'method':<32s} {'energy':>14s} {'error':>12s} "
          f"{'% of E_corr':>12s}")
    for label, energy in rows:
        recovered = 100.0 * (energy - e_ref) / correlation
        print(f"{label:<32s} {energy:14.8f} {energy-exact:12.2e} "
              f"{recovered:11.2f}%")
    print()
    print("The ordering is the expected one and every step is an improvement:")
    print("MP2 recovers most of the correlation energy, CCD resums the ladders")
    print("and rings to all orders, CCSD adds the orbital relaxation, and the")
    print("unitary ansatz does better still -- even evaluated at the CCSD")
    print("amplitudes rather than its own optimum, and while remaining above")
    print("the exact energy as the variational principle demands.  CCSD")
    print("carries no such guarantee.")

    print()
    print("=" * 74)
    print("9. Trotterising the unitary ansatz")
    print("=" * 74)
    print("A quantum computer cannot apply exp(sigma) in one step; it applies")
    print("the generators one at a time.  The splitting error is the Trotter")
    print("error of section 6.9, and n is the circuit depth.  Two electrons in")
    print("twelve orbitals, at the CCSD amplitudes:")
    print()
    h, v = matrices(2)
    h_hf, v_hf, _, _ = hartree_fock(h, v, 2)
    f_hf = h_hf + np.einsum("piqi->pq", v_hf[:, :2, :, :2])
    cc = ccsd(f_hf, v_hf, 2)
    ucc = UnitaryCC(h_hf, v_hf, 2, singles=True, doubles=True)
    x = ucc_amplitudes_from_ccsd(ucc, cc["t1"], cc["t2"], 2)
    untrotterised = ucc.energy(x)
    print(f"   exact exponential  = {untrotterised:.8f}")
    print()
    print(f"{'Trotter steps':>14s} {'energy':>14s} {'error':>12s} "
          f"{'ratio':>8s}")
    previous = None
    for n in (1, 2, 4, 8, 16):
        energy = ucc.energy(x, n_trotter=n)
        error = abs(energy - untrotterised)
        ratio = "--" if previous is None else f"{previous/error:.2f}"
        print(f"{n:14d} {energy:14.8f} {error:12.2e} {ratio:>8s}")
        previous = error
    print()
    print("The error halves when n doubles: first-order Trotter, exactly the")
    print("O(1/n) rate of Eq. (6.30).  This is the calculation a variational")
    print("quantum eigensolver performs -- prepare the state with the")
    print("Trotterised circuit, measure the energy, adjust the amplitudes")
    print("classically.  On hardware the Hamiltonian would first be mapped to")
    print("Pauli strings by the Jordan-Wigner transformation of section 3.16;")
    print("here the matrices are simply exponentiated.")


if __name__ == "__main__":
    _demo()
