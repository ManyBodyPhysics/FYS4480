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
    h, v = dot.matrices(progress)
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
    dot = QuantumDot(2, shells, hw)
    h, v = dot.matrices()
    energy, dimension = two_electron_fci(h, v)
    return energy, dot.n_orbitals, dimension


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
    print("1. The single-particle basis")
    print("=" * 74)
    dot = QuantumDot(2, shells=5)
    print(f"shells 0-{dot.shells}: {len(dot.spatial)} spatial states, "
          f"{dot.n_orbitals} spin-orbitals")
    print()
    print(f"{'N_s':>4s} {'(n, m)':>28s} {'degeneracy':>12s} "
          f"{'cumulative N':>14s}")
    total = 0
    for n_s in range(dot.shells + 1):
        states = [s for k, s in enumerate(dot.spatial)
                  if dot.shell_of(k) == n_s]
        total += 2 * len(states)
        label = ", ".join(f"({n},{m:+d})" if m else f"({n},0)"
                          for n, m in states)
        print(f"{n_s:4d} {label:>28s} {2*len(states):12d} {total:14d}")
    print()
    print("The cumulative column is the sequence of magic numbers,")
    print("2, 6, 12, 20, 30, 42: the closed-shell electron numbers.")

    print()
    print("=" * 74)
    print("2. The Coulomb matrix element against known values")
    print("=" * 74)
    element = coulomb_ho(1.0, 0, 0, 0, 0, 0, 0, 0, 0)
    print(f"   <00;00|v|00;00>        = {element:.10f}")
    print(f"   sqrt(pi/2)             = {math.sqrt(math.pi/2):.10f}")
    print(f"   agree: {abs(element - math.sqrt(math.pi/2)) < 1e-12}")
    print()
    print("   scaling in hw (the element must go as sqrt(hw)):")
    for hw in (0.25, 0.5, 1.0, 2.0, 4.0):
        value = coulomb_ho(hw, 0, 0, 0, 0, 0, 0, 0, 0)
        print(f"      hw = {hw:5.2f}:  {value:.8f}   "
              f"ratio to sqrt(hw) = {value/math.sqrt(hw):.8f}")
    print()
    print("   angular-momentum conservation, m_p + m_q = m_r + m_s:")
    for label, args in (
            ("<0,1;0,-1|v|0,0;0,0>  (allowed)", (0, 1, 0, -1, 0, 0, 0, 0)),
            ("<0,1;0,0|v|0,0;0,0>   (forbidden)", (0, 1, 0, 0, 0, 0, 0, 0)),
            ("<0,2;0,-1|v|0,1;0,0>  (allowed)", (0, 2, 0, -1, 0, 1, 0, 0))):
        print(f"      {label:36s} = {coulomb_ho(1.0, *args):+.8f}")
    print()
    print("   symmetries of the spatial element:")
    dot = QuantumDot(2, shells=2)
    quads = [(0, 1, 2, 3), (1, 3, 0, 2), (2, 2, 1, 1)]
    for p, q, r, s in quads:
        direct = dot.spatial_element(p, q, r, s)
        swapped = dot.spatial_element(q, p, s, r)
        conjugate = dot.spatial_element(r, s, p, q)
        print(f"      ({p}{q}|{r}{s}): {direct:+.8f}   "
              f"particle swap {abs(direct-swapped):.1e}   "
              f"hermiticity {abs(direct-conjugate):.1e}")

    print()
    print("=" * 74)
    print("3. Hartree-Fock: convergence with the size of the basis")
    print("=" * 74)
    print("Adding shells can only lower the energy, since the Hartree-Fock")
    print("state gains variational freedom.")
    print()
    print(f"{'shells':>7s} {'orbitals':>9s} {'N=2 E_HF':>14s} {'dE':>11s} "
          f"{'N=6 E_HF':>14s} {'dE':>11s}")
    previous2 = previous6 = None
    for shells in range(0, 6):
        row = [f"{shells:7d}", f"{2*len(spatial_basis(shells)):9d}"]
        for particles, previous in ((2, previous2), (6, previous6)):
            if 2 * len(spatial_basis(shells)) < particles:
                row += [f"{'--':>14s}", f"{'--':>11s}"]
                continue
            result = solve(particles, shells, methods=())
            energy = result["E_HF"]
            row.append(f"{energy:14.8f}")
            row.append(f"{'--':>11s}" if previous is None
                       else f"{energy-previous:+11.6f}")
            if particles == 2:
                previous2 = energy
            else:
                previous6 = energy
        print(" ".join(row))

    print()
    print("=" * 74)
    print("4. The hierarchy of methods at hw = 1, shells 0-5")
    print("=" * 74)
    results = {}
    for particles in (2, 6):
        print(f"   solving N = {particles} ...", flush=True)
        results[particles] = solve(particles, 5)
    print()
    print(f"{'quantity':<34s} {'N = 2':>16s} {'N = 6':>16s}")
    print("-" * 68)
    rows = [
        ("non-interacting  E_0", "E0"),
        ("Hartree-Fock  E_HF", "E_HF"),
        ("  interaction  E_HF - E_0", None),
        ("MP2 correlation  E_MP2", "E_MP2"),
        ("  total  E_HF + E_MP2", None),
        ("CCD correlation  E_CCD", "E_CCD"),
        ("  total  E_HF + E_CCD", None),
        ("CCSD correlation  E_CCSD", "E_CCSD"),
        ("  total  E_HF + E_CCSD", None),
    ]
    for label, key in rows:
        values = []
        for particles in (2, 6):
            r = results[particles]
            if key is not None:
                values.append(r[key])
            elif "interaction" in label:
                values.append(r["E_HF"] - r["E0"])
            elif "MP2" in label:
                values.append(r["E_HF"] + r["E_MP2"])
            elif "CCD" in label:
                values.append(r["E_HF"] + r["E_CCD"])
            else:
                values.append(r["E_HF"] + r["E_CCSD"])
        print(f"{label:<34s} {values[0]:16.8f} {values[1]:16.8f}")
    print()
    for particles in (2, 6):
        r = results[particles]
        print(f"   N = {particles}: SCF {r['iterations']} iterations, "
              f"CCD {r['CCD_iterations']}, CCSD {r['CCSD_iterations']}; "
              f"|E_CCSD - E_CCD| = {abs(r['E_CCSD']-r['E_CCD']):.2e}")
    print()
    print("The last number is the check that matters.  For a closed-shell")
    print("system in a canonical Hartree-Fock basis Brillouin's theorem makes")
    print("the singles amplitudes vanish, so CCSD must reduce to CCD exactly.")
    print("A difference larger than the convergence tolerance would mean one")
    print("of the two solvers is wrong.")

    print()
    print("=" * 74)
    print("5. Two electrons against the exact answer")
    print("=" * 74)
    print("With N = 2 the determinant space is small enough to diagonalise,")
    print("so every approximation can be measured against the exact energy of")
    print("the same basis, and the basis series against Taut's analytic")
    print("result E = 3 at hw = 1.")
    print()
    print(f"{'shells':>7s} {'orbitals':>9s} {'E_HF':>12s} {'E_MP2':>12s} "
          f"{'E_CCSD':>12s} {'E_FCI':>12s}")
    for shells in (1, 2, 3, 4):
        result = solve(2, shells)
        exact, n_orbitals = exact_two_electron(shells)
        print(f"{shells:7d} {n_orbitals:9d} {result['E_HF']:12.8f} "
              f"{result['E_HF']+result['E_MP2']:12.8f} "
              f"{result['E_HF']+result['E_CCSD']:12.8f} {exact:12.8f}")
    print(f"{'exact':>7s} {'infinite':>9s} {'':12s} {'':12s} {'':12s} "
          f"{taut_energy():12.8f}")
    print()
    print("For two electrons CCSD is exact within the given basis -- singles")
    print("and doubles exhaust the excitations available to two particles --")
    print("so the CCSD and FCI columns must agree to machine precision.  They")
    print("do, which validates the coupled-cluster solver on a real")
    print("interaction.  The remaining gap to 3 is basis truncation, and it")
    print("closes slowly: the Coulomb cusp is hard for an oscillator basis.")

    print()
    print("=" * 74)
    print("6. Dependence on the confinement strength")
    print("=" * 74)
    print("Kinetic energy scales as hw and the Coulomb energy as sqrt(hw), so")
    print("a shallow trap is strongly correlated and a steep one is nearly")
    print("free.  The ratio E_corr/(E_HF - E_0) measures that directly.")
    print()
    for particles in (2, 6):
        print(f"   N = {particles}")
        print(f"{'hw':>8s} {'E_0':>10s} {'E_HF':>13s} {'E_MP2':>12s} "
              f"{'E_CCSD':>12s} {'ratio':>9s}")
        for hw in (0.28, 0.5, 1.0, 2.0, 4.0):
            r = solve(particles, 4, hw=hw)
            ratio = r["E_CCSD"] / (r["E_HF"] - r["E0"])
            print(f"{hw:8.2f} {r['E0']:10.4f} {r['E_HF']:13.8f} "
                  f"{r['E_MP2']:12.8f} {r['E_CCSD']:12.8f} {ratio:9.4f}")
        print()

    print("=" * 74)
    print("7. Unitary coupled cluster")
    print("=" * 74)
    print("The unitary ansatz of chapter 10 applied to a real interaction.")
    print("With sigma = T - T^dagger the expansion no longer terminates, so")
    print("the energy is evaluated by exponentiating sigma directly in the")
    print("truncated space, and the result is variational: it lies above the")
    print("exact energy, whereas CCSD carries no such guarantee.")
    print()
    print(f"{'shells':>7s} {'orbitals':>9s} {'E_CCSD':>13s} {'E_UCCSD':>13s} "
          f"{'E_FCI':>13s}")
    for shells in (1, 2):
        dot = QuantumDot(2, shells)
        h, v = dot.matrices()
        h_hf, v_hf, _, _ = hartree_fock(h, v, 2)
        f_hf = h_hf + np.einsum("piqi->pq", v_hf[:, :2, :, :2])
        e_ref = reference_energy(h_hf, v_hf, 2)
        e_ccsd = ccsd(f_hf, v_hf, 2)["energy"]
        ucc = UnitaryCC(h_hf, v_hf, 2, singles=True, doubles=True)
        e_ucc = ucc.optimise()[0]
        exact = fci_energy(h, v, 2)
        print(f"{shells:7d} {dot.n_orbitals:9d} {e_ref+e_ccsd:13.8f} "
              f"{e_ucc:13.8f} {exact:13.8f}")
    print()
    print("For two electrons all three coincide, because singles and doubles")
    print("span the whole space: the unitary and the non-unitary ansatz then")
    print("describe the same state.  The difference between them appears only")
    print("when the truncation actually bites, which for this system means")
    print("six electrons or more.")


if __name__ == "__main__":
    _demo()
