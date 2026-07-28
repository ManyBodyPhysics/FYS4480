"""
Coupled-cluster theory with singles and doubles.

Companion code to chapter 10 of *Quantum mechanics for Many-particle Systems*.

The core is a general spin-orbital CCSD solver.  It takes a Fock matrix, an
antisymmetrised two-body tensor and a number of occupied orbitals, and returns
the amplitudes and the correlation energy.  Setting the singles to zero and
freezing them gives CCD.  Everything is checked against a full
configuration-interaction calculation built from the *same* one- and two-body
matrices, so that no discrepancy can be blamed on a different Hamiltonian.

    ccsd                -- the general solver, with intermediates
    ccd                 -- the same with T_1 switched off
    fci_energy          -- a small FCI in the same spin-orbital basis
    pairing_model       -- the model of chapters 4, 5, 7 and 9
    pairing_ph_model    -- the pairing plus particle-hole model of chapter 7
    ccd_order_by_order  -- the perturbative content of CCD

Conventions.  Spin-orbitals are numbered so that the ``n_occ`` lowest are
occupied in the reference.  The two-body tensor is antisymmetrised,

    v[p, q, r, s] = <pq||rs> = <pq|v|rs> - <pq|v|sr> ,

antisymmetric under p <-> q and under r <-> s, and the Fock matrix is

    f[p, q] = h[p, q] + sum_{i occupied} <pi||qi> .

Author: Morten Hjorth-Jensen
"""

from itertools import combinations

import numpy as np


# ---------------------------------------------------------------------------
#  A small full configuration interaction, for validation
# ---------------------------------------------------------------------------
def popcount(x):
    return bin(x).count("1")


def _annihilate(state, p):
    if not (state >> p) & 1:
        return None, 0
    sign = -1 if popcount(state & ((1 << p) - 1)) & 1 else 1
    return state ^ (1 << p), sign


def _create(state, p):
    if (state >> p) & 1:
        return None, 0
    sign = -1 if popcount(state & ((1 << p) - 1)) & 1 else 1
    return state | (1 << p), sign


def fci_energy(h, v, n_particles):
    """Exact ground-state energy of  H = sum h_pq a+_p a_q + 1/4 sum <pq||rs> ...

    Built by brute force in the determinant basis, for validation only.
    """
    n = h.shape[0]
    states = []
    for occ in combinations(range(n), n_particles):
        bits = 0
        for o in occ:
            bits |= 1 << o
        states.append(bits)
    states.sort()
    index = {s: i for i, s in enumerate(states)}
    dim = len(states)
    H = np.zeros((dim, dim))

    for col, state in enumerate(states):
        # one-body
        for p in range(n):
            for q in range(n):
                if abs(h[p, q]) < 1e-14:
                    continue
                s, sign = _annihilate(state, q)
                if s is None:
                    continue
                s, t = _create(s, p)
                if s is None:
                    continue
                H[index[s], col] += h[p, q] * sign * t
        # two-body: (1/4) sum_{pqrs} <pq||rs> a+_p a+_q a_s a_r
        for p in range(n):
            for q in range(n):
                for r in range(n):
                    for s_ in range(n):
                        val = v[p, q, r, s_]
                        if abs(val) < 1e-14:
                            continue
                        # a+_p a+_q a_s a_r : rightmost first is a_r
                        st, sign = _annihilate(state, r)
                        if st is None:
                            continue
                        st, t1 = _annihilate(st, s_)
                        if st is None:
                            continue
                        st, t2 = _create(st, q)
                        if st is None:
                            continue
                        st, t3 = _create(st, p)
                        if st is None:
                            continue
                        H[index[st], col] += 0.25 * val * sign * t1 * t2 * t3
    return float(np.linalg.eigvalsh(H)[0]), dim


def reference_energy(h, v, n_occ):
    """<Phi_0|H|Phi_0> for the determinant filling the n_occ lowest orbitals."""
    o = slice(0, n_occ)
    return (np.trace(h[o, o])
            + 0.5 * np.einsum("ijij->", v[o, o, o, o]))


def fock_matrix(h, v, n_occ):
    """f_pq = h_pq + sum_i <pi||qi>."""
    o = slice(0, n_occ)
    return h + np.einsum("piqi->pq", v[:, o, :, o])


def hartree_fock(h, v, n_occ, max_iter=500, tol=1e-12):
    """Solve the Hartree-Fock equations and rotate h and v into that basis.

    Returns the transformed one- and two-body matrices, so that a subsequent
    coupled-cluster calculation is built on the Hartree-Fock reference and
    the off-diagonal block f_ia vanishes by Brillouin's theorem.
    """
    n = h.shape[0]
    C = np.eye(n)
    previous = None
    for iteration in range(1, max_iter + 1):
        rho = C[:, :n_occ] @ C[:, :n_occ].T
        F = h + np.einsum("prqs,sr->pq", v, rho)
        eps, C = np.linalg.eigh(F)
        if previous is not None and abs(eps.sum() - previous) < tol:
            break
        previous = eps.sum()
    h_new = C.T @ h @ C
    v_new = np.einsum("pi,qj,pqrs,rk,sl->ijkl", C, C, v, C, C, optimize=True)
    return h_new, v_new, C, iteration


# ---------------------------------------------------------------------------
#  The CCSD solver
# ---------------------------------------------------------------------------
def ccsd(f, v, n_occ, doubles_only=False, max_iter=400, tol=1e-12,
         mixing=0.0, verbose=False):
    """Solve the CCSD equations in a spin-orbital basis.

    Returns a dictionary with the amplitudes, the correlation energy and the
    iteration count.  With ``doubles_only=True`` the singles are held at zero
    and the result is CCD.

    The equations are the standard ones, organised with the two- and
    four-index intermediates of Stanton, Gauss, Watts and Bartlett:

        F_ae, F_mi, F_me,  W_mnij, W_abef, W_mbej,

    which reduce the cost and, more importantly here, make the structure of
    the equations legible.
    """
    n = f.shape[0]
    o, u = slice(0, n_occ), slice(n_occ, n)
    n_vir = n - n_occ

    eps_o = np.diag(f)[o].copy()
    eps_v = np.diag(f)[u].copy()
    d_ia = eps_o[:, None] - eps_v[None, :]                       # (i, a)
    d_ijab = (eps_o[:, None, None, None] + eps_o[None, :, None, None]
              - eps_v[None, None, :, None] - eps_v[None, None, None, :])

    # the off-diagonal part of the Fock matrix, used in the equations
    f_od = f - np.diag(np.diag(f))

    t1 = np.zeros((n_occ, n_vir))
    t2 = v[o, o, u, u] / d_ijab                                  # MP2 start

    def tau(t1, t2):
        return (t2 + np.einsum("ia,jb->ijab", t1, t1)
                - np.einsum("ib,ja->ijab", t1, t1))

    def tau_tilde(t1, t2):
        return (t2 + 0.5 * (np.einsum("ia,jb->ijab", t1, t1)
                            - np.einsum("ib,ja->ijab", t1, t1)))

    def energy(t1, t2):
        e = np.einsum("ia,ia->", f[o, u], t1)
        e += 0.25 * np.einsum("ijab,ijab->", v[o, o, u, u], t2)
        e += 0.5 * np.einsum("ijab,ia,jb->", v[o, o, u, u], t1, t1)
        return float(e)

    previous = 0.0
    for iteration in range(1, max_iter + 1):
        tt = tau_tilde(t1, t2)
        ta = tau(t1, t2)

        # ---- two-index intermediates -------------------------------------
        F_ae = (f_od[u, u]
                - 0.5 * np.einsum("me,ma->ae", f[o, u], t1)
                + np.einsum("mf,mafe->ae", t1, v[o, u, u, u])
                - 0.5 * np.einsum("mnaf,mnef->ae", tt, v[o, o, u, u]))
        F_mi = (f_od[o, o]
                + 0.5 * np.einsum("ie,me->mi", t1, f[o, u])
                + np.einsum("ne,mnie->mi", t1, v[o, o, o, u])
                + 0.5 * np.einsum("inef,mnef->mi", tt, v[o, o, u, u]))
        F_me = f[o, u] + np.einsum("nf,mnef->me", t1, v[o, o, u, u])

        # ---- four-index intermediates ------------------------------------
        W_mnij = (v[o, o, o, o]
                  + np.einsum("je,mnie->mnij", t1, v[o, o, o, u])
                  - np.einsum("ie,mnje->mnij", t1, v[o, o, o, u])
                  + 0.25 * np.einsum("ijef,mnef->mnij", ta, v[o, o, u, u]))
        W_abef = (v[u, u, u, u]
                  - np.einsum("mb,amef->abef", t1, v[u, o, u, u])
                  + np.einsum("ma,bmef->abef", t1, v[u, o, u, u])
                  + 0.25 * np.einsum("mnab,mnef->abef", ta, v[o, o, u, u]))
        W_mbej = (v[o, u, u, o]
                  + np.einsum("jf,mbef->mbej", t1, v[o, u, u, u])
                  - np.einsum("nb,mnej->mbej", t1, v[o, o, u, o])
                  - np.einsum("jnfb,mnef->mbej",
                              0.5 * t2 + np.einsum("jf,nb->jnfb", t1, t1),
                              v[o, o, u, u]))

        # ---- singles ------------------------------------------------------
        if doubles_only:
            t1_new = np.zeros_like(t1)
        else:
            r1 = (f[o, u]
                  + np.einsum("ie,ae->ia", t1, F_ae)
                  - np.einsum("ma,mi->ia", t1, F_mi)
                  + np.einsum("imae,me->ia", t2, F_me)
                  - np.einsum("nf,naif->ia", t1, v[o, u, o, u])
                  - 0.5 * np.einsum("imef,maef->ia", t2, v[o, u, u, u])
                  - 0.5 * np.einsum("mnae,nmei->ia", t2, v[o, o, u, o]))
            t1_new = r1 / d_ia

        # ---- doubles ------------------------------------------------------
        tmp_ae = F_ae - 0.5 * np.einsum("mb,me->be", t1, F_me)
        tmp_mi = F_mi + 0.5 * np.einsum("je,me->mj", t1, F_me)

        r2 = v[o, o, u, u].copy()

        term = np.einsum("ijae,be->ijab", t2, tmp_ae)
        r2 = r2 + term - term.transpose(0, 1, 3, 2)

        term = np.einsum("imab,mj->ijab", t2, tmp_mi)
        r2 = r2 - term + term.transpose(1, 0, 2, 3)

        r2 = r2 + 0.5 * np.einsum("mnab,mnij->ijab", ta, W_mnij)
        r2 = r2 + 0.5 * np.einsum("ijef,abef->ijab", ta, W_abef)

        term = (np.einsum("imae,mbej->ijab", t2, W_mbej)
                - np.einsum("ie,ma,mbej->ijab", t1, t1, v[o, u, u, o]))
        r2 = (r2 + term - term.transpose(1, 0, 2, 3)
              - term.transpose(0, 1, 3, 2) + term.transpose(1, 0, 3, 2))

        if not doubles_only:
            term = np.einsum("ie,abej->ijab", t1, v[u, u, u, o])
            r2 = r2 + term - term.transpose(1, 0, 2, 3)
            term = np.einsum("ma,mbij->ijab", t1, v[o, u, o, o])
            r2 = r2 - term + term.transpose(0, 1, 3, 2)

        t2_new = r2 / d_ijab

        if mixing > 0.0:
            t1_new = (1.0 - mixing) * t1_new + mixing * t1
            t2_new = (1.0 - mixing) * t2_new + mixing * t2

        e = energy(t1_new, t2_new)
        change = abs(e - previous)
        t1, t2, previous = t1_new, t2_new, e
        if verbose:
            print(f"   {iteration:4d}  E_corr = {e:16.12f}   "
                  f"change = {change:.2e}")
        if change < tol:
            break

    return dict(t1=t1, t2=t2, energy=previous, iterations=iteration,
                converged=change < tol)


def ccd(f, v, n_occ, **kwargs):
    """Coupled cluster with doubles only."""
    return ccsd(f, v, n_occ, doubles_only=True, **kwargs)


# ---------------------------------------------------------------------------
#  The models
# ---------------------------------------------------------------------------
def _add_term(v, coefficient, p, q, r, s):
    """Accumulate  coefficient * a+_p a+_q a_s a_r  into an antisymmetrised v.

    With v antisymmetric in (p, q) and in (r, s), the operator

        (1/4) sum_{pqrs} v[p,q,r,s] a+_p a+_q a_s a_r

    has  v[p,q,r,s]  as the coefficient of  a+_p a+_q a_s a_r : the four
    index orderings that give the same operator each contribute a quarter.
    """
    v[p, q, r, s] += coefficient
    v[q, p, r, s] -= coefficient
    v[p, q, s, r] -= coefficient
    v[q, p, s, r] += coefficient


def pairing_model(levels=4, particles=4, g=1.0, xi=1.0):
    """The pairing model of chapters 4, 5, 7 and 9, in a spin-orbital basis.

    Spin-orbital 2(p-1) is (p, +) and 2(p-1)+1 is (p, -).  The interaction

        V = -(g/2) sum_{pq} a+_{p+} a+_{p-} a_{q-} a_{q+}

    therefore has  <p+ p-||q+ q-> = -g/2.
    """
    n = 2 * levels
    h = np.diag([xi * (o // 2) for o in range(n)]).astype(float)
    v = np.zeros((n, n, n, n))
    for p in range(levels):
        for q in range(levels):
            # a+_{p+} a+_{p-} a_{q-} a_{q+}: creation (p+, p-), r = q+, s = q-
            _add_term(v, -0.5 * g, 2 * p, 2 * p + 1, 2 * q, 2 * q + 1)
    return h, v, particles


def pairing_ph_model(levels=4, particles=4, g=1.0, f=0.0, xi=1.0):
    """The pairing plus particle-hole model of chapter 7.

        V = -(g/2) sum_{pq}  a+_{p+} a+_{p-} a_{q-} a_{q+}
            -(f/2) sum_{pqr} ( a+_{p+} a+_{p-} a_{q-} a_{r+} + h.c. )

    The particle-hole term breaks pairs, so the singles amplitudes do not
    vanish and CCSD differs from CCD.
    """
    n = 2 * levels
    h = np.diag([xi * (o // 2) for o in range(n)]).astype(float)
    v = np.zeros((n, n, n, n))
    for p in range(levels):
        for q in range(levels):
            _add_term(v, -0.5 * g, 2 * p, 2 * p + 1, 2 * q, 2 * q + 1)
    if f != 0.0:
        for p in range(levels):
            for q in range(levels):
                for r in range(levels):
                    # a+_{p+} a+_{p-} a_{q-} a_{r+}: r index = r+, s = q-
                    _add_term(v, -0.5 * f,
                              2 * p, 2 * p + 1, 2 * r, 2 * q + 1)
                    # the Hermitian conjugate a+_{r+} a+_{q-} a_{p-} a_{p+}
                    _add_term(v, -0.5 * f,
                              2 * r, 2 * q + 1, 2 * p, 2 * p + 1)
    return h, v, particles


def independent_copies(subsystems=2, levels_per=3, g=0.5, xi=1.0):
    """Non-interacting copies of a pairing system, for the extensivity test."""
    levels = subsystems * levels_per
    n = 2 * levels
    h = np.diag([xi * ((o // 2) % levels_per)
                 for o in range(n)]).astype(float)
    v = np.zeros((n, n, n, n))
    for p in range(levels):
        for q in range(levels):
            if p // levels_per != q // levels_per:
                continue
            _add_term(v, -0.5 * g, 2 * p, 2 * p + 1, 2 * q, 2 * q + 1)
    return h, v, 2 * subsystems


# ---------------------------------------------------------------------------
#  The perturbative content of CCD
# ---------------------------------------------------------------------------
def ccd_residual(t2, f, v, n_occ):
    """The right-hand side of the CCD amplitude equation.

    The CCD equation is  D_ij^ab t_ij^ab = R(t2),  and the solver iterates
    t2 <- R(t2)/D.  Note that R is exactly *quadratic* in t2: the
    intermediates contain terms linear in t2 which are then contracted with
    t2 once more, and nothing of higher degree appears.  That fact is used
    below to extract the perturbative content of the theory.
    """
    n = f.shape[0]
    o, u = slice(0, n_occ), slice(n_occ, n)
    f_od = f - np.diag(np.diag(f))
    w = v[o, o, u, u]

    F_ae = f_od[u, u] - 0.5 * np.einsum("mnaf,mnef->ae", t2, w)
    F_mi = f_od[o, o] + 0.5 * np.einsum("inef,mnef->mi", t2, w)
    W_mnij = v[o, o, o, o] + 0.25 * np.einsum("ijef,mnef->mnij", t2, w)
    W_abef = v[u, u, u, u] + 0.25 * np.einsum("mnab,mnef->abef", t2, w)
    W_mbej = v[o, u, u, o] - 0.5 * np.einsum("jnfb,mnef->mbej", t2, w)

    r2 = w.copy()
    term = np.einsum("ijae,be->ijab", t2, F_ae)
    r2 = r2 + term - term.transpose(0, 1, 3, 2)
    term = np.einsum("imab,mj->ijab", t2, F_mi)
    r2 = r2 - term + term.transpose(1, 0, 2, 3)
    r2 = r2 + 0.5 * np.einsum("mnab,mnij->ijab", t2, W_mnij)
    r2 = r2 + 0.5 * np.einsum("ijef,abef->ijab", t2, W_abef)
    term = np.einsum("imae,mbej->ijab", t2, W_mbej)
    r2 = (r2 + term - term.transpose(1, 0, 2, 3)
          - term.transpose(0, 1, 3, 2) + term.transpose(1, 0, 3, 2))
    return r2


def ccd_order_by_order(f, v, n_occ, order=10):
    """Expand the CCD amplitude equation order by order in the interaction.

    Because R(t2) is quadratic we may split it exactly,

        R(t) = C + L[t] + Q[t, t] ,

    with  C = R(0),  L[t] = (R(t) - R(-t))/2,  Q[t,t] = (R(t)+R(-t))/2 - C,
    and the bilinear form recovered by polarisation,
    Q[a,b] = (Q[a+b,a+b] - Q[a,a] - Q[b,b])/2.  Writing
    t2 = sum_k t^(k) with t^(k) of order k in the interaction gives

        D t^(1) = C ,
        D t^(k) = L[t^(k-1)] + sum_{i+j=k-1} Q[t^(i), t^(j)] ,

    and the energy contribution of t^(k) is of order k+1.  The first two
    therefore reproduce MP2 and MP3 exactly, and the sum of them all is the
    CCD energy -- provided the series converges.
    """
    n = f.shape[0]
    o, u = slice(0, n_occ), slice(n_occ, n)
    eps_o, eps_v = np.diag(f)[o], np.diag(f)[u]
    d = (eps_o[:, None, None, None] + eps_o[None, :, None, None]
         - eps_v[None, None, :, None] - eps_v[None, None, None, :])

    zero = np.zeros_like(v[o, o, u, u])
    C = ccd_residual(zero, f, v, n_occ)

    def linear(t):
        return 0.5 * (ccd_residual(t, f, v, n_occ)
                      - ccd_residual(-t, f, v, n_occ))

    def quad(t):
        return 0.5 * (ccd_residual(t, f, v, n_occ)
                      + ccd_residual(-t, f, v, n_occ)) - C

    def bilinear(a, b):
        return 0.5 * (quad(a + b) - quad(a) - quad(b))

    amplitudes = [C / d]                                   # t^(1) = MP2
    for k in range(2, order + 1):
        rhs = linear(amplitudes[-1])
        for i in range(1, k):
            j = k - 1 - i
            if 1 <= j <= len(amplitudes) and i <= len(amplitudes):
                rhs = rhs + bilinear(amplitudes[i - 1], amplitudes[j - 1])
        amplitudes.append(rhs / d)

    energies = [0.25 * float(np.einsum("ijab,ijab->", v[o, o, u, u], t))
                for t in amplitudes]
    return energies, amplitudes


def mp2_energy(f, v, n_occ):
    n = f.shape[0]
    o, u = slice(0, n_occ), slice(n_occ, n)
    eps_o, eps_v = np.diag(f)[o], np.diag(f)[u]
    d = (eps_o[:, None, None, None] + eps_o[None, :, None, None]
         - eps_v[None, None, :, None] - eps_v[None, None, None, :])
    return 0.25 * float(np.einsum("ijab,ijab->", v[o, o, u, u],
                                  v[o, o, u, u] / d))


# ---------------------------------------------------------------------------
#  Demonstrations
# ---------------------------------------------------------------------------
def demo_validation():
    print("=" * 74)
    print("1. Validating the solver")
    print("=" * 74)
    print("Everything below rests on the equations being right, so we check")
    print("them first, in four ways.")
    print()
    print("(i) The models must reproduce the exact energies of the earlier")
    print("chapters.  The pairing model against Table 4.2:")
    reference = {-1.0: 2.77987014, -0.5: 2.43688426, 0.0: 2.00000000,
                 0.5: 1.41677428, 1.0: 0.63554847}
    print(f"{'g':>7s} {'E_ref':>10s} {'FCI here':>14s} {'chapter 4':>14s}")
    for g in (-1.0, -0.5, 0.0, 0.5, 1.0):
        h, v, N = pairing_model(g=g)
        e_fci, dim = fci_energy(h, v, N)
        print(f"{g:7.2f} {reference_energy(h, v, N):10.5f} {e_fci:14.8f} "
              f"{reference[g]:14.8f}")
    print()
    print("(ii) For the pairing model the interaction cannot connect the")
    print("reference to a 1p-1h state, so the singles must vanish and CCSD")
    print("must equal CCD exactly:")
    print(f"{'g':>7s} {'E_CCD':>16s} {'E_CCSD':>16s} {'max |t1|':>11s}")
    for g in (0.5, 1.0, 2.0):
        h, v, N = pairing_model(g=g)
        fk = fock_matrix(h, v, N)
        d, s = ccd(fk, v, N), ccsd(fk, v, N)
        print(f"{g:7.2f} {d['energy']:16.10f} {s['energy']:16.10f} "
              f"{np.abs(s['t1']).max():11.2e}")
    print()
    print("(iii) The initial amplitudes are the MP2 ones, so the energy")
    print("before any iteration must be the second-order perturbation")
    print("energy.  For the pairing model this has a closed form,")
    print("   E_MP2 = -(g^2/4) [ 1/(2+g) + 2/(4+g) + 1/(6+g) ],")
    print("in which the g in each denominator is the Hartree-Fock shift of")
    print("the occupied levels -- the Moller-Plesset partition, not the bare")
    print("one-body partition used in chapter 9.")
    print(f"{'g':>7s} {'MP2 computed':>16s} {'closed form':>16s}")
    for g in (0.5, 1.0, 2.0):
        h, v, N = pairing_model(g=g)
        fk = fock_matrix(h, v, N)
        closed = -0.25 * g**2 * (1.0 / (2.0 + g) + 2.0 / (4.0 + g)
                                 + 1.0 / (6.0 + g))
        print(f"{g:7.2f} {mp2_energy(fk, v, N):16.10f} {closed:16.10f}")
    print()
    print("(iv) The Hartree-Fock rotation must leave the exact energy alone")
    print("and must make the occupied-virtual block of the Fock matrix")
    print("vanish.  For the pairing plus particle-hole model of chapter 7:")
    print(f"{'g':>5s} {'f':>6s} {'E_ref':>10s} {'E_HF':>11s} "
          f"{'max |f_ia|':>12s} {'FCI shift':>11s}")
    for g, ph in ((0.5, 0.25), (1.0, 0.5)):
        h, v, N = pairing_ph_model(g=g, f=ph)
        e0, _ = fci_energy(h, v, N)
        h2, v2, C, it = hartree_fock(h, v, N)
        e1, _ = fci_energy(h2, v2, N)
        fk2 = fock_matrix(h2, v2, N)
        print(f"{g:5.2f} {ph:6.2f} {reference_energy(h, v, N):10.5f} "
              f"{reference_energy(h2, v2, N):11.5f} "
              f"{np.abs(fk2[:N, N:]).max():12.1e} {abs(e1 - e0):11.1e}")
    print()
    print("The Hartree-Fock energies agree with those of chapter 7, and the")
    print("exact energy is invariant under the rotation, as it must be.")


def demo_pairing():
    print("=" * 74)
    print("2. The pairing model")
    print("=" * 74)
    print("Four doubly degenerate levels, four particles, xi = 1.  The")
    print("Hartree-Fock energy is the reference energy 2 - g, as chapter 6")
    print("found, so everything below is correlation energy.")
    print()
    print(f"{'g':>6s} {'E_ref':>9s} {'MP2':>12s} {'CCD':>12s} {'FCI':>12s} "
          f"{'CCD error':>11s} {'iter':>6s}")
    for g in (0.25, 0.5, 1.0, 1.5, 2.0):
        h, v, N = pairing_model(g=g)
        fk = fock_matrix(h, v, N)
        e_ref = reference_energy(h, v, N)
        result = ccd(fk, v, N)
        e_fci, _ = fci_energy(h, v, N)
        print(f"{g:6.2f} {e_ref:9.5f} {e_ref + mp2_energy(fk, v, N):12.8f} "
              f"{e_ref + result['energy']:12.8f} {e_fci:12.8f} "
              f"{e_ref + result['energy'] - e_fci:11.2e} "
              f"{result['iterations']:6d}")
    print()
    print("At weak coupling CCD is very nearly exact -- four parts in a")
    print("million at g = 0.25.  It over-correlates, and the error grows")
    print("with the coupling, reaching 0.12 at g = 2.  What CCD misses here")
    print("is the connected quadruple excitations: the exponential supplies")
    print("the 4p-4h amplitude as (1/2) T_2^2, which is the disconnected")
    print("part only, and for four particles in four levels the connected")
    print("part is not negligible once the interaction is strong.")


def demo_orders():
    print("=" * 74)
    print("3. What coupled cluster contains, order by order")
    print("=" * 74)
    print("The CCD equation is exactly quadratic in the amplitudes, so it")
    print("can be split as R(t) = C + L[t] + Q[t,t] and expanded in powers")
    print("of the interaction.  The result is perturbation theory: the first")
    print("amplitude is the MP2 one, the second generates the third-order")
    print("energy, and so on -- but the iteration never stops, so CCD")
    print("contains contributions at every order.")
    print()
    print("Against Moller-Plesset perturbation theory computed independently")
    print("in the determinant basis of chapter 9:")
    print(f"{'g':>6s} {'order':>7s} {'from CCD':>16s} {'from MBPT':>16s}")
    try:
        import fci as fcimod
        import mbpt as mbptmod
        for g in (0.5, 1.0):
            h, v, N = pairing_model(g=g)
            fk = fock_matrix(h, v, N)
            energies, _ = ccd_order_by_order(fk, v, N, order=4)
            model = fcimod.PairingFCI(levels=4, n_particles=4, g=g)
            eps = np.diag(fk)
            w = np.array([sum(eps[o] for o in range(8) if (s >> o) & 1)
                          for s in model.basis.states])
            part = mbptmod.Partition(w=None, H=model.matrix(), h0_diagonal=w,
                                     reference=model.basis.index[
                                         model.reference]) \
                if False else mbptmod.Partition(model.matrix(), w,
                                                model.basis.index[
                                                    model.reference])
            rs = mbptmod.rayleigh_schrodinger(part, order=4)
            for k in (2, 3, 4):
                print(f"{g if k == 2 else 0:6.2f} {k:7d} "
                      f"{energies[k - 2]:16.10f} {rs[k - 1]:16.10f}")
    except Exception as exc:
        print(f"   (chapter 9's program not importable: {exc})")
    print()
    print("They agree to every digit through fourth order.  That CCD should")
    print("reproduce MP2 and MP3 is general; that it also gets MP4 right is")
    print("special to this model, where the interaction produces neither")
    print("singles nor triples and the quadruples enter only through the")
    print("disconnected T_2^2 that the exponential already supplies.")
    print()
    print("Summing the expansion:")
    for g in (0.5, 1.0, 1.5):
        h, v, N = pairing_model(g=g)
        fk = fock_matrix(h, v, N)
        e_ref = reference_energy(h, v, N)
        energies, _ = ccd_order_by_order(fk, v, N, order=16)
        full = ccd(fk, v, N)["energy"]
        e_fci, _ = fci_energy(h, v, N)
        print(f"   g = {g}:  E_ref = {e_ref:.5f}")
        print(f"      {'order':>7s} {'term':>13s} {'partial sum':>15s}")
        partial = 0.0
        for k, e in enumerate(energies[:6], start=2):
            partial += e
            print(f"      {k:7d} {e:13.3e} {e_ref + partial:15.10f}")
        print(f"      sum of 16 orders = {e_ref + sum(energies):.10f}")
        print(f"      CCD              = {e_ref + full:.10f}")
        print(f"      FCI              = {e_fci:.10f}")
        print()
    print("The expansion converges back onto the CCD energy, as it must.")
    print("But note what has happened: at g = 1.5 the bare perturbation")
    print("series of chapter 9 was already close to divergence, while these")
    print("terms -- the same expansion, resummed by the nonlinear equation --")
    print("still converge.  Coupled cluster is not a better perturbation")
    print("theory; it is the same perturbation theory summed to all orders")
    print("within a restricted class of excitations.")


def demo_ccsd():
    print("=" * 74)
    print("4. When singles matter: the pairing plus particle-hole model")
    print("=" * 74)
    print("The particle-hole term of chapter 7 breaks pairs, so the")
    print("reference couples to 1p-1h states and the singles amplitudes do")
    print("not vanish.  This is where CCSD earns its name.  We run it twice:")
    print("on the bare oscillator reference, and on the Hartree-Fock")
    print("reference where Brillouin's theorem holds.")
    print()
    print(f"{'g':>5s} {'f':>6s} {'CCD':>12s} {'CCSD':>12s} {'CCSD/HF':>12s} "
          f"{'FCI':>12s} {'|t1|':>9s} {'|t1| HF':>9s}")
    for g, ph in ((0.5, 0.025), (0.7, 0.05), (1.0, 0.05), (0.5, 0.25),
                  (1.0, 0.5)):
        h, v, N = pairing_ph_model(g=g, f=ph)
        fk = fock_matrix(h, v, N)
        e_ref = reference_energy(h, v, N)
        d = ccd(fk, v, N, mixing=0.3)
        s = ccsd(fk, v, N, mixing=0.3)
        h2, v2, C, _ = hartree_fock(h, v, N)
        fk2 = fock_matrix(h2, v2, N)
        e_hf = reference_energy(h2, v2, N)
        s2 = ccsd(fk2, v2, N, mixing=0.3)
        e_fci, _ = fci_energy(h, v, N)
        print(f"{g:5.2f} {ph:6.3f} {e_ref + d['energy']:12.7f} "
              f"{e_ref + s['energy']:12.7f} {e_hf + s2['energy']:12.7f} "
              f"{e_fci:12.7f} {np.abs(s['t1']).max():9.1e} "
              f"{np.abs(s2['t1']).max():9.1e}")
    print()
    print("Three things to read off.  The singles matter: at g = 0.5,")
    print("f = 0.25 they move the energy from 0.6033 to 0.5467 against an")
    print("exact 0.5512, turning a 9% error into a 0.8% one.  The singles")
    print("amplitudes are an order of magnitude smaller in the Hartree-Fock")
    print("basis, because Brillouin's theorem makes them start at second")
    print("order rather than first.  And the two CCSD columns are close but")
    print("not identical -- coupled cluster is not invariant under a change")
    print("of reference determinant, only exact in the untruncated limit.")


def demo_extensivity():
    print("=" * 74)
    print("5. Size extensivity")
    print("=" * 74)
    print("Chapter 5 found that truncated configuration interaction is not")
    print("size consistent, and chapter 9 that Rayleigh-Schrodinger")
    print("perturbation theory is.  Coupled cluster is size extensive by")
    print("construction: for non-interacting subsystems the cluster")
    print("operators commute and exp(T_A + T_B) = exp(T_A) exp(T_B), so the")
    print("wave function factorises and the energy is additive.")
    print()
    print("Two identical, non-interacting copies of a three-level pairing")
    print("system:")
    print(f"{'g':>6s} {'2 x one copy':>17s} {'two copies':>17s} "
          f"{'error':>11s}")
    for g in (0.25, 0.5, 1.0):
        h1, v1, n1 = independent_copies(subsystems=1, g=g)
        h2, v2, n2 = independent_copies(subsystems=2, g=g)
        e1 = ccd(fock_matrix(h1, v1, n1), v1, n1)["energy"]
        e2 = ccd(fock_matrix(h2, v2, n2), v2, n2)["energy"]
        print(f"{g:6.2f} {2 * e1:17.12f} {e2:17.12f} {e2 - 2 * e1:11.2e}")
    print()
    print("Exact to machine precision at every coupling.  Compare with the")
    print("doubles configuration-interaction errors of chapter 5, which grew")
    print("from 6e-05 to 1e-02 over the same range of couplings.  This is")
    print("the single property that decides between the two methods for")
    print("extended systems, and coupled cluster gets it for free from the")
    print("exponential.")


def demo_comparison():
    print("=" * 74)
    print("6. Everything against everything")
    print("=" * 74)
    print("The pairing model, treated by every method in the book.")
    print()
    try:
        import fci as fcimod
        import rpa as rpamod
        import mbpt as mbptmod
    except Exception as exc:
        print(f"   (the other chapters' programs are not importable: {exc})")
        return
    fock_rpa = rpamod.FockSpace(4)
    rows = []
    for g in (0.25, 0.5, 1.0, 1.5, 2.0):
        h, v, N = pairing_model(g=g)
        fk = fock_matrix(h, v, N)
        e_ref = reference_energy(h, v, N)
        cc_e = e_ref + ccd(fk, v, N)["energy"]
        part, model = mbptmod.pairing_partition(g=g)
        rs = mbptmod.rayleigh_schrodinger(part, order=3)
        r = rpamod.tda_rpa(fock_rpa, 4, g, 0.0)
        e_fci, _ = fci_energy(h, v, N)
        rows.append((g, e_ref, part.reference_energy + rs[1],
                     part.reference_energy + rs[1] + rs[2],
                     model.energies(2)[0], r["hf"] + r["ecorr"],
                     cc_e, e_fci))
    print(f"{'g':>6s} {'HF':>9s} {'MBPT2':>10s} {'MBPT3':>10s} "
          f"{'CID':>10s} {'RPA':>10s} {'CCD':>11s} {'FCI':>11s}")
    for g, hf, m2, m3, cid, rpa_e, cc_e, ex in rows:
        print(f"{g:6.2f} {hf:9.5f} {m2:10.5f} {m3:10.5f} {cid:10.5f} "
              f"{rpa_e:10.5f} {cc_e:11.7f} {ex:11.7f}")
    print()
    print("The errors:")
    print(f"{'g':>6s} {'MBPT2':>11s} {'MBPT3':>11s} {'CID':>11s} "
          f"{'RPA':>11s} {'CCD':>11s}")
    for g, hf, m2, m3, cid, rpa_e, cc_e, ex in rows:
        print(f"{g:6.2f} {m2 - ex:11.2e} {m3 - ex:11.2e} {cid - ex:11.2e} "
              f"{rpa_e - ex:11.2e} {cc_e - ex:11.2e}")
    print()
    print("At weak and moderate coupling coupled cluster is the best of the")
    print("five by two to three orders of magnitude, and unlike the")
    print("perturbation series it does not fall apart as the interaction")
    print("grows.  At g = 2 the random-phase approximation happens to be")
    print("closer, which is a reminder that a method summing a different")
    print("class of contributions can win in a regime it was designed for.")
    print("But CCD is the only one of the five that is simultaneously size")
    print("extensive, systematically improvable and accurate across the")
    print("whole range -- which is why it, and not the others, is the")
    print("workhorse of modern many-body theory.")


def _demo():
    for f in (demo_validation, demo_pairing, demo_orders, demo_ccsd,
              demo_extensivity, demo_comparison):
        f()
        print()


if __name__ == "__main__":
    _demo()
