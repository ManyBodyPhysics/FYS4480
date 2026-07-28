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
                        st, sign = _annihilate(state, s_)
                        if st is None:
                            continue
                        st, t1 = _annihilate(st, r)
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
def pairing_model(levels=4, particles=4, g=1.0, xi=1.0):
    """The pairing model of chapters 4, 5, 7 and 9, in a spin-orbital basis.

    Spin-orbital 2(p-1) is (p, +) and 2(p-1)+1 is (p, -).  The interaction is

        V = -(g/2) sum_{pq} a+_{p+} a+_{p-} a_{q-} a_{q+} ,

    which in antisymmetrised form has <p+ p-||q+ q-> = -g.
    """
    n = 2 * levels
    h = np.diag([xi * (o // 2) for o in range(n)])
    v = np.zeros((n, n, n, n))
    for p in range(levels):
        for q in range(levels):
            a, b = 2 * p, 2 * p + 1
            c, d = 2 * q, 2 * q + 1
            # <ab||cd> = -g, with the four antisymmetric orderings
            v[a, b, c, d] += -g
            v[b, a, c, d] += +g
            v[a, b, d, c] += +g
            v[b, a, d, c] += -g
    return h, v, particles


def pairing_ph_model(levels=4, particles=4, g=1.0, f=0.0, xi=1.0):
    """The pairing plus particle-hole model of chapter 7.

        V = -(g/2) sum_{pq}  a+_{p+} a+_{p-} a_{q-} a_{q+}
            -(f/2) sum_{pqr} ( a+_{p+} a+_{p-} a_{q-} a_{r+} + h.c. )

    The particle-hole term breaks pairs, so the singles amplitudes do not
    vanish and CCSD differs from CCD.
    """
    n = 2 * levels
    h = np.diag([xi * (o // 2) for o in range(n)])
    w = np.zeros((n, n, n, n))          # non-antisymmetrised <pq|v|rs>
    for p in range(levels):
        a, b = 2 * p, 2 * p + 1         # (p,+), (p,-)
        for q in range(levels):
            c, d = 2 * q, 2 * q + 1     # (q,+), (q,-)
            w[a, b, c, d] += -0.5 * g
        for q in range(levels):
            for r in range(levels):
                c = 2 * q + 1           # (q,-)
                d = 2 * r               # (r,+)
                # -(f/2) a+_{p+} a+_{p-} a_{q-} a_{r+} : <p+ p-|v|r+ q->
                w[a, b, d, c] += -0.5 * f
                # the Hermitian conjugate
                w[d, c, a, b] += -0.5 * f
    # antisymmetrise
    v = w - w.transpose(0, 1, 3, 2)
    v = 0.5 * (v - v.transpose(1, 0, 2, 3))
    v = 0.5 * (v + v.transpose(2, 3, 0, 1))
    return h, v, particles


def independent_copies(subsystems=2, levels_per=3, g=0.5, xi=1.0):
    """Non-interacting copies of a pairing system, for the extensivity test."""
    levels = subsystems * levels_per
    n = 2 * levels
    h = np.diag([xi * ((o // 2) % levels_per) for o in range(n)])
    v = np.zeros((n, n, n, n))
    for p in range(levels):
        for q in range(levels):
            if p // levels_per != q // levels_per:
                continue
            a, b = 2 * p, 2 * p + 1
            c, d = 2 * q, 2 * q + 1
            v[a, b, c, d] += -g
            v[b, a, c, d] += +g
            v[a, b, d, c] += +g
            v[b, a, d, c] += -g
    return h, v, 2 * subsystems


# ---------------------------------------------------------------------------
#  The perturbative content of CCD
# ---------------------------------------------------------------------------
def ccd_order_by_order(f, v, n_occ, order=8):
    """Expand the CCD amplitude equation order by order in the interaction.

    Writing t2 = t2^(1) + t2^(2) + ... and collecting powers, the linear
    terms of the CCD equation generate one new order at a time and the
    quadratic term couples lower orders.  The energy through order n is
    then the sum of the corresponding energy contributions, and the first
    two reproduce MP2 and MP3 exactly.
    """
    n = f.shape[0]
    o, u = slice(0, n_occ), slice(n_occ, n)
    eps_o, eps_v = np.diag(f)[o], np.diag(f)[u]
    d = (eps_o[:, None, None, None] + eps_o[None, :, None, None]
         - eps_v[None, None, :, None] - eps_v[None, None, None, :])
    f_od = f - np.diag(np.diag(f))

    def linear(t2):
        """The terms of the CCD equation that are linear in t2."""
        term = np.einsum("ijae,be->ijab", t2, f_od[u, u])
        out = term - term.transpose(0, 1, 3, 2)
        term = np.einsum("imab,mj->ijab", t2, f_od[o, o])
        out = out - term + term.transpose(1, 0, 2, 3)
        out = out + 0.5 * np.einsum("mnab,mnij->ijab", t2, v[o, o, o, o])
        out = out + 0.5 * np.einsum("ijef,abef->ijab", t2, v[u, u, u, u])
        term = np.einsum("imae,mbej->ijab", t2, v[o, u, u, o])
        out = (out + term - term.transpose(1, 0, 2, 3)
               - term.transpose(0, 1, 3, 2) + term.transpose(1, 0, 3, 2))
        return out

    def quadratic(ta, tb):
        """The terms quadratic in t2, symmetrised between the two arguments."""
        w = v[o, o, u, u]
        out = 0.25 * np.einsum("ijef,mnef,mnab->ijab", ta, w, tb)
        term = 0.5 * np.einsum("imae,mnef,njbf->ijab", ta, w, tb)
        out = (out + term - term.transpose(1, 0, 2, 3)
               - term.transpose(0, 1, 3, 2) + term.transpose(1, 0, 3, 2))
        term = -0.5 * np.einsum("imef,mnef,njab->ijab", ta, w, tb)
        out = out + term - term.transpose(1, 0, 2, 3)
        term = -0.5 * np.einsum("mnae,mnef,ijbf->ijab", ta, w, tb)
        out = out + term - term.transpose(0, 1, 3, 2)
        return out

    amplitudes = [v[o, o, u, u] / d]                # t2^(1), the MP2 amplitude
    for k in range(2, order + 1):
        rhs = linear(amplitudes[-1])
        for i in range(len(amplitudes)):
            for j in range(len(amplitudes)):
                if i + j + 2 == k:                  # orders (i+1) and (j+1)
                    rhs = rhs + 0.5 * (quadratic(amplitudes[i], amplitudes[j]))
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
    print("Everything below rests on the CCSD equations being right, so we")
    print("check them first.  Three tests.")
    print()
    print("(i) For the pairing model the interaction cannot connect the")
    print("reference to a 1p-1h state, so the singles must vanish and CCSD")
    print("must equal CCD:")
    print(f"{'g':>6s} {'E_CCD':>15s} {'E_CCSD':>15s} {'max |t1|':>11s}")
    for g in (0.5, 1.0, 2.0):
        h, v, N = pairing_model(g=g)
        fk = fock_matrix(h, v, N)
        d = ccd(fk, v, N)
        s = ccsd(fk, v, N)
        print(f"{g:6.2f} {d['energy']:15.10f} {s['energy']:15.10f} "
              f"{np.abs(s['t1']).max():11.2e}")
    print()
    print("(ii) The first iteration, started from zero amplitudes, must")
    print("reproduce the second-order perturbation energy of chapter 9 in")
    print("the Hartree-Fock partition:")
    print(f"{'g':>6s} {'first iteration':>17s} {'MP2':>15s} "
          f"{'closed form':>15s}")
    for g in (0.5, 1.0, 2.0):
        h, v, N = pairing_model(g=g)
        fk = fock_matrix(h, v, N)
        first = ccd(fk, v, N, max_iter=1)
        closed = -0.25 * g**2 * (1.0 / (2.0 + g) + 2.0 / (4.0 + g)
                                 + 1.0 / (6.0 + g))
        print(f"{g:6.2f} {first['energy']:17.10f} {mp2_energy(fk, v, N):15.10f} "
              f"{closed:15.10f}")
    print()
    print("(iii) For a system small enough to diagonalise, CCSD must lie")
    print("between the Hartree-Fock and the exact energies, and CCSD with a")
    print("complete cluster operator would be exact.  We compare against a")
    print("full configuration-interaction calculation built from the same")
    print("one- and two-body matrices:")
    print(f"{'g':>6s} {'f':>6s} {'E_ref':>12s} {'E_CCSD':>14s} "
          f"{'E_FCI':>14s} {'error':>11s}")
    for g, ph in ((0.5, 0.0), (1.0, 0.0), (0.5, 0.25), (1.0, 0.5)):
        h, v, N = pairing_ph_model(g=g, f=ph)
        fk = fock_matrix(h, v, N)
        e_ref = reference_energy(h, v, N)
        s = ccsd(fk, v, N, mixing=0.3)
        e_fci, dim = fci_energy(h, v, N)
        print(f"{g:6.2f} {ph:6.2f} {e_ref:12.6f} "
              f"{e_ref + s['energy']:14.8f} {e_fci:14.8f} "
              f"{e_ref + s['energy'] - e_fci:11.2e}")
    print()
    print("The CCSD energies sit within a few parts in ten thousand of the")
    print("exact ones, which for a four-particle system in an eight-orbital")
    print("space is about as much as singles and doubles can do.")


def demo_pairing():
    print("=" * 74)
    print("2. The pairing model")
    print("=" * 74)
    print("Four doubly degenerate levels, four particles, xi = 1.  The")
    print("Hartree-Fock energy is the reference energy 2 - g, as chapter 6")
    print("found, so everything below is correlation energy.")
    print()
    print(f"{'g':>6s} {'E_ref':>9s} {'MP2':>12s} {'CCD':>12s} "
          f"{'FCI':>12s} {'CCD error':>12s} {'iterations':>11s}")
    for g in (0.25, 0.5, 1.0, 1.5, 2.0):
        h, v, N = pairing_model(g=g)
        fk = fock_matrix(h, v, N)
        e_ref = reference_energy(h, v, N)
        result = ccd(fk, v, N)
        e_fci, _ = fci_energy(h, v, N)
        print(f"{g:6.2f} {e_ref:9.5f} {e_ref + mp2_energy(fk, v, N):12.8f} "
              f"{e_ref + result['energy']:12.8f} {e_fci:12.8f} "
              f"{e_ref + result['energy'] - e_fci:12.2e} "
              f"{result['iterations']:11d}")
    print()
    print("CCD is essentially exact for this model, to five or six digits")
    print("even at g = 2 where the perturbation series of chapter 9 has")
    print("long since diverged.  The reason is that the pairing interaction")
    print("moves pairs and nothing else: the exact ground state lives in the")
    print("seniority-zero space, and the doubles cluster operator can reach")
    print("all of it.  What CCD misses is only the disconnected quadruples")
    print("that T_2^2 does not generate correctly, which for four particles")
    print("in four levels is very little.")


def demo_orders():
    print("=" * 74)
    print("3. What CCD contains, order by order")
    print("=" * 74)
    print("Expanding the CCD amplitude equation in powers of the interaction")
    print("gives back perturbation theory: the first amplitude is the MP2")
    print("one, the second generates the third-order energy, and so on.  But")
    print("the iteration does not stop, so CCD contains contributions at")
    print("every order.")
    print()
    for g in (0.5, 1.0):
        h, v, N = pairing_model(g=g)
        fk = fock_matrix(h, v, N)
        e_ref = reference_energy(h, v, N)
        energies, _ = ccd_order_by_order(fk, v, N, order=10)
        full = ccd(fk, v, N)
        e_fci, _ = fci_energy(h, v, N)
        print(f"g = {g}:  E_ref = {e_ref:.5f}, "
              f"CCD = {e_ref + full['energy']:.8f}, "
              f"FCI = {e_fci:.8f}")
        print(f"   {'order':>7s} {'contribution':>16s} {'partial sum':>16s}")
        partial = 0.0
        for k, e in enumerate(energies[:8], start=2):
            partial += e
            print(f"   {k:7d} {e:16.10f} {e_ref + partial:16.10f}")
        print()
    print("The order-two and order-three entries are the MP2 and MP3")
    print("energies -- in the Hartree-Fock partition, so they differ from")
    print("the numbers of chapter 9, which used the bare one-body operator")
    print("as H_0.  Beyond that the two schemes part company: the")
    print("perturbation series adds one order at a time and eventually")
    print("diverges, while CCD sums this same expansion to convergence.")


def demo_ccsd():
    print("=" * 74)
    print("4. When singles matter: the pairing plus particle-hole model")
    print("=" * 74)
    print("The particle-hole term of chapter 7 breaks pairs, so the")
    print("reference couples to 1p-1h states and the singles amplitudes are")
    print("no longer zero.  This is where CCSD earns its name.")
    print()
    print(f"{'g':>5s} {'f':>6s} {'E_ref':>10s} {'CCD':>13s} {'CCSD':>13s} "
          f"{'FCI':>13s} {'max |t1|':>10s}")
    for g, ph in ((0.5, 0.025), (0.7, 0.05), (1.0, 0.05), (0.5, 0.25),
                  (1.0, 0.5)):
        h, v, N = pairing_ph_model(g=g, f=ph)
        fk = fock_matrix(h, v, N)
        e_ref = reference_energy(h, v, N)
        d = ccd(fk, v, N, mixing=0.3)
        s = ccsd(fk, v, N, mixing=0.3)
        e_fci, _ = fci_energy(h, v, N)
        print(f"{g:5.2f} {ph:6.3f} {e_ref:10.5f} {e_ref + d['energy']:13.8f} "
              f"{e_ref + s['energy']:13.8f} {e_fci:13.8f} "
              f"{np.abs(s['t1']).max():10.2e}")
    print()
    print("The singles are small when f is small, and they matter more as")
    print("the particle-hole term grows.  Note that the reference here is")
    print("*not* the Hartree-Fock determinant -- we have used the bare")
    print("oscillator basis -- so Brillouin's theorem does not hold and the")
    print("singles have real work to do.  In a Hartree-Fock basis they would")
    print("start at second order in the interaction instead of first.")


def demo_extensivity():
    print("=" * 74)
    print("5. Size extensivity")
    print("=" * 74)
    print("Chapter 5 found that truncated configuration interaction is not")
    print("size consistent, and chapter 9 that Rayleigh-Schrodinger")
    print("perturbation theory is.  Coupled cluster is size extensive by")
    print("construction, because exp(T_A + T_B) = exp(T_A) exp(T_B) for")
    print("non-interacting subsystems.  The test:")
    print()
    print(f"{'g':>6s} {'2 x one copy':>16s} {'two copies':>16s} "
          f"{'error':>11s}")
    for g in (0.25, 0.5, 1.0):
        h1, v1, n1 = independent_copies(subsystems=1, g=g)
        h2, v2, n2 = independent_copies(subsystems=2, g=g)
        e1 = ccd(fock_matrix(h1, v1, n1), v1, n1)["energy"]
        e2 = ccd(fock_matrix(h2, v2, n2), v2, n2)["energy"]
        print(f"{g:6.2f} {2 * e1:16.10f} {e2:16.10f} {e2 - 2 * e1:11.2e}")
    print()
    print("Exact to machine precision, at every coupling.  This is the")
    print("property that truncated configuration interaction cannot have and")
    print("that coupled cluster gets for free from the exponential.")


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
        have = True
    except Exception:
        have = False
    if not have:
        print("(the other chapters' programs are not importable here)")
        return
    fock_rpa = rpamod.FockSpace(4)
    print(f"{'g':>6s} {'HF':>9s} {'MBPT2':>10s} {'MBPT3':>10s} "
          f"{'CID':>10s} {'RPA':>10s} {'CCD':>11s} {'FCI':>11s}")
    for g in (0.25, 0.5, 1.0, 1.5, 2.0):
        h, v, N = pairing_model(g=g)
        fk = fock_matrix(h, v, N)
        e_ref = reference_energy(h, v, N)
        cc = ccd(fk, v, N)
        part, model = mbptmod.pairing_partition(g=g)
        rs = mbptmod.rayleigh_schrodinger(part, order=3)
        r = rpamod.tda_rpa(fock_rpa, 4, g, 0.0)
        e_fci, _ = fci_energy(h, v, N)
        print(f"{g:6.2f} {e_ref:9.5f} "
              f"{part.reference_energy + rs[1]:10.5f} "
              f"{part.reference_energy + rs[1] + rs[2]:10.5f} "
              f"{model.energies(2)[0]:10.5f} "
              f"{r['hf'] + r['ecorr']:10.5f} "
              f"{e_ref + cc['energy']:11.7f} {e_fci:11.7f}")
    print()
    print("The errors:")
    print(f"{'g':>6s} {'MBPT2':>11s} {'MBPT3':>11s} {'CID':>11s} "
          f"{'RPA':>11s} {'CCD':>11s}")
    for g in (0.25, 0.5, 1.0, 1.5, 2.0):
        h, v, N = pairing_model(g=g)
        fk = fock_matrix(h, v, N)
        e_ref = reference_energy(h, v, N)
        cc = ccd(fk, v, N)
        part, model = mbptmod.pairing_partition(g=g)
        rs = mbptmod.rayleigh_schrodinger(part, order=3)
        r = rpamod.tda_rpa(fock_rpa, 4, g, 0.0)
        e_fci, _ = fci_energy(h, v, N)
        print(f"{g:6.2f} "
              f"{part.reference_energy + rs[1] - e_fci:11.2e} "
              f"{part.reference_energy + rs[1] + rs[2] - e_fci:11.2e} "
              f"{model.energies(2)[0] - e_fci:11.2e} "
              f"{r['hf'] + r['ecorr'] - e_fci:11.2e} "
              f"{e_ref + cc['energy'] - e_fci:11.2e}")
    print()
    print("Coupled cluster with doubles is better than everything else by")
    print("three to five orders of magnitude, at every coupling, and it does")
    print("not deteriorate as the interaction grows.  It is not variational")
    print("-- the errors have no fixed sign -- and it is not exact, but for")
    print("this model it is as close to exact as makes no practical")
    print("difference, at a cost that would scale as n^6 rather than")
    print("exponentially.")


def _demo():
    for f in (demo_validation, demo_pairing, demo_orders, demo_ccsd,
              demo_extensivity, demo_comparison):
        f()
        print()


if __name__ == "__main__":
    _demo()
