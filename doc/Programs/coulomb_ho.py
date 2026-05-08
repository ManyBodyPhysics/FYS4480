"""
coulomb_ho.py
=============
Python translation of QD_Coulomb.cpp / Coulomb_Functions.cpp / Coulomb_Functions.hpp

Computes the two-body Coulomb matrix element in a 2-D harmonic-oscillator
basis (quantum dot):

    < ni, mi ; nj, mj || V || nk, mk ; nl, ml >

where n  = radial quantum number  (0, 1, 2, …)
      ml = angular-momentum projection (integer, positive or negative)
      hw = oscillator energy in whatever units the caller uses

The matrix element is proportional to sqrt(hw), and vanishes unless
    mi + mj = mk + ml   (angular-momentum conservation).

Usage as a script
-----------------
    python coulomb_ho.py  hw  n1 ml1  n2 ml2  n3 ml3  n4 ml4

Usage as a module
-----------------
    from coulomb_ho import Coulomb_HO
    value = Coulomb_HO(hw, n1, ml1, n2, ml2, n3, ml3, n4, ml4)

Notes on the translation
------------------------
* Every C++ function is translated 1-to-1.
* std::lgamma(x) → math.lgamma(x)   (identical for x > 0)
* logfac(n)  uses the same explicit loop as the C++ version for
  exact correspondence; lgamma(n+1) would also work.
* The exchange term (exch) remains zero (the corresponding block is
  commented out in the C++ source).
* All arguments are plain Python ints/floats – no pass-by-reference
  needed.
"""

import math
import sys


# ── helper: log factorial ──────────────────────────────────────────────────────

def logfac(n: int) -> float:
    """Return ln(n!).

    Mirrors the C++ implementation exactly: sum of log(a) for a in [2, n].
    logfac(0) = logfac(1) = 0.0   (empty sum).
    """
    if n < 0:
        raise ValueError(f"{n}: logfac input must be >= 0")
    fac = 0.0
    for a in range(2, n + 1):
        fac += math.log(a)
    return fac


# ── helper: log-ratio of four factorials (negative sum) ───────────────────────

def logratio1(int1: int, int2: int, int3: int, int4: int) -> float:
    """Return  -ln(int1!) - ln(int2!) - ln(int3!) - ln(int4!)."""
    return -(logfac(int1) + logfac(int2) + logfac(int3) + logfac(int4))


# ── helper: -½(G+1)·ln 2 ──────────────────────────────────────────────────────

def logratio2(G: int) -> float:
    """Return  -½(G+1) · ln(2)."""
    return -0.5 * (G + 1) * math.log(2)


# ── helper: normalisation prefactor (product1) ────────────────────────────────

def product1(n1: int, m1: int, n2: int, m2: int,
             n3: int, m3: int, n4: int, m4: int) -> float:
    """Return  exp{ ½[ ln(n1!)+ln(n2!)+ln(n3!)+ln(n4!)
                      - ln((n1+|m1|)!) - ln((n2+|m2|)!)
                      - ln((n3+|m3|)!) - ln((n4+|m4|)!) ] }.

    This is the radial normalisation prefactor of the HO wave functions.
    """
    prod  = logfac(n1) + logfac(n2) + logfac(n3) + logfac(n4)
    prod -= (logfac(n1 + abs(m1)) + logfac(n2 + abs(m2))
           + logfac(n3 + abs(m3)) + logfac(n4 + abs(m4)))
    return math.exp(0.5 * prod)


# ── helper: log of the j-dependent binomial product ───────────────────────────

def logproduct2(n1: int, m1: int, n2: int, m2: int,
                n3: int, m3: int, n4: int, m4: int,
                j1: int, j2: int, j3: int, j4: int) -> float:
    """Return  ln[ C(n1+|m1|, j1+|m1|) · C(n2+|m2|, j2+|m2|)
                 · C(n3+|m3|, j3+|m3|) · C(n4+|m4|, j4+|m4|) ]

    where C(a, b) = a! / (b! · (a-b)!), written in log space as
        ln(a!) - ln((a-j)!) - ln((j+|m|)!)
    """
    arg1 = n1 + abs(m1);  arg2 = n2 + abs(m2)
    arg3 = n3 + abs(m3);  arg4 = n4 + abs(m4)

    narg1 = n1 - j1;  narg2 = n2 - j2
    narg3 = n3 - j3;  narg4 = n4 - j4

    jarg1 = j1 + abs(m1);  jarg2 = j2 + abs(m2)
    jarg3 = j3 + abs(m3);  jarg4 = j4 + abs(m4)

    prod  = (logfac(arg1) + logfac(arg2) + logfac(arg3) + logfac(arg4))
    prod -= (logfac(narg1) + logfac(narg2) + logfac(narg3) + logfac(narg4))
    prod -= (logfac(jarg1) + logfac(jarg2) + logfac(jarg3) + logfac(jarg4))
    return prod


# ── helper: log of the l-dependent binomial product ───────────────────────────

def logproduct3(l1: int, l2: int, l3: int, l4: int,
                g1: int, g2: int, g3: int, g4: int) -> float:
    """Return  ln[ C(g1,l1) · C(g2,l2) · C(g3,l3) · C(g4,l4) ]

    where C(gi, li) = gi! / (li! · (gi-li)!).
    """
    garg1 = g1 - l1;  garg2 = g2 - l2
    garg3 = g3 - l3;  garg4 = g4 - l4

    prod  = (logfac(g1) + logfac(g2) + logfac(g3) + logfac(g4))
    prod -= (logfac(l1) + logfac(l2) + logfac(l3) + logfac(l4))
    prod -= (logfac(garg1) + logfac(garg2) + logfac(garg3) + logfac(garg4))
    return prod


# ── main function ──────────────────────────────────────────────────────────────

def Coulomb_HO(hw: float,
               ni: int, mi: int,
               nj: int, mj: int,
               nk: int, mk: int,
               nl: int, ml: int) -> float:
    """Compute the two-body Coulomb matrix element in a 2-D HO basis.

    Parameters
    ----------
    hw          : oscillator energy (determines the overall scale ∝ sqrt(hw))
    ni, mi      : radial and angular-momentum quantum numbers of state i (bra, left)
    nj, mj      : quantum numbers of state j (bra, right)
    nk, mk      : quantum numbers of state k (ket, left)
    nl, ml      : quantum numbers of state l (ket, right)

    Returns
    -------
    float
        The antisymmetrised two-body matrix element
            < ni,mi ; nj,mj || V || nk,mk ; nl,ml >
        in units where V = e²/(4πε₀) and lengths are in oscillator units.

    Notes
    -----
    * Returns 0 immediately if angular-momentum is not conserved:
          mi + mj ≠ mk + ml
    * The exchange contribution (exch) is currently zero because the
      corresponding code block is commented out in the original C++ source.
    * The inner logproduct2 call passes arguments in the order
          (ni,mi, nj,mj, nl,ml, nk,mk, j1,j2,j3,j4)
      — note nl/ml before nk/mk — reproducing the C++ source exactly.
    """
    # Angular-momentum conservation
    if mi + mj != mk + ml:
        return 0.0

    dir  = 0.0
    exch = 0.0   # exchange term is commented out in C++, stays zero

    # ── outer four-fold loop over j indices ───────────────────────────────────
    for j1 in range(ni + 1):
        for j2 in range(nj + 1):
            for j3 in range(nl + 1):          # note: j3 runs to nl (not nk)
                for j4 in range(nk + 1):      # note: j4 runs to nk (not nl)

                    # ── compute the g values ──────────────────────────────────
                    # Each g combines the j summation index with a
                    # step-function on the angular-momentum projection
                    # that selects the positive/negative part of m.
                    g1 = int(j1 + j4
                             + 0.5 * (abs(mi) + mi)
                             + 0.5 * (abs(mk) - mk))
                    g2 = int(j2 + j3
                             + 0.5 * (abs(mj) + mj)
                             + 0.5 * (abs(ml) - ml))
                    g3 = int(j3 + j2
                             + 0.5 * (abs(ml) + ml)
                             + 0.5 * (abs(mj) - mj))
                    g4 = int(j4 + j1
                             + 0.5 * (abs(mk) + mk)
                             + 0.5 * (abs(mi) - mi))
                    G = g1 + g2 + g3 + g4

                    log_ratio1 = logratio1(j1, j2, j3, j4)
                    # Note: logproduct2 receives (nl,ml,nk,mk) — swapped
                    # relative to the outer function's (nk,mk,nl,ml) order.
                    log_prod2  = logproduct2(ni, mi, nj, mj,
                                            nl, ml, nk, mk,
                                            j1, j2, j3, j4)
                    log_ratio2 = logratio2(G)

                    # ── inner four-fold loop over l indices ───────────────────
                    temp = 0.0
                    for l1 in range(g1 + 1):
                        for l2 in range(g2 + 1):
                            for l3 in range(g3 + 1):
                                for l4 in range(g4 + 1):

                                    # Angular-momentum conservation
                                    # in the l-sum
                                    if l1 + l2 != l3 + l4:
                                        continue

                                    L = l1 + l2 + l3 + l4

                                    # Sign factor:
                                    #   (-1)^(g2 + g3 - l2 - l3)
                                    sign_l = (-2 * ((g2 + g3 - l2 - l3) % 2)
                                              + 1)

                                    log_p3 = logproduct3(l1, l2, l3, l4,
                                                         g1, g2, g3, g4)
                                    temp += (sign_l
                                             * math.exp(
                                                 log_p3
                                                 + math.lgamma(1.0 + 0.5 * L)
                                                 + math.lgamma(0.5 * (G - L + 1.0))
                                             ))

                    # Sign factor: (-1)^(j1 + j2 + j3 + j4)
                    sign_j = -2 * ((j1 + j2 + j3 + j4) % 2) + 1

                    dir += (sign_j
                            * math.exp(log_ratio1 + log_prod2 + log_ratio2)
                            * temp)

    # Multiply by the normalisation prefactor
    # Note: product1 also receives (nl,ml,nk,mk) — matching the C++ source.
    dir *= product1(ni, mi, nj, mj, nl, ml, nk, mk)

    # exchange term is zero (commented out in C++)
    return math.sqrt(hw) * (dir - exch)


# ── command-line interface ─────────────────────────────────────────────────────

def _print_usage():
    print(
        "Usage: python coulomb_ho.py  hw  n1 ml1  n2 ml2  n3 ml3  n4 ml4",
        file=sys.stderr,
    )


if __name__ == "__main__":
    if len(sys.argv) != 10:
        _print_usage()
        sys.exit(1)

    try:
        hw_  = float(sys.argv[1])
        n1_  = int(sys.argv[2]);  ml1_ = int(sys.argv[3])
        n2_  = int(sys.argv[4]);  ml2_ = int(sys.argv[5])
        n3_  = int(sys.argv[6]);  ml3_ = int(sys.argv[7])
        n4_  = int(sys.argv[8]);  ml4_ = int(sys.argv[9])
    except ValueError:
        _print_usage()
        sys.exit(1)

    tbme = Coulomb_HO(hw_, n1_, ml1_, n2_, ml2_, n3_, ml3_, n4_, ml4_)

    print(
        f"< {n1_},{ml1_} ; {n2_},{ml2_} || V || "
        f"{n3_},{ml3_} ; {n4_},{ml4_} > = {tbme:.12f}"
    )
