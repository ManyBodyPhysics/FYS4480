"""
Wick's theorem, verified by brute force.

Companion code to chapter 3 of *Quantum mechanics for Many-particle Systems*.

A vacuum expectation value of a string of fermion creation and annihilation
operators can be evaluated in two entirely different ways:

  1. by repeatedly using the anticommutation relation
     {a_p, a_q^dagger} = delta_pq to push annihilation operators to the right
     until they hit the vacuum -- the elementary but exponentially branching
     route of the worked examples in the text;

  2. by Wick's theorem, as a signed sum over all ways of pairing the operators
     into contractions -- one term per perfect matching, and only the fully
     contracted terms contribute.

This program implements both and checks that they agree.  That check *is* the
content of the theorem, and it is worth running once before trusting any
hand calculation.

Operators are written as (index, dagger) pairs, so ('p', True) is
a_p^dagger and ('p', False) is a_p.  Indices are symbolic: the result of a
vacuum expectation value is a linear combination of products of Kronecker
deltas, represented as a dictionary mapping a frozenset of delta pairs to a
coefficient.

Author: Morten Hjorth-Jensen
"""

from itertools import combinations


# ---------------------------------------------------------------------------
#  Symbolic products of Kronecker deltas
# ---------------------------------------------------------------------------
class DeltaSum:
    """A linear combination of products of Kronecker deltas.

    Stored as {term: coefficient}, where a term is a sorted tuple of pairs
    (p, q) standing for delta_pq.  Deltas with equal indices are dropped,
    since delta_pp = 1.
    """

    def __init__(self, terms=None):
        self.terms = dict(terms) if terms else {}

    @staticmethod
    def zero():
        return DeltaSum()

    @staticmethod
    def one():
        return DeltaSum({(): 1})

    @staticmethod
    def delta(p, q):
        if p == q:
            return DeltaSum.one()
        return DeltaSum({((p, q) if p < q else (q, p),): 1})

    # -- algebra ---------------------------------------------------------
    def __add__(self, other):
        out = dict(self.terms)
        for term, coeff in other.terms.items():
            out[term] = out.get(term, 0) + coeff
            if out[term] == 0:
                del out[term]
        return DeltaSum(out)

    def __mul__(self, other):
        if isinstance(other, int):
            return DeltaSum({t: c * other for t, c in self.terms.items()
                             if c * other != 0})
        out = {}
        for t1, c1 in self.terms.items():
            for t2, c2 in other.terms.items():
                term = tuple(sorted(set(t1) | set(t2)))
                out[term] = out.get(term, 0) + c1 * c2
                if out[term] == 0:
                    del out[term]
        return DeltaSum(out)

    __rmul__ = __mul__

    def __eq__(self, other):
        return self.terms == other.terms

    def __bool__(self):
        return bool(self.terms)

    def __repr__(self):
        if not self.terms:
            return "0"
        pieces = []
        for term, coeff in sorted(self.terms.items()):
            body = "".join(f"d({p}{q})" for p, q in term) or "1"
            sign = "+" if coeff > 0 else "-"
            mag = "" if abs(coeff) == 1 else str(abs(coeff))
            pieces.append(f"{sign} {mag}{body}")
        text = " ".join(pieces)
        return text[2:] if text.startswith("+ ") else text

    def substitute(self, values):
        """Evaluate with concrete index values, returning a number."""
        total = 0
        for term, coeff in self.terms.items():
            product = 1
            for p, q in term:
                product *= 1 if values.get(p, p) == values.get(q, q) else 0
                if product == 0:
                    break
            total += coeff * product
        return total


# ---------------------------------------------------------------------------
#  Route 1: brute-force anticommutation
# ---------------------------------------------------------------------------
_call_count = [0]


def vacuum_expectation_bruteforce(ops):
    """<0| ops |0> by repeated use of {a_p, a_q^dagger} = delta_pq.

    The recursion is the one carried out by hand in the text: find the
    leftmost place where an annihilation operator stands immediately to the
    left of a creation operator, and replace

        a_p a_q^dagger  ->  delta_pq - a_q^dagger a_p .

    When no such place remains, the string is normal-ordered and its vacuum
    expectation value vanishes unless the string is empty.
    """
    ops = tuple(ops)
    _call_count[0] += 1

    if len(ops) == 0:
        return DeltaSum.one()
    if len(ops) % 2 == 1:              # an odd string can never pair up
        return DeltaSum.zero()
    if ops[0][1]:                      # starts with a creation operator
        return DeltaSum.zero()         #   -> kills the bra vacuum
    if not ops[-1][1]:                 # ends with an annihilation operator
        return DeltaSum.zero()         #   -> kills the ket vacuum

    for k in range(len(ops) - 1):
        left, right = ops[k], ops[k + 1]
        if not left[1] and right[1]:               # a_p a_q^dagger
            rest = ops[:k] + ops[k + 2:]
            swapped = ops[:k] + (right, left) + ops[k + 2:]
            return (DeltaSum.delta(left[0], right[0])
                    * vacuum_expectation_bruteforce(rest)
                    + vacuum_expectation_bruteforce(swapped) * (-1))
    return DeltaSum.zero()             # already normal-ordered, and non-empty


# ---------------------------------------------------------------------------
#  Route 2: Wick's theorem
# ---------------------------------------------------------------------------
def contraction(a, b):
    """The elementary contraction of two operators, <0| a b |0>.

    Only a_p a_q^dagger survives; the other three combinations vanish, which
    is the statement that the vacuum contains no particles.
    """
    if (not a[1]) and b[1]:
        return DeltaSum.delta(a[0], b[0])
    return DeltaSum.zero()


def perfect_matchings(indices):
    """All ways of pairing up a list of positions, as lists of index pairs."""
    if not indices:
        yield []
        return
    first, rest = indices[0], indices[1:]
    for k in range(len(rest)):
        partner = rest[k]
        remainder = rest[:k] + rest[k + 1:]
        for tail in perfect_matchings(remainder):
            yield [(first, partner)] + tail


def matching_sign(pairs):
    """The sign of a set of contractions: (-1)^(number of crossings).

    Each crossing of two contraction lines costs one interchange of adjacent
    fermion operators, hence one minus sign.
    """
    crossings = 0
    for (a, b), (c, d) in combinations(pairs, 2):
        if a < c < b < d or c < a < d < b:
            crossings += 1
    return (-1) ** crossings


def vacuum_expectation_wick(ops):
    """<0| ops |0> from Wick's theorem: a signed sum over perfect matchings.

    Only fully contracted terms survive the vacuum expectation value, so the
    sum runs over the (M-1)!! ways of pairing M operators.
    """
    ops = tuple(ops)
    if len(ops) % 2 == 1:
        return DeltaSum.zero()

    total = DeltaSum.zero()
    for pairs in perfect_matchings(list(range(len(ops)))):
        value = DeltaSum.one()
        for i, j in pairs:
            value = value * contraction(ops[i], ops[j])
            if not value:
                break
        if value:
            total = total + value * matching_sign(pairs)
    return total


def bruteforce_cost(ops):
    """Number of recursive calls the anticommutation route actually makes."""
    _call_count[0] = 0
    vacuum_expectation_bruteforce(ops)
    return _call_count[0]


def double_factorial(m):
    """(m-1)!! -- the number of perfect matchings of m objects."""
    out, k = 1, m - 1
    while k > 1:
        out *= k
        k -= 2
    return out


# ---------------------------------------------------------------------------
#  Convenience constructors
# ---------------------------------------------------------------------------
def c(index):
    """A creation operator a_index^dagger."""
    return (index, True)


def a(index):
    """An annihilation operator a_index."""
    return (index, False)


# ---------------------------------------------------------------------------
def _demo():
    print("=" * 74)
    print("1. Two operators")
    print("=" * 74)
    for ops, label in (((a("p"), c("q")), "<0| a_p a_q^+ |0>"),
                       ((a("p"), a("q")), "<0| a_p a_q |0>"),
                       ((c("p"), a("q")), "<0| a_p^+ a_q |0>"),
                       ((c("p"), c("q")), "<0| a_p^+ a_q^+ |0>")):
        brute = vacuum_expectation_bruteforce(ops)
        wick = vacuum_expectation_wick(ops)
        print(f"{label:24s} = {str(brute):28s} "
              f"(Wick: {wick},  agree: {brute == wick})")

    print()
    print("=" * 74)
    print("2. Three operators: an odd string always vanishes")
    print("=" * 74)
    for ops, label in (((a("p"), c("q"), c("r")), "<0| a_p a_q^+ a_r^+ |0>"),
                       ((a("p"), a("q"), c("r")), "<0| a_p a_q a_r^+ |0>")):
        brute = vacuum_expectation_bruteforce(ops)
        wick = vacuum_expectation_wick(ops)
        print(f"{label:26s} = {str(brute):6s} "
              f"(Wick: {wick},  agree: {brute == wick})")
    print("With an odd number of operators at least one is always left")
    print("uncontracted, and <0|a|0> = <0|a^+|0> = 0.")

    print()
    print("=" * 74)
    print("3. Four operators: the overlap <rs|pq>")
    print("=" * 74)
    ops = (a("s"), a("r"), c("p"), c("q"))
    brute = vacuum_expectation_bruteforce(ops)
    wick = vacuum_expectation_wick(ops)
    print("<0| a_s a_r a_p^+ a_q^+ |0>")
    print(f"  by anticommutation : {brute}")
    print(f"  by Wick's theorem  : {wick}")
    print(f"  agree              : {brute == wick}")
    print("  which is <rs|pq> = d(rp)d(sq) - d(sp)d(rq), the antisymmetrised")
    print("  overlap of two two-particle states.")

    print()
    print("=" * 74)
    print("4. Six and eight operators")
    print("=" * 74)
    for names in (("f", "e", "d", "c", "b", "a"),
                  ("h", "g", "f", "e", "d", "c", "b", "a")):
        half = len(names) // 2
        ops = tuple(a(n) for n in names[:half]) + \
              tuple(c(n) for n in names[half:])
        brute = vacuum_expectation_bruteforce(ops)
        wick = vacuum_expectation_wick(ops)
        print(f"M = {len(names)}: {len(brute.terms):3d} terms, "
              f"(M-1)!! = {double_factorial(len(names)):3d} matchings, "
              f"agree: {brute == wick}")

    print()
    print("=" * 74)
    print("5. The number operator, <12| N |12>")
    print("=" * 74)
    # |12> = a_1^+ a_2^+ |0>,  N = sum_i a_i^+ a_i
    total = 0
    for i in ("1", "2", "3"):
        ops = (a("2"), a("1"), c(i), a(i), c("1"), c("2"))
        value = vacuum_expectation_wick(ops).substitute({})
        total += value
        print(f"  i = {i}: <12| a_{i}^+ a_{i} |12> = {value}")
    print(f"  sum over i          = {total}   (= N, the particle number)")

    print()
    print("=" * 74)
    print("6. The two-body interaction, <ij| H_I |ij>")
    print("=" * 74)
    print("H_I = (1/2) sum_pqrs <pq|v|rs> a_p^+ a_q^+ a_s a_r, and")
    print("|ij> = a_i^+ a_j^+ |0>.  We evaluate the operator string for each")
    print("(p,q,r,s) and collect the surviving matrix elements.")
    print()
    surviving = {}
    labels = ("i", "j")
    for p in labels:
        for q in labels:
            for r in labels:
                for s in labels:
                    ops = (a("j"), a("i"), c(p), c(q), a(s), a(r),
                           c("i"), c("j"))
                    value = vacuum_expectation_wick(ops).substitute({})
                    if value:
                        surviving[(p, q, r, s)] = value
    for key, value in sorted(surviving.items()):
        p, q, r, s = key
        print(f"  <{p}{q}|v|{r}{s}> coefficient {value:+d}")
    direct = sum(v for k, v in surviving.items() if k in (("i", "j", "i", "j"),
                                                          ("j", "i", "j", "i")))
    exchange = sum(v for k, v in surviving.items()
                   if k in (("i", "j", "j", "i"), ("j", "i", "i", "j")))
    print()
    print(f"  half the direct terms   : {direct // 2:+d} <ij|v|ij>")
    print(f"  half the exchange terms : {exchange // 2:+d} <ij|v|ji>")
    print("  so <ij|H_I|ij> = <ij|v|ij> - <ij|v|ji> = <ij|v|ij>_AS,")
    print("  the antisymmetrised matrix element of the text.")

    print()
    print("=" * 74)
    print("7. Cost of the two routes")
    print("=" * 74)
    print(f"{'M':>4s} {'matchings (M-1)!!':>20s} "
          f"{'recursive calls':>18s} {'surviving terms':>18s}")
    alphabet = "abcdefghijklmnop"
    for m in (2, 4, 6, 8, 10):
        half = m // 2
        names = alphabet[:m]
        ops = tuple(a(n) for n in names[:half]) + \
              tuple(c(n) for n in names[half:])
        calls = bruteforce_cost(ops)
        terms = len(vacuum_expectation_wick(ops).terms)
        print(f"{m:4d} {double_factorial(m):20d} {calls:18d} {terms:18d}")
    print()
    print("Both routes grow quickly, but Wick's theorem grows in a way we")
    print("can reason about: the terms are enumerated in advance, most vanish")
    print("on inspection because a contraction is zero, and the survivors can")
    print("be written down directly.  The anticommutation recursion offers no")
    print("such structure -- which is why it is unusable beyond four or six")
    print("operators by hand.")


if __name__ == "__main__":
    _demo()
