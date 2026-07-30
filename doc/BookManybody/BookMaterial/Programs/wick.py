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
     contracted terms contribute;

  3. by Wick's generalised theorem, when the string comes to us as a product
     of groups that are separately normal-ordered -- the bra, the operator
     and the ket -- in which case only the contractions joining two different
     groups survive.

This program implements all three and checks that they agree.  That check
*is* the content of the theorems, and it is worth running once before
trusting any hand calculation.  The generalised theorem is checked twice
over: under a vacuum expectation value, and as an identity between operators
represented as matrices in a small Fock space.

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
#  Route 3: Wick's generalised theorem
# ---------------------------------------------------------------------------
#  The generalised theorem applies to a product of groups, each of which is
#  already normal-ordered,
#
#      N[A_1 A_2 ...] N[B_1 B_2 ...] N[C_1 C_2 ...]
#          = N[A_1 A_2 ... B_1 B_2 ... C_1 C_2 ...]
#            + sum over all contractions BETWEEN DIFFERENT GROUPS.
#
#  Contractions inside a group are absent because a group in normal-ordered
#  form has all its creation operators on the left, and a non-vanishing
#  contraction needs an annihilation operator to the left of a creation
#  operator.  The functions below check that statement in two ways: as an
#  identity between vacuum expectation values, and -- more stringently -- as
#  an identity between operators, represented as matrices in a finite Fock
#  space.
# ---------------------------------------------------------------------------
def normal_order(ops):
    """Bring a string into normal-ordered form.

    Returns (sign, reordered), where reordered has all creation operators
    first, the relative order inside each class being preserved, and sign is
    the (-1)^kappa of the definition of the normal-ordered product: one
    factor of -1 for every interchange of neighbouring fermion operators
    needed to get there.
    """
    ops = tuple(ops)
    creators = [o for o in ops if o[1]]
    annihilators = [o for o in ops if not o[1]]
    swaps = sum(1 for i, o in enumerate(ops) if not o[1]
                for p in ops[i + 1:] if p[1])
    return (-1) ** swaps, tuple(creators + annihilators)


def is_normal_ordered(ops):
    """True if no annihilation operator stands to the left of a creator."""
    seen_annihilator = False
    for _, dagger in ops:
        if dagger and seen_annihilator:
            return False
        if not dagger:
            seen_annihilator = True
    return True


def _inversions(sequence):
    """Number of out-of-order pairs, i.e. the parity of a permutation."""
    return sum(1 for i, x in enumerate(sequence)
               for y in sequence[i + 1:] if y < x)


def contraction_sign(pairs, n_ops):
    """Sign of a (possibly partial) set of contractions.

    The rule is the one of the text: bring the contracted operators next to
    each other, leaving the uncontracted ones in their original relative
    order, and count the interchanges.  For a fully contracted term this
    reduces to (-1)^(number of crossings), matching_sign() above.
    """
    used = {i for pair in pairs for i in pair}
    order = [i for pair in pairs for i in pair]
    order += [k for k in range(n_ops) if k not in used]
    return (-1) ** _inversions(order)


def contraction_sets(n_ops, allowed=None):
    """All sets of disjoint contractions, from the empty set upwards.

    ``allowed(i, j)`` filters the admissible pairs; with allowed=None every
    pair i < j is admissible, which is ordinary Wick.  Passing a filter that
    rejects pairs inside the same group gives the generalised theorem.
    """
    def build(decided):
        free = [k for k in range(n_ops) if k not in decided]
        if not free:
            yield []
            return
        first, rest = free[0], free[1:]
        for tail in build(decided | {first}):       # first stays uncontracted
            yield tail
        for partner in rest:                        # or it is contracted
            if allowed and not allowed(first, partner):
                continue
            for tail in build(decided | {first, partner}):
                yield [(first, partner)] + tail
    return build(frozenset())


def _flatten(groups):
    """Concatenate the groups and record which group each operator came from."""
    ops, labels = [], []
    for number, group in enumerate(groups):
        for operator in group:
            ops.append(operator)
            labels.append(number)
    return tuple(ops), labels


def generalised_wick_vev(groups):
    """<0| N[A...] N[B...] ... |0> from the generalised theorem.

    Only the fully contracted terms survive the vacuum expectation value, and
    every contraction line must join two different groups.
    """
    ops, labels = _flatten(groups)
    n_ops = len(ops)
    total = DeltaSum.zero()
    for pairs in contraction_sets(n_ops,
                                  lambda i, j: labels[i] != labels[j]):
        if 2 * len(pairs) != n_ops:            # not fully contracted
            continue
        value = DeltaSum.one()
        for i, j in pairs:
            value = value * contraction(ops[i], ops[j])
            if not value:
                break
        if value:
            total = total + value * contraction_sign(pairs, n_ops)
    return total


def internal_contractions(group):
    """The contractions of a group with itself, all of which should vanish."""
    return [contraction(group[i], group[j])
            for i in range(len(group)) for j in range(i + 1, len(group))]


def count_matchings(sizes):
    """(all perfect matchings, those with no line inside a group)."""
    labels = [number for number, size in enumerate(sizes)
              for _ in range(size)]
    n_ops = len(labels)
    every = inter = 0
    for pairs in perfect_matchings(list(range(n_ops))):
        every += 1
        if all(labels[i] != labels[j] for i, j in pairs):
            inter += 1
    return every, inter


# ---------------------------------------------------------------------------
#  A finite Fock space, so that the theorem can be checked as an operator
#  identity and not only under a vacuum expectation value
# ---------------------------------------------------------------------------
class FockSpace:
    """Fermion operators as matrices on the 2^n states of a Fock space.

    Basis states are bit patterns, orbital p occupying bit p, and the
    Jordan-Wigner sign (-1)^(number of occupied orbitals below p) keeps the
    anticommutation relations exact.
    """

    def __init__(self, n_orbitals):
        self.n_orbitals = n_orbitals
        self.dim = 1 << n_orbitals

    def annihilate(self, p):
        import numpy as np
        matrix = np.zeros((self.dim, self.dim))
        for state in range(self.dim):
            if state & (1 << p):
                sign = (-1) ** bin(state & ((1 << p) - 1)).count("1")
                matrix[state ^ (1 << p), state] = sign
        return matrix

    def create(self, p):
        return self.annihilate(p).T

    def operator(self, op):
        index, dagger = op
        return self.create(index) if dagger else self.annihilate(index)

    def string(self, ops):
        import numpy as np
        matrix = np.eye(self.dim)
        for op in ops:
            matrix = matrix @ self.operator(op)
        return matrix

    def normal_product(self, ops):
        """N[ops] as a matrix: reorder, and carry the sign along."""
        sign, reordered = normal_order(ops)
        return sign * self.string(reordered)

    def check_generalised(self, groups):
        """max |LHS - RHS| of the generalised theorem, as an operator identity."""
        import numpy as np
        lhs = np.eye(self.dim)
        for group in groups:
            lhs = lhs @ self.normal_product(group)

        ops, labels = _flatten(groups)
        n_ops = len(ops)
        rhs = np.zeros((self.dim, self.dim))
        for pairs in contraction_sets(n_ops,
                                      lambda i, j: labels[i] != labels[j]):
            value = 1
            for i, j in pairs:
                value *= contraction(ops[i], ops[j]).substitute({})
                if value == 0:
                    break
            if value == 0:
                continue
            used = {i for pair in pairs for i in pair}
            rest = [ops[k] for k in range(n_ops) if k not in used]
            remainder = (np.eye(self.dim) if not rest
                         else self.normal_product(rest))
            rhs = rhs + contraction_sign(pairs, n_ops) * value * remainder
        return float(np.abs(lhs - rhs).max())


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

    print()
    print("=" * 74)
    print("8. Wick's generalised theorem")
    print("=" * 74)
    one_body = [[a("j"), a("i")], [c("p"), a("q")], [c("k"), c("l")]]
    two_body = [[a("j"), a("i")], [c("p"), c("q"), a("s"), a("r")],
                [c("k"), c("l")]]

    print("The groups are the bra, the operator and the ket, each of them")
    print("already normal-ordered.  Their internal contractions therefore")
    print("vanish one by one:")
    for name, groups in (("one-body", one_body), ("two-body", two_body)):
        ordered = all(is_normal_ordered(g) for g in groups)
        internal = [x for g in groups for x in internal_contractions(g) if x]
        print(f"  {name}: groups normal-ordered: {ordered}, "
              f"non-vanishing internal contractions: {len(internal)}")

    print()
    print("Only lines joining different groups are left, and there are far")
    print("fewer of them than of pairings in general:")
    print(f"{'groups':>16s} {'M':>4s} {'(M-1)!! matchings':>19s} "
          f"{'between groups':>16s}")
    for name, groups in (("one-body", one_body), ("two-body", two_body)):
        sizes = [len(g) for g in groups]
        every, inter = count_matchings(sizes)
        label = "+".join(str(s) for s in sizes)
        print(f"{label:>16s} {sum(sizes):4d} {every:19d} {inter:16d}")

    print()
    print("and the three routes agree on the matrix elements themselves:")
    for name, groups in (("one-body", one_body), ("two-body", two_body)):
        ops, _ = _flatten(groups)
        brute = vacuum_expectation_bruteforce(ops)
        full = vacuum_expectation_wick(ops)
        general = generalised_wick_vev(groups)
        print(f"  {name}: anticommutation == Wick: {brute == full}, "
              f"Wick == generalised: {full == general}")
        print(f"     {general}")

    print()
    print("=" * 74)
    print("9. The generalised theorem as an operator identity")
    print("=" * 74)
    print("The checks above hold under a vacuum expectation value.  The")
    print("theorem is stronger than that: it is an identity between")
    print("operators.  We verify it by representing the operators as")
    print("matrices in a small Fock space and comparing")
    print("  N[A...] N[B...] ...   against   N[A...B...] + all cross terms.")
    print()
    space = FockSpace(4)
    cases = [
        ("N[a_1 a_0] N[a_0^+ a_2] N[a_2^+ a_3^+]",
         [[a(1), a(0)], [c(0), a(2)], [c(2), c(3)]]),
        ("N[a_1 a_0] N[a_0^+ a_1^+ a_3 a_2] N[a_2^+ a_3^+]",
         [[a(1), a(0)], [c(0), c(1), a(3), a(2)], [c(2), c(3)]]),
        ("N[a_0^+ a_1] N[a_1^+ a_0]",
         [[c(0), a(1)], [c(1), a(0)]]),
        ("N[a_0^+ a_1^+ a_1 a_0] N[a_2^+ a_3] N[a_3^+ a_2]",
         [[c(0), c(1), a(1), a(0)], [c(2), a(3)], [c(3), a(2)]]),
    ]
    for label, groups in cases:
        print(f"  max |LHS - RHS| = {space.check_generalised(groups):.1e}"
              f"   for {label}")
    print()
    print("Ordinary Wick's theorem is the special case in which every group")
    print("holds a single operator, so that every contraction is allowed:")
    single = [a(0), c(1), a(2), c(0), c(2), a(1)]
    residual = space.check_generalised([[op] for op in single])
    print(f"  max |LHS - RHS| = {residual:.1e}   for a string of "
          f"{len(single)} operators")


if __name__ == "__main__":
    _demo()
