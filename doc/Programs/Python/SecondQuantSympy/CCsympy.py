import sys

from sympy.physics.secondquant import (AntiSymmetricTensor, wicks,
        F, Fd, NO, evaluate_deltas, substitute_dummies, Commutator,
        simplify_index_permutations, PermutationOperator, contraction)
from sympy import (
    symbols, expand, pprint, Rational, latex, Dummy
)

pretty_dummies_dict = {
    'above': 'cdefgh',
    'below': 'klmno',
    'general': 'pqrstu'
}

def get_CC_operators():
	
	"""
	Returns a tuple (T1,T2) of unique operators.
	"""
	
	i = symbols('i', below_fermi=True, cls=Dummy)
	a = symbols('a', above_fermi=True, cls=Dummy)
	t_ai = AntiSymmetricTensor('t', (a,), (i,))
	ai = NO(Fd(a)*F(i))
	
	i, j = symbols('i,j', below_fermi=True, cls=Dummy)
	a, b = symbols('a,b', above_fermi=True, cls=Dummy)
	t_abij = AntiSymmetricTensor('t', (a, b), (i, j))
	abji = NO(Fd(a)*Fd(b)*F(j)*F(i))

	T1 = t_ai*ai
	T2 = Rational(1, 4)*t_abij*abji
	
	return (T1, T2)

p, q, r, s = symbols('p,q,r,s', cls=Dummy)
#Setup Hamiltonian on normal ordered form
E0 = symbols('Eref', real = True, constant = True) #Reference energy

f    = AntiSymmetricTensor('f', (p,), (q,))
pq   = NO((Fd(p)*F(q)))
Fock = f*pq #F is reserved by sympy

V    = AntiSymmetricTensor('V', (p, q), (r, s))
pqsr = NO(Fd(p)*Fd(q)*F(s)*F(r)) 

HI   = Rational(1, 4)*V*pqsr

HN   = E0+F+HI

i, j, k, l = symbols('i,j,k,l', below_fermi=True)
a, b, c, d = symbols('a,b,c,d', above_fermi=True)

print("Compute CCD energy and amplitude equations term by term")
#BCH: HN + [HN,T] + 0.5*[[HN,T],T] + 1/6 * [[[HN,T],T],T] + 1/24 * [[[[HN,T],T],T],T]

print("zero-order term")
print("<Phi|HN|Phi> = ", latex(wicks(HN, simplify_dummies=True,keep_only_fully_contracted=True)))
eqT2 = wicks(NO(Fd(i)*Fd(j)*F(b)*F(a))*HN, simplify_dummies=True, keep_only_fully_contracted=True, simplify_kronecker_deltas=True)
P    = PermutationOperator
eqT2 = simplify_index_permutations(eqT2, [P(a, b), P(i, j)])
print("<Phi|HN|Phi_ij^ab = ", latex(eqT2))
print()


C = Commutator
T1, T2 = get_CC_operators()
T = T2

print("[HN,T]-term")
comm1 = wicks(C(HN, T))
comm1 = evaluate_deltas(comm1)
comm1 = substitute_dummies(comm1,new_indices=True, pretty_indices=pretty_dummies_dict)
print("<Phi|[HN,T]|Phi> = ", latex(wicks(comm1, simplify_dummies=True,keep_only_fully_contracted=True)))
eqT2 = wicks(NO(Fd(i)*Fd(j)*F(b)*F(a))*comm1, simplify_dummies=True, keep_only_fully_contracted=True, simplify_kronecker_deltas=True)
P    = PermutationOperator
eqT2 = simplify_index_permutations(eqT2, [P(a, b), P(i, j)])
print("<Phi|[HN,T]|Phi_ij^ab = ", latex(eqT2))
print()

T1, T2 = get_CC_operators()
T = T2
print("[[HN,T],T]-term")
comm2 = wicks(C(comm1, T))
comm2 = evaluate_deltas(comm2)/2
comm2 = substitute_dummies(comm2,new_indices=True, pretty_indices=pretty_dummies_dict)
print("0.5*<Phi|[[HN,T],T]|Phi> = ", latex(wicks(comm2, simplify_dummies=True,keep_only_fully_contracted=True)))
eqT2 = wicks(NO(Fd(i)*Fd(j)*F(b)*F(a))*comm2, simplify_dummies=True, keep_only_fully_contracted=True, simplify_kronecker_deltas=True)
P    = PermutationOperator
eqT2 = simplify_index_permutations(eqT2, [P(a, b), P(i, j)])
print("0.5*<Phi|[[HN,T],T]|Phi_ij^ab = ", latex(eqT2))
print()

T1, T2 = get_CC_operators()
T = T2
print("[[[HN,T],T],T]-term")
comm3 = wicks(C(comm2, T))
comm3 = evaluate_deltas(comm3)/6
comm3 = substitute_dummies(comm3,new_indices=True, pretty_indices=pretty_dummies_dict)
print("1/6*<Phi|[[[HN,T],T],T]|Phi> = ", latex(wicks(comm3, simplify_dummies=True,keep_only_fully_contracted=True)))
eqT2 = wicks(NO(Fd(i)*Fd(j)*F(b)*F(a))*comm3, simplify_dummies=True, keep_only_fully_contracted=True, simplify_kronecker_deltas=True)
P    = PermutationOperator
eqT2 = simplify_index_permutations(eqT2, [P(a, b), P(i, j)])
print("1/6*<Phi|[[[HN,T],T],T]|Phi_ij^ab = ", latex(eqT2))
print()




