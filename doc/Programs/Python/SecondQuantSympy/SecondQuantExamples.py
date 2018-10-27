import sys

from sympy.physics.secondquant import (AntiSymmetricTensor, wicks,
        F, Fd, NO, evaluate_deltas, substitute_dummies, Commutator,
        simplify_index_permutations, PermutationOperator, contraction)
from sympy import (
    symbols, expand, pprint, Rational, latex, Dummy
)

#Link to Second quant documentation: https://docs.sympy.org/latest/modules/physics/secondquant.html


pretty_dummies_dict = {
    'above': 'abcdefgh',
    'below': 'ijklmno',
    'general': 'pqrstu'
}


p, q, r, s = symbols('p,q,r,s', cls=Dummy)

#Setup creation and annihilation operators
ap_dagger  = Fd(p)
aq		   = F(q)
#Perform a contraction
contr      = evaluate_deltas(contraction(ap_dagger,aq))
print("Example outputs")
print()
print("contraction(a_p^\dagger a_q): ", latex(contr))
print()

#Setup Hamiltonian, not on normal order form
h  = AntiSymmetricTensor('h', (p,), (q,))
pq = ap_dagger*aq

V    = AntiSymmetricTensor('V', (p, q), (r, s))
pqsr = Fd(p)*Fd(q)*F(s)*F(r)

H0 = h*pq
HI = Rational(1, 4)*V*pqsr

H = H0+HI

#Compute the normal ordered form of the Hamiltonian
#sympy.physics.secondquant.wicks(e, **kw_args)[source]
#Returns the normal ordered equivalent of an expression using Wicks Theorem
H_N  = evaluate_deltas(wicks(H))
H_N  = substitute_dummies(H_N, new_indices=True, pretty_indices=pretty_dummies_dict)

Eref = evaluate_deltas(wicks(H, keep_only_fully_contracted=True))
Eref = substitute_dummies(Eref, new_indices=True, pretty_indices=pretty_dummies_dict)

print("Eref: ",latex(Eref))
print()
print("Normal ordered Hamiltonian") 
print(latex(H_N))
print()

#Setup Hamiltonian on normal ordered form

E0 = symbols('Eref', real = True, constant = True) #Reference energy

f  = AntiSymmetricTensor('f', (p,), (q,))
pq = NO(ap_dagger*aq)

V    = AntiSymmetricTensor('V', (p, q), (r, s))
pqsr = NO(Fd(p)*Fd(q)*F(s)*F(r)) 

HI = Rational(1, 4)*V*pqsr
Fock = f*pq #F is reserved by sympy

HN = E0+Fock+HI

#Compute <c|F|Phi_i^a> 

#Define indices above and below Fermi level
i, j, k, l = symbols('i,j,k,l', below_fermi=True)
a, b, c, d = symbols('a,b,c,d', above_fermi=True)

#Compute <c|HN|c>
cHc = wicks(HN, simplify_kronecker_deltas=True, keep_only_fully_contracted=True)
print("<c|Hnormal|c> = ", latex(cHc))

#Compute <c|F|Phi_i^a>
cFphi_ia = wicks((Fock)*NO(Fd(a)*F(i)), simplify_kronecker_deltas=True, keep_only_fully_contracted=True)
cFphi_ia = substitute_dummies(cFphi_ia)

print("<c|F|Phi_i^a> = ", latex(cFphi_ia))

