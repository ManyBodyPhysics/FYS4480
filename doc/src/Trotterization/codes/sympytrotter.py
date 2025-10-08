from sympy.physics.quantum import Commutator, Operator
from sympy import Rational, expand
X, Y = Operator('X'), Operator('Y')
# BCH series up to third order:
Z =X+Y+Rational(1,2)*Commutator(X, Y)+ Rational(1,12)*(Commutator(X, Commutator(X,Y))+ Commutator(Y, Commutator(Y,X)))
print(Z.expand(commutator=True))   
