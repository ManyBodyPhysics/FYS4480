import numpy as np
from numpy.linalg import norm
from scipy.linalg import expm  # matrix exponential
A = np.array([[0, 0.1],[0, 0   ]])
B = np.array([[0,   0 ],[0.1, 0 ]])

#Compute exponentials:
U = expm(A) @ expm(B)          # e^A e^B
U_direct = expm(A + B)         # e^{A+B}
U_BCH = expm(A + B + 0.5*(A@B - B@A))  # e^{A+B+0.5[A,B]}
print('(e^A e^B - e^{A+B}=)', norm(U - U_direct))
print('(e^A e^B - e^{A+B+0.5[A,B]}=)', norm(U - U_BCH))
