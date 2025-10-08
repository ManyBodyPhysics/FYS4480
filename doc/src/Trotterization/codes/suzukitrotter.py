import numpy as np
from numpy.linalg import norm
from scipy.linalg import expm

#Define Pauli matrices
X = np.array([[0, 1],[1, 0]])
Z = np.array([[1, 0],[0,-1]])
H = X + Z
t = 1.0
N = 16
dt = t/N
#First-order Trotter approximation
U_trot = np.eye(2)
for k in range(N):
    U_trot = expm(-1j * X * dt) @ expm(-1j * Z * dt) @ U_trot
#Exact evolution
U_exact = expm(-1j * H * t)
error = norm(U_trot - U_exact)
print(error)




 









