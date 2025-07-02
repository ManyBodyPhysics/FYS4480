import numpy as np

# 1. System parameters (adjustable for flexibility)
N = 6         # Number of orbitals (basis functions)
N_e = 3       # Number of electrons (occupied spin-orbitals, <= N)
dt = 0.01     # Time step for propagation
n_steps = 100 # Number of time steps to propagate

# 2. Generate random one-electron (h) and two-electron (g) integrals
np.random.seed(42)  # for reproducibility (optional)
# Generate a random symmetric matrix for one-electron integrals
h = np.random.rand(N, N) - 0.5
h = 0.5 * (h + h.T)  # make it Hermitian (symmetric real matrix)

# Generate a random 4D array for two-electron integrals and enforce symmetry
g = np.random.rand(N, N, N, N) - 0.5
# Impose physicochemical symmetries: symmetric in (p<->q), (r<->s), and (pq<->rs)
g = 0.5 * (g + g.transpose(1, 0, 2, 3))   # swap p and q
g = 0.5 * (g + g.transpose(0, 1, 3, 2))   # swap r and s
g = 0.5 * (g + g.transpose(2, 3, 0, 1))   # swap pair (pr) with (qs)
# Now g[p,q,r,s] = g[q,p,r,s] = g[p,q,s,r] = g[r,s,p,q], ensuring Hermitian properties in Fock build.

# 3. Initialize density matrix from Hartree-Fock-like ground state (core Hamiltonian guess)
# Diagonalize the one-electron Hamiltonian h to get molecular orbital energies and orbitals
eigvals, eigvecs = np.linalg.eigh(h)            # eigen-decomposition of h (since h is symmetric)
idx = np.argsort(eigvals)                       # sort eigenvalues (ascending)
eigvecs = eigvecs[:, idx]                       # sort eigenvectors correspondingly
# Build density matrix by occupying the lowest N_e eigenstates
P = np.zeros((N, N), dtype=float)
for i in range(N_e):
    v = eigvecs[:, i]
    P += np.outer(v, v)       # add |v><v| for each occupied orbital
# P is now an idempotent projector of rank N_e (approximately the HF ground-state density guess).

# 4. Define a function to build the Fock matrix F(P) given density P
def build_fock(P):
    """Construct Fock matrix F = h + J - K from density matrix P."""
    # Coulomb term: J_{pq} = sum_{r,s} P_{rs} * g_{p,q,r,s}
    J = np.einsum('pqrs,rs->pq', g, P, optimize=True)
    # Exchange term: K_{pq} = sum_{r,s} P_{rs} * g_{p,r,s,q}
    K = np.einsum('prsq,rs->pq', g, P, optimize=True)
    F = h + J - K
    return F

# 5. Function to compute total energy for monitoring (E = 0.5 * Tr[P*(h+F)])
def compute_energy(P, F):
    # Use trace formula: E = 0.5 * trace(P (h + F))
    return 0.5 * np.trace(P.dot(h + F)).real

# 6. Time-propagation using 4th-order Runge-Kutta (RK4)
# Convert P to complex type for time evolution (to allow complex off-diagonals)
P = P.astype(complex)
energies = []  # to record total energy at each time step
times = []     # to record time points

for step in range(n_steps):
    t = step * dt
    # Compute Fock matrix and energy at current time
    F = build_fock(P)
    E = compute_energy(P, F)
    energies.append(E)
    times.append(t)
    # Print or log the energy (optional, can be commented out to reduce output)
    print(f"Step {step:3d}, Time {t:.3f}: Total Energy = {E:.6f}")
    # Evaluate RK4 derivatives
    # k1
    dP1 = -1j * (F.dot(P) - P.dot(F))
    # k2: evaluate at P + 0.5*dt*k1
    P2 = P + 0.5 * dt * dP1
    F2 = build_fock(P2)
    dP2 = -1j * (F2.dot(P2) - P2.dot(F2))
    # k3: evaluate at P + 0.5*dt*k2
    P3 = P + 0.5 * dt * dP2
    F3 = build_fock(P3)
    dP3 = -1j * (F3.dot(P3) - P3.dot(F3))
    # k4: evaluate at P + dt*k3
    P4 = P + dt * dP3
    F4 = build_fock(P4)
    dP4 = -1j * (F4.dot(P4) - P4.dot(F4))
    # RK4 update: combine derivatives to update P
    P = P + (dt/6.0) * (dP1 + 2*dP2 + 2*dP3 + dP4)

# After propagation, print final results
F_final = build_fock(P)
E_final = compute_energy(P, F_final)
print(f"Final Time {n_steps*dt:.3f}: Total Energy = {E_final:.6f}")
# (Energy should remain nearly constant over time if the integration is accurate and no external perturbation.)
