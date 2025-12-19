import numpy as np
from numpy.polynomial.hermite import hermgauss

# User-specified parameters
n_max = 5          # maximum basis function index (e.g., 5 means 0 through 5)
k     = 1.0        # Yukawa screening parameter
b     = 1.0        # oscillator length (can be adjusted)

# Choose number of Gauss-Hermite points for integration
# We use an even number to avoid sampling x=0 directly.
M = max(2*(n_max + 3), 20)   # e.g., for safety, at least 20 points or 2*(n_max+3)
if M % 2 == 1:
    M += 1  # ensure even
# Get Gauss-Hermite quadrature nodes and weights
nodes, weights = hermgauss(M)  # nodes, weights for ∫ e^{-x^2} f(x) dx

# Scale nodes for our weight e^{-y^2}: (Already suited, nodes correspond to Hermite poly roots)
x_points = nodes  # rename for clarity
w_points = weights

# Evaluate Hermite polynomials H_n(x) for n=0..n_max at all nodes
# We can use a recurrence: H_0 = 1, H_1(x) = 2x, H_{n+1}(x) = 2x H_n(x) - 2n H_{n-1}(x)
Hvals = np.zeros((n_max+1, M))
Hvals[0, :] = 1.0
if n_max >= 1:
    Hvals[1, :] = 2 * x_points
for n in range(2, n_max+1):
    Hvals[n, :] = 2 * x_points * Hvals[n-1, :] - 2*(n-1) * Hvals[n-2, :]

# Precompute normalization factors N_n for each basis function
N = np.zeros(n_max+1)
for n in range(n_max+1):
    N[n] = 1.0 / np.sqrt((2**n) * np.math.factorial(n) * np.sqrt(np.pi) * b)

# Initialize tensor for two-particle integrals
size = n_max + 1
V = np.zeros((size, size, size, size))

# Compute the two-body integral tensor
for i in range(size):
    for j in range(size):
        # We will exploit symmetry in k,l vs i,j to reduce computations
        for k_idx in range(size):
            for l_idx in range(size):
                # Compute double sum for this combination
                total_sum = 0.0
                for p in range(M):
                    # Values that depend on p (x1) only:
                    H_ip = Hvals[i, p]  # H_i(x_p)
                    H_kp = Hvals[k_idx, p]  # H_k(x_p)
                    # Combine those once for efficiency
                    poly_factor_p = H_ip * H_kp
                    # Pre-factor including weight for p:
                    wp = w_points[p]
                    x_p = x_points[p]
                    for q in range(M):
                        # Compute polynomial part for q
                        H_jq = Hvals[j, q]      # H_j(x_q)
                        H_lq = Hvals[l_idx, q]  # H_l(x_q)
                        poly_factor_q = H_jq * H_lq
                        # Potential part:
                        x_q = x_points[q]
                        # Avoid division by zero (though with even M, x_p == x_q == 0 shouldn't occur)
                        dx = x_p - x_q
                        if dx == 0.0:
                            # Skip or continue (the contribution at exactly dx=0 is undefined, 
                            # but it's a zero-measure point; we handle it by skipping)
                            continue
                        V_pq = np.exp(-k * b * abs(dx)) / abs(dx)
                        # Accumulate contribution
                        total_sum += wp * w_points[q] * poly_factor_p * poly_factor_q * V_pq
                # Multiply by overall constants: b factor and normalization constants
                V[i, j, k_idx, l_idx] = b * N[i]*N[j]*N[k_idx]*N[l_idx] * total_sum

# V now contains the two-particle Yukawa integrals tensor.
# (Indices correspond to basis functions 0..n_max for i,j,k,l respectively.)


