import numpy as np
from numpy.polynomial.hermite import hermgauss

def compute_yukawa_integrals(n, b_x, b_y, k, M=40):
    """
    Compute the 4D tensor of two-particle Yukawa integrals I_{ij,kl} 
    for 2D harmonic oscillator basis functions up to quantum number n in x and y.
    
    Parameters:
        n   (int): Maximum HO quantum number in each dimension (0..n in x and y).
        b_x (float): HO length (oscillator width) in the x-direction.
        b_y (float): HO length in the y-direction.
        k   (float): Yukawa screening parameter.
        M   (int): Number of Gauss-Hermite quadrature points (per dimension).
        
    Returns:
        integrals (ndarray): A 4D array of shape [N, N, N, N] where N=(n+1)^2 
                             (number of 2D basis functions). The element [i,j,k,l] 
                             corresponds to basis functions i,j,k,l in lexicographic order.
    """
    # Number of basis functions in 2D
    n_basis_1d = n + 1
    N_basis = n_basis_1d ** 2  # total number of 2D basis states
    
    # Prepare an index mapping from a single index to (i_x, i_y)
    # Here we choose lexicographic: index = i_x*(n+1) + i_y
    idx_to_pair = [(ix, iy) for ix in range(n_basis_1d) for iy in range(n_basis_1d)]
    
    # Gauss-Hermite nodes and weights for integration (for weight function exp(-u^2))
    nodes, weights = hermgauss(M)
    # Convert to numpy arrays for broadcasting
    nodes = np.array(nodes)
    weights = np.array(weights)
    
    # Precompute Hermite polynomial values H_m(x) for m=0..n at all node points
    # Using the recursion: H_0 = 1, H_1 = 2x, H_{m+1} = 2x H_m - 2m H_{m-1}
    H_vals = np.zeros((n_basis_1d, M))
    H_vals[0, :] = 1.0
    if n >= 1:
        H_vals[1, :] = 2 * nodes
    for m in range(1, n):
        H_vals[m+1, :] = 2 * nodes * H_vals[m, :] - 2 * m * H_vals[m-1, :]
    
    # Precompute normalization constants for each 1D HO state (for x and y)
    # N_n = 1 / sqrt(2^n n! sqrt(pi) * b) 
    # (We use math.factorial or scipy.special.factorial for n!, but n is small here)
    from math import factorial, sqrt, pi
    norm_x = [1.0/np.sqrt((np.sqrt(pi) * (2**m) * factorial(m)) * 1.0) / np.sqrt(b_x) 
              for m in range(n_basis_1d)]
    norm_y = [1.0/np.sqrt((np.sqrt(pi) * (2**m) * factorial(m)) * 1.0) / np.sqrt(b_y) 
              for m in range(n_basis_1d)]
    
    # Initialize the result tensor
    integrals = np.zeros((N_basis, N_basis, N_basis, N_basis))
    
    # Iterate over all basis function combinations
    # We will map each basis index to (i_x,i_y) using idx_to_pair.
    # For clarity, use four nested loops over the combined index. (This is O(N_basis^4) loops,
    # which can be heavy for large n, but we assume n is small. One can use symmetries to reduce.)
    for idx_i, (i_x, i_y) in enumerate(idx_to_pair):
        # Use symmetry to skip some computations if desired (not done here for clarity)
        for idx_j, (j_x, j_y) in enumerate(idx_to_pair):
            for idx_k, (k_x, k_y) in enumerate(idx_to_pair):
                for idx_l, (l_x, l_y) in enumerate(idx_to_pair):
                    # Compute the integral (i,j | k,l)
                    
                    # We perform a double sum over Gauss-Hermite nodes for (x1,x2) and (y1,y2).
                    # To optimize, do the inner two sums (over b and d indices) vectorized.
                    total_val = 0.0
                    # Loop over 'a' and 'c' (nodes for x1 and y1 integrals)
                    for a in range(M):
                        # Factors from x1 and y1 integrals for this node 'a' and 'c'
                        H_ia = H_vals[i_x, a] * H_vals[k_x, a]  # H_{i_x}(u_a)*H_{k_x}(u_a)
                        if H_ia == 0:  # skip if zero (saves work)
                            continue
                        # x1 weight factor:
                        W_a = weights[a] * H_ia
                        u_a = nodes[a]
                        for c in range(M):
                            H_ic = H_vals[i_y, c] * H_vals[k_y, c]  # H_{i_y}(v_c)*H_{k_y}(v_c)
                            if H_ic == 0:
                                continue
                            # y1 weight factor:
                            W_c = weights[c] * H_ic
                            # Combine x1,y1 factors:
                            outer_weight = W_a * W_c
                            
                            # Prepare vectorized inner sum over b (x2) and d (y2) nodes:
                            u_diff = u_a - nodes       # vector of length M: (x_a - x_b) for all b
                            v_diff = nodes[c] - nodes  # (v_c - v_d) for all d (here nodes used for y as well)
                            
                            # Compute distance matrix R for combinations of b and d:
                            # R[b,d] = sqrt( b_x^2 * (u_a - u_b)^2  +  b_y^2 * (v_c - v_d)^2 ).
                            # We can compute this efficiently using broadcasting.
                            X_diff_sq = (b_x * u_diff)**2                        # shape (M,)
                            Y_diff_sq = (b_y * v_diff)**2                        # shape (M,)
                            # Use broadcasting to get an MxM matrix of R^2:
                            R_sq_matrix = X_diff_sq[:, None] + Y_diff_sq[None, :]
                            R_matrix = np.sqrt(R_sq_matrix)  # elementwise sqrt
                            
                            # Compute Yukawa potential values for each (b,d) pair:
                            # Avoid division by zero by setting those entries manually.
                            with np.errstate(divide='ignore', invalid='ignore'):
                                pot_matrix = np.exp(-k * R_matrix) / R_matrix   # MxM matrix
                            # Remove singular point if present (when b=a and d=c -> R=0):
                            if b_x == 0 and b_y == 0:
                                # trivial case (won't happen here as b_x,b_y > 0)
                                pass
                            else:
                                pot_matrix[a, c] = 0.0
                            
                            # Now incorporate the x2, y2 basis factors and weights:
                            # Compute vectors for b and d:
                            H_jb = H_vals[j_x, :] * H_vals[l_x, :]    # length M: H_{j_x}(u_b)*H_{l_x}(u_b) for each b
                            H_jd = H_vals[j_y, :] * H_vals[l_y, :]    # length M: H_{j_y}(v_d)*H_{l_y}(v_d) for each d
                            # Weight vectors include Gauss-Hermite weights:
                            W_b_vec = weights * H_jb   # length M
                            W_d_vec = weights * H_jd   # length M
                            
                            # Perform the double sum over b,d as a matrix contraction:
                            # inner_sum = sum_{b,d} W_b[b] * pot_matrix[b,d] * W_d[d]
                            inner_sum = W_b_vec @ (pot_matrix @ W_d_vec)
                            
                            total_val += outer_weight * inner_sum
                    # Multiply by normalization constants for i,j,k,l and Jacobian
                    total_val *= (norm_x[i_x]*norm_y[i_y] * norm_x[j_x]*norm_y[j_y] *
                                  norm_x[k_x]*norm_y[k_y] * norm_x[l_x]*norm_y[l_y] *
                                  (b_x**2 * b_y**2))
                    integrals[idx_i, idx_j, idx_k, idx_l] = total_val
    return integrals

# Example usage:
integrals = compute_yukawa_integrals(n=1, b_x=1.0, b_y=1.0, k=0.5, M=40)
print("Shape of integrals tensor:", integrals.shape)
# The element [0,0,0,0] corresponds to (i_x=i_y=0, j_x=j_y=0, k_x=k_y=0, l_x=l_y=0).
print("Integral (0,0;0,0 | 0,0;0,0) =", integrals[0, 0, 0, 0])

