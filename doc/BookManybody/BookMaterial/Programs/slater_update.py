"""
Efficient evaluation and updating of Slater determinants.

Companion code to chapter 2 of *Quantum mechanics for Many-particle Systems*.

In a variational or diffusion Monte Carlo calculation the Slater determinant
is not evaluated once but millions of times, and together with its gradient
and Laplacian.  Recomputing it from scratch at every step costs O(d N^4) per
sweep.  The algorithms implemented here reduce that to O(d N^2):

    SlaterMatrix    -- the matrix d_ij = phi_j(r_i), its determinant and
                       inverse from an LU factorisation
    SlaterUpdater   -- the ratio R when one particle moves, the
                       Sherman-Morrison update of the inverse, and the
                       gradient and Laplacian ratios
    SpinFactorized  -- the factorisation |D| = |D_up| |D_down| and the
                       operation counts that follow from it

The single-particle basis is the two-dimensional harmonic oscillator, the
standard setting for quantum-dot Monte Carlo calculations.

Author: Morten Hjorth-Jensen
"""

import numpy as np


# ---------------------------------------------------------------------------
class HarmonicOscillator2D:
    """Two-dimensional harmonic oscillator orbitals and their derivatives.

    phi_{nx,ny}(x,y) = H_{nx}(sqrt(w) x) H_{ny}(sqrt(w) y)
                       exp(-w (x^2 + y^2) / 2)

    The gradient and the Laplacian are returned in closed form.  As the text
    stresses, one should always differentiate the orbitals analytically and
    write separate functions for the derivatives; finite differences here
    would cost accuracy and speed at once.
    """

    def __init__(self, n_orb, omega=1.0):
        self.omega = omega
        self.quantum_numbers = self._shells(n_orb)
        self.n_orb = len(self.quantum_numbers)

    @staticmethod
    def _shells(n_orb):
        """The first n_orb (nx, ny) pairs, ordered by shell nx + ny."""
        pairs = []
        shell = 0
        while len(pairs) < n_orb:
            for nx in range(shell + 1):
                pairs.append((nx, shell - nx))
                if len(pairs) == n_orb:
                    break
            shell += 1
        return pairs

    # -- Hermite polynomials and their first two derivatives ---------------
    @staticmethod
    def _hermite(n, z):
        H_prev, H = np.zeros_like(z), np.ones_like(z)
        for k in range(n):
            H_prev, H = H, 2.0 * z * H - 2.0 * k * H_prev
        return H

    def _hermite_d(self, n, z):
        return 2.0 * n * self._hermite(n - 1, z) if n > 0 else np.zeros_like(z)

    def _hermite_dd(self, n, z):
        return (4.0 * n * (n - 1) * self._hermite(n - 2, z) if n > 1
                else np.zeros_like(z))

    # ------------------------------------------------------------------
    def value(self, j, r):
        """phi_j evaluated at the position r = (x, y)."""
        nx, ny = self.quantum_numbers[j]
        s = np.sqrt(self.omega)
        x, y = np.atleast_1d(r[..., 0]), np.atleast_1d(r[..., 1])
        return (self._hermite(nx, s * x) * self._hermite(ny, s * y)
                * np.exp(-0.5 * self.omega * (x**2 + y**2)))

    def gradient(self, j, r):
        """The two components of grad phi_j at r."""
        nx, ny = self.quantum_numbers[j]
        s = np.sqrt(self.omega)
        x, y = np.atleast_1d(r[..., 0]), np.atleast_1d(r[..., 1])
        Hx, Hy = self._hermite(nx, s * x), self._hermite(ny, s * y)
        dHx = s * self._hermite_d(nx, s * x)
        dHy = s * self._hermite_d(ny, s * y)
        g = np.exp(-0.5 * self.omega * (x**2 + y**2))
        return np.stack([(dHx - self.omega * x * Hx) * Hy * g,
                         (dHy - self.omega * y * Hy) * Hx * g], axis=-1)

    def laplacian(self, j, r):
        """The Laplacian of phi_j at r."""
        nx, ny = self.quantum_numbers[j]
        w, s = self.omega, np.sqrt(self.omega)
        x, y = np.atleast_1d(r[..., 0]), np.atleast_1d(r[..., 1])
        Hx, Hy = self._hermite(nx, s * x), self._hermite(ny, s * y)
        dHx, dHy = s * self._hermite_d(nx, s*x), s * self._hermite_d(ny, s*y)
        ddHx = w * self._hermite_dd(nx, s * x)
        ddHy = w * self._hermite_dd(ny, s * y)
        g = np.exp(-0.5 * w * (x**2 + y**2))
        term_x = (ddHx - 2.0 * w * x * dHx + (w**2 * x**2 - w) * Hx) * Hy
        term_y = (ddHy - 2.0 * w * y * dHy + (w**2 * y**2 - w) * Hy) * Hx
        return (term_x + term_y) * g


# ---------------------------------------------------------------------------
class SlaterMatrix:
    """The Slater matrix d_ij = phi_j(r_i), with its determinant and inverse.

    Rows are labelled by particles and columns by the single-particle states,
    following the convention of the text.  The determinant and the inverse are
    obtained once, from an LU factorisation: the determinant as the product of
    the diagonal elements of U (Section 1.15), the inverse as N back
    substitutions.  Both cost O(N^3), and both are done once per simulation
    rather than once per move.
    """

    def __init__(self, basis, positions):
        self.basis = basis
        self.positions = np.array(positions, dtype=float)
        self.N = len(self.positions)
        self.D = self._build(self.positions)
        self.Dinv = np.linalg.inv(self.D)

    def _build(self, positions):
        N = len(positions)
        D = np.zeros((N, N))
        for i, r in enumerate(positions):
            for j in range(N):
                D[i, j] = self.basis.value(j, r)[0]
        return D

    # ------------------------------------------------------------------
    def determinant(self):
        """det(D) from the LU factorisation."""
        return np.linalg.det(self.D)

    def log_determinant(self):
        """(sign, log|det|): the only safe form for large N.

        The determinant of an N x N Slater matrix underflows or overflows
        long before N becomes interesting, so production codes carry the
        logarithm.  It comes straight out of the LU factors as the sum of
        log|u_ii|.
        """
        return np.linalg.slogdet(self.D)

    def determinant_from_svd(self):
        """|det(D)| = prod_i sigma_i, from the singular values.

        Slower than LU, but it returns the condition number along the way,
        which is what one wants near the nodal surface where det(D) -> 0.
        """
        sigma = np.linalg.svd(self.D, compute_uv=False)
        return np.prod(sigma), sigma[0] / sigma[-1]

    def condition_number(self):
        sigma = np.linalg.svd(self.D, compute_uv=False)
        return sigma[0] / sigma[-1]


# ---------------------------------------------------------------------------
class SlaterUpdater:
    """Ratios and inverse updates when a single particle moves.

    Moving particle i changes only row i of the Slater matrix.  The ratio of
    the new determinant to the old one is then the single dot product

        R = sum_j phi_j(r_i^new) d^{-1}_{ji}(r^old),

    an O(N) operation, and the inverse is repaired by the Sherman-Morrison
    formula in O(N^2).  Neither step ever recomputes a determinant.
    """

    def __init__(self, slater):
        self.s = slater
        self.refreshes = 0

    # ------------------------------------------------------------------
    def ratio(self, i, r_new):
        """R = |D(r_new)| / |D(r_old)|, in O(N) operations."""
        row = np.array([self.s.basis.value(j, r_new)[0]
                        for j in range(self.s.N)])
        return float(row @ self.s.Dinv[:, i]), row

    def accept(self, i, r_new, R=None, row=None):
        """Update D and D^{-1} after the move of particle i is accepted."""
        if row is None:
            R, row = self.ratio(i, r_new)
        Dinv = self.s.Dinv

        # S_j = (D(new) D^{-1}(old))_{ij} for every column j
        S = row @ Dinv                       # length N, S[i] == R
        new_inv = Dinv.copy()
        for j in range(self.s.N):
            if j != i:
                new_inv[:, j] = Dinv[:, j] - (S[j] / R) * Dinv[:, i]
        new_inv[:, i] = Dinv[:, i] / R

        self.s.Dinv = new_inv
        self.s.D[i, :] = row
        self.s.positions[i] = r_new

    def refresh(self):
        """Recompute the inverse from scratch by LU.

        The Sherman-Morrison updates accumulate rounding error, and near the
        nodal surface, where the Slater matrix becomes ill conditioned, they
        do so quickly.  Every production code therefore refreshes the inverse
        periodically; this is the O(N^3) operation that the updating scheme
        was designed to avoid, so it must stay rare.
        """
        self.s.Dinv = np.linalg.inv(self.s.D)
        self.refreshes += 1

    def inverse_error(self):
        """|D D^{-1} - I|, the drift accumulated by the updates."""
        return np.linalg.norm(self.s.D @ self.s.Dinv - np.eye(self.s.N))

    # ------------------------------------------------------------------
    def gradient_ratio(self, i):
        """grad_i |D| / |D| = sum_j grad_i phi_j(r_i) d^{-1}_{ji}."""
        r = self.s.positions[i]
        grads = np.array([self.s.basis.gradient(j, r)[0]
                          for j in range(self.s.N)])
        return grads.T @ self.s.Dinv[:, i]

    def laplacian_ratio(self, i):
        """lap_i |D| / |D| = sum_j lap_i phi_j(r_i) d^{-1}_{ji}."""
        r = self.s.positions[i]
        laps = np.array([self.s.basis.laplacian(j, r)[0]
                         for j in range(self.s.N)])
        return float(laps @ self.s.Dinv[:, i])

    def quantum_force(self):
        """F_i = 2 grad_i |D| / |D|, the drift term of the Langevin move."""
        return 2.0 * np.array([self.gradient_ratio(i)
                               for i in range(self.s.N)])


# ---------------------------------------------------------------------------
class SpinFactorized:
    """The factorisation |D| = |D_up| |D_down|.

    For a spin-independent Hamiltonian the full Slater determinant may be
    replaced by the product of one determinant for the spin-up particles and
    one for the spin-down particles without changing the energy.  Each factor
    is half the size, and a move of one particle touches only one of them.
    """

    def __init__(self, basis, positions_up, positions_down):
        self.up = SlaterMatrix(basis, positions_up)
        self.down = SlaterMatrix(basis, positions_down)

    def log_determinant(self):
        s1, l1 = self.up.log_determinant()
        s2, l2 = self.down.log_determinant()
        return s1 * s2, l1 + l2

    @staticmethod
    def operation_counts(N):
        """FLOP counts for one full sweep over N particles.

        Without factorisation, one move costs O(N) for the ratio and O(N^2)
        for the inverse update, so a sweep costs O(N^2) and O(N^3).  With two
        determinants of size N/2 the same sweep costs O(N^2/2) and O(N^3/4).
        """
        return {"full":       (N * N, N**3),
                "factorised": (N * N // 2, N**3 // 4)}


# ---------------------------------------------------------------------------
def _demo():
    rng = np.random.default_rng(2026)
    N = 6
    basis = HarmonicOscillator2D(n_orb=N, omega=1.0)
    positions = rng.normal(scale=1.0, size=(N, 2))

    print("=" * 74)
    print("1. The Slater matrix, its determinant and its inverse")
    print("=" * 74)
    slater = SlaterMatrix(basis, positions)
    print(f"N = {N} particles in two dimensions, quantum numbers "
          f"{basis.quantum_numbers}")
    print(f"det(D) from LU            = {slater.determinant():+.12e}")
    sign, logdet = slater.log_determinant()
    print(f"sign, log|det| from LU    = {sign:+.0f}, {logdet:.12f}")
    prod_sigma, kappa = slater.determinant_from_svd()
    print(f"prod of singular values   = {prod_sigma:.12e}")
    print(f"|det| from LU             = {abs(slater.determinant()):.12e}")
    print(f"condition number          = {kappa:.3e}")

    print()
    print("=" * 74)
    print("2. The ratio R: Sherman-Morrison against brute force")
    print("=" * 74)
    updater = SlaterUpdater(slater)
    print(f"{'move':>6s} {'R (one dot product)':>22s} "
          f"{'R (two determinants)':>22s} {'difference':>13s}")
    for k in range(4):
        i = rng.integers(N)
        r_new = slater.positions[i] + 0.3 * rng.normal(size=2)
        R, row = updater.ratio(i, r_new)

        trial = slater.positions.copy()
        trial[i] = r_new
        R_brute = (np.linalg.det(SlaterMatrix(basis, trial).D)
                   / np.linalg.det(slater.D))
        print(f"{k:6d} {R:22.14f} {R_brute:22.14f} {abs(R-R_brute):13.2e}")
        updater.accept(i, r_new, R, row)

    print()
    print("=" * 74)
    print("3. Gradient and Laplacian ratios against finite differences")
    print("=" * 74)

    def logdet_at(pos):
        return np.linalg.slogdet(SlaterMatrix(basis, pos).D)[1]

    i, eps = 2, 1.0e-5
    grad = updater.gradient_ratio(i)
    lap = updater.laplacian_ratio(i)

    fd_grad, fd_lap = np.zeros(2), 0.0
    base = logdet_at(slater.positions)
    for d in range(2):
        plus, minus = slater.positions.copy(), slater.positions.copy()
        plus[i, d] += eps
        minus[i, d] -= eps
        fd_grad[d] = (logdet_at(plus) - logdet_at(minus)) / (2 * eps)
        fd_lap += (logdet_at(plus) - 2 * base + logdet_at(minus)) / eps**2

    # the analytic Laplacian ratio is lap|D|/|D|; the finite difference of
    # log|D| gives lap|D|/|D| - |grad|D|/|D||^2, so add the square back
    fd_lap += fd_grad @ fd_grad

    print(f"grad|D|/|D|   analytic = [{grad[0]:+.10f}, {grad[1]:+.10f}]")
    print(f"              numerical= [{fd_grad[0]:+.10f}, {fd_grad[1]:+.10f}]")
    print(f"              |diff|   = {np.linalg.norm(grad-fd_grad):.2e}")
    print(f"lap|D|/|D|    analytic = {lap:+.8f}")
    print(f"              numerical= {fd_lap:+.8f}")
    print(f"              |diff|   = {abs(lap-fd_lap):.2e}")

    print()
    print("=" * 74)
    print("4. Drift of the updated inverse, and the LU refresh")
    print("=" * 74)
    print("(a) Metropolis sampling of |Phi|^2.  The nodal surface is a set of")
    print("    measure zero where the sampling density vanishes, so the walk")
    print("    almost never goes near it and the updates stay accurate.")
    walk = np.random.default_rng(11)
    s = SlaterMatrix(basis, walk.normal(size=(N, 2)))
    u = SlaterUpdater(s)
    worst, smallest_R, accepted = 0.0, np.inf, 0
    for step in range(5000):
        i = int(walk.integers(N))
        r_new = s.positions[i] + 0.4 * walk.normal(size=2)
        R, row = u.ratio(i, r_new)
        if walk.random() < min(1.0, R * R):        # Metropolis on |Phi|^2
            u.accept(i, r_new, R, row)
            accepted += 1
            smallest_R = min(smallest_R, abs(R))
        worst = max(worst, u.inverse_error())
    print(f"    5000 moves, {accepted} accepted, smallest accepted |R| = "
          f"{smallest_R:.3f}")
    print(f"    worst |D D^-1 - I| = {worst:.2e}")

    print()
    print("(b) The rare near-nodal move is what costs digits.  Each pass")
    print("    close to a node multiplies the error in the inverse by")
    print("    roughly 1/R, and the damage is permanent until the inverse is")
    print("    rebuilt.")
    s = SlaterMatrix(basis, np.random.default_rng(5).normal(size=(N, 2)))
    u = SlaterUpdater(s)
    print(f"\n    {'R of the move':>16s} {'|D D^-1 - I| after':>21s}")
    print(f"    {'(start)':>16s} {u.inverse_error():21.2e}")
    for sep in (1.0e-2, 1.0e-4, 1.0e-6, 1.0e-8):
        r_node = s.positions[0] + np.array([sep, 0.0])
        R, row = u.ratio(1, r_node)
        u.accept(1, r_node, R, row)                # step onto the node
        r_back = s.positions[0] + np.array([1.0, 0.7])
        R2, row2 = u.ratio(1, r_back)
        u.accept(1, r_back, R2, row2)              # and step away again
        print(f"    {R:16.2e} {u.inverse_error():21.2e}")
    u.refresh()
    print(f"    {'after LU refresh':>16s} {u.inverse_error():21.2e}")
    print()
    print("    This is why production codes rebuild the inverse from an LU")
    print("    factorisation every few hundred moves: the O(N^3) cost is")
    print("    amortised to nothing, and the accumulated error is reset.")

    print()
    print("=" * 74)
    print("5. Approaching the nodal surface")
    print("=" * 74)
    print("As two particles approach the same point the Slater matrix becomes")
    print("singular, the determinant vanishes and the inverse blows up.")
    print(f"\n{'separation':>12s} {'|det D|':>14s} {'kappa_2(D)':>14s} "
          f"{'|D D^-1 - I|':>15s}")
    base_pos = rng.normal(size=(N, 2))
    for sep in (1.0, 1.0e-1, 1.0e-2, 1.0e-4, 1.0e-6):
        pos = base_pos.copy()
        pos[1] = pos[0] + np.array([sep, 0.0])
        s = SlaterMatrix(basis, pos)
        err = np.linalg.norm(s.D @ s.Dinv - np.eye(N))
        print(f"{sep:12.0e} {abs(np.linalg.det(s.D)):14.3e} "
              f"{s.condition_number():14.3e} {err:15.2e}")
    print("\nThe condition number of Section 1.13 is the diagnostic: it is the")
    print("ratio of the largest to the smallest singular value, and it tells")
    print("us how many digits the inverse has lost.")

    print()
    print("=" * 74)
    print("6. Spin factorisation and the operation counts")
    print("=" * 74)
    Nup = Ndown = 3
    b = HarmonicOscillator2D(n_orb=max(Nup, Ndown))
    sf = SpinFactorized(b, rng.normal(size=(Nup, 2)),
                        rng.normal(size=(Ndown, 2)))
    sign, logdet = sf.log_determinant()
    print(f"log|D_up| + log|D_down| = {logdet:.10f}   (sign {sign:+.0f})")
    print()
    print(f"{'N':>5s} {'ratios, full':>14s} {'inverse, full':>15s}"
          f" {'ratios, split':>15s} {'inverse, split':>16s}")
    for n in (10, 50, 100):
        counts = SpinFactorized.operation_counts(n)
        rf, invf = counts["full"]
        rs, invs = counts["factorised"]
        print(f"{n:5d} {rf:14.1e} {invf:15.1e} {rs:15.1e} {invs:16.1e}")


if __name__ == "__main__":
    _demo()
