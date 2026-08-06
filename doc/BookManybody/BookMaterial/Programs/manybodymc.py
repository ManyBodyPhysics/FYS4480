"""
Variational and diffusion Monte Carlo for N electrons in a two-dimensional dot.

Companion code to chapter 16 of *Quantum Mechanics for Many-particle Systems*.

Chapters 13 to 15 solved the two-electron dot.  Two electrons in a singlet need
no determinant: the spatial wave function is symmetric, nodeless and positive,
and a Gaussian times a Jastrow factor is the whole ansatz.  Six electrons need
all of it.  The trial function is a Slater-Jastrow product

    Psi_T(R, sigma) = |D_up| |D_down| exp( sum_{i<j} a_ij r_ij / (1 + b r_ij) )

with |D_up| and |D_down| determinants of scaled oscillator orbitals over the
spin-up and spin-down electrons respectively.  Three pieces of machinery that
were unnecessary for two electrons are now essential, and all three come from
Chapter 2:

  * the ratio of two determinants after a single-particle move, obtained from
    the inverse Slater matrix in O(N) rather than O(N^3) operations
    (section 2.11, Sherman-Morrison);

  * the gradient and Laplacian ratios, which supply the quantum force and the
    kinetic part of the local energy, again from the inverse (section 2.12);

  * the factorisation |D| = |D_up| |D_down|, which halves each determinant
    (section 2.13, the Moskowitz-Kalos ansatz).

The cusp coefficients a_ij are fixed, not variational, and they are *spin
dependent*: in two dimensions a = 1 for an antiparallel pair and a = 1/3 for a
parallel one, because the determinant already supplies a node between parallel
electrons and the Jastrow factor must not supply a second one.  This is what
makes the spin assignment matter dynamically rather than only in bookkeeping.

Spin sampling.  The Moskowitz-Kalos argument says the energy does not depend on
*which* electrons are assigned to which determinant, so a production code fixes
the assignment once and never touches it.  This program instead treats the spin
configuration as a stochastic variable and proposes exchanges of an unlike
pair, which conserves S_z and which -- because a_ij is spin dependent -- has a
non-trivial acceptance rate.  That turns the Moskowitz-Kalos claim from an
assertion into something measured, and it adds a move that relabels electrons
without moving them, which no amount of diffusion can accomplish.

Author: Morten Hjorth-Jensen
"""

import math

import numpy as np

from slater_update import HarmonicOscillator2D

# ---------------------------------------------------------------------------
#  The cusp coefficients in two dimensions
# ---------------------------------------------------------------------------
#  Requiring the local energy to stay finite as r_ij -> 0 gives
#  a = 1/(2 l + d - 1), with l = 0 for an antiparallel pair (the spatial
#  function is finite at contact) and l = 1 for a parallel one (the
#  determinant already vanishes linearly).  In d = 2 that is 1 and 1/3.
# ---------------------------------------------------------------------------
A_ANTIPARALLEL = 1.0
A_PARALLEL = 1.0 / 3.0

#  Reference values from chapter 11, six oscillator shells (42 spin-orbitals),
#  hbar omega = 1.
REFERENCE = {
    2: dict(hf=3.16192140, mp2=3.02703812, ccsd=3.01362613, exact=3.0),
    6: dict(hf=20.72025707, mp2=20.30256126, ccsd=20.27324663, exact=None),
}


# ---------------------------------------------------------------------------
def cusp_matrix(spins):
    """a_ij for every pair, from the spins (+1 up, -1 down)."""
    spins = np.asarray(spins)
    parallel = spins[:, None] == spins[None, :]
    return np.where(parallel, A_PARALLEL, A_ANTIPARALLEL)


def safe_inverse(D):
    """Inverse of one matrix or a stack of them, tolerant of nodal walkers.

    A configuration sitting on a nodal surface of the determinant has a
    singular Slater matrix.  The Metropolis test rejects such a proposal
    anyway, since |Psi_T|^2 vanishes there, but the inverse is needed to
    evaluate the proposal before it can be rejected, so it must not raise.
    The pseudo-inverse is used for the offending walkers and nothing else
    changes.  This is unnecessary for two electrons -- the singlet trial
    function of Chapter 13 has no nodes at all -- and unavoidable here.
    """
    try:
        return np.linalg.inv(D)
    except np.linalg.LinAlgError:
        stack = np.atleast_3d(D) if D.ndim == 3 else D[None, :, :]
        out = np.empty_like(stack)
        for m in range(stack.shape[0]):
            try:
                out[m] = np.linalg.inv(stack[m])
            except np.linalg.LinAlgError:
                out[m] = np.linalg.pinv(stack[m])
        return out if D.ndim == 3 else out[0]


# ---------------------------------------------------------------------------
#  A single walker, with the full Chapter-2 machinery
# ---------------------------------------------------------------------------
class SlaterJastrow:
    """One configuration (positions and spins) and everything derived from it.

    The two spin blocks are held as separate Slater matrices with their own
    inverses, exactly as Eq. (2.51) prescribes.  A move of one particle touches
    one row of one block, so the ratio is a dot product and the inverse is
    repaired by Sherman-Morrison.
    """

    def __init__(self, n_particles, alpha, beta, omega=1.0,
                 positions=None, spins=None, rng=None):
        if n_particles % 2:
            raise ValueError("closed shells only: N must be even")
        self.N = n_particles
        self.alpha, self.beta, self.omega = alpha, beta, omega
        # The orbitals are oscillator functions at the *scaled* frequency
        # alpha*omega.  At alpha = 1 they are the exact single-particle
        # eigenfunctions of the trap.
        self.basis = HarmonicOscillator2D(n_orb=n_particles // 2,
                                          omega=alpha * omega)
        if rng is None:
            rng = np.random.default_rng(2024)
        if positions is None:
            positions = rng.normal(0.0, 1.0, (n_particles, 2))
        if spins is None:
            spins = np.array([1] * (n_particles // 2)
                             + [-1] * (n_particles // 2))
        self.positions = np.array(positions, dtype=float)
        self.spins = np.array(spins, dtype=int)
        self.refresh()

    # -- construction and repair -------------------------------------------
    def refresh(self):
        """Rebuild both Slater matrices and both inverses from scratch."""
        self.up = np.flatnonzero(self.spins > 0)
        self.down = np.flatnonzero(self.spins < 0)
        self.D_up, self.Dinv_up = self._block(self.up)
        self.D_down, self.Dinv_down = self._block(self.down)
        self.a = cusp_matrix(self.spins)

    def _block(self, index):
        n = len(index)
        D = np.empty((n, n))
        for row, i in enumerate(index):
            for j in range(n):
                D[row, j] = float(self.basis.value(j, self.positions[i])[0])
        return D, safe_inverse(D)

    def _locate(self, i):
        """(block index array, inverse, row of particle i in that block)."""
        if self.spins[i] > 0:
            return self.up, self.Dinv_up, int(np.flatnonzero(self.up == i)[0])
        return self.down, self.Dinv_down, int(np.flatnonzero(self.down == i)[0])

    # -- Eq. (2.44): the determinant ratio in O(N) --------------------------
    def ratio(self, i, r_new):
        """|D(r_new)| / |D(r_old)| for a move of particle i alone."""
        _, Dinv, row = self._locate(i)
        n = Dinv.shape[0]
        values = np.array([float(self.basis.value(j, r_new)[0])
                           for j in range(n)])
        return float(values @ Dinv[:, row]), values

    def accept(self, i, r_new, R, values):
        """Sherman-Morrison repair of the inverse, Eqs. (2.47) and (2.48)."""
        index, Dinv, row = self._locate(i)
        S = values @ Dinv
        updated = Dinv.copy()
        for j in range(Dinv.shape[0]):
            if j != row:
                updated[:, j] = Dinv[:, j] - (S[j] / R) * Dinv[:, row]
        updated[:, row] = Dinv[:, row] / R
        if self.spins[i] > 0:
            self.Dinv_up, self.D_up[row, :] = updated, values
        else:
            self.Dinv_down, self.D_down[row, :] = updated, values
        self.positions[i] = r_new

    # -- Eqs. (2.49) and (2.50): derivatives of the determinant -------------
    def determinant_gradient(self, i):
        """grad_i |D| / |D|, a two-vector."""
        _, Dinv, row = self._locate(i)
        n = Dinv.shape[0]
        out = np.zeros(2)
        for j in range(n):
            out += (self.basis.gradient(j, self.positions[i])[0]
                    * Dinv[j, row])
        return out

    def determinant_laplacian(self, i):
        """laplacian_i |D| / |D|, a scalar."""
        _, Dinv, row = self._locate(i)
        n = Dinv.shape[0]
        return float(sum(float(self.basis.laplacian(j, self.positions[i])[0])
                         * Dinv[j, row] for j in range(n)))

    # -- the Jastrow factor and its derivatives ----------------------------
    def _separations(self):
        diff = self.positions[:, None, :] - self.positions[None, :, :]
        dist = np.linalg.norm(diff, axis=-1)
        np.fill_diagonal(dist, np.inf)
        return diff, dist

    def jastrow_exponent(self):
        _, dist = self._separations()
        upper = np.triu_indices(self.N, 1)
        r = dist[upper]
        return float(np.sum(self.a[upper] * r / (1.0 + self.beta * r)))

    def jastrow_gradient(self, i):
        """grad_i U, with U the exponent of the Jastrow factor."""
        diff, dist = self._separations()
        d = 1.0 / (1.0 + self.beta * dist[i])
        weight = self.a[i] * d * d / dist[i]
        weight[i] = 0.0
        return np.sum(weight[:, None] * diff[i], axis=0)

    def jastrow_laplacian(self, i):
        """laplacian_i U in two dimensions."""
        _, dist = self._separations()
        r = dist[i]
        d = 1.0 / (1.0 + self.beta * r)
        terms = self.a[i] * (d * d / r - 2.0 * self.beta * d**3)
        terms[i] = 0.0
        return float(np.sum(terms))

    # -- the quantities the sampling needs ---------------------------------
    def log_psi(self):
        s_up, l_up = np.linalg.slogdet(self.D_up)
        s_dn, l_dn = np.linalg.slogdet(self.D_down)
        return l_up + l_dn + self.jastrow_exponent(), s_up * s_dn

    def quantum_force(self):
        """F_i = 2 grad_i ln Psi_T, one row per particle."""
        force = np.empty((self.N, 2))
        for i in range(self.N):
            force[i] = 2.0 * (self.determinant_gradient(i)
                              + self.jastrow_gradient(i))
        return force

    def local_energy(self):
        """E_L = H Psi_T / Psi_T.

        With Psi_T = |D| exp(U),

            lap Psi / Psi = lap|D|/|D| + 2 (grad|D|/|D|) . grad U
                            + lap U + |grad U|^2 ,

        so the determinant and the Jastrow factor never have to be
        differentiated together -- only their separate derivatives are needed,
        and each of those is an O(N) dot product with the inverse.
        """
        kinetic = 0.0
        for i in range(self.N):
            g_det = self.determinant_gradient(i)
            g_jas = self.jastrow_gradient(i)
            kinetic += (self.determinant_laplacian(i)
                        + self.jastrow_laplacian(i)
                        + float(g_jas @ g_jas)
                        + 2.0 * float(g_det @ g_jas))
        trap = 0.5 * self.omega**2 * float(np.sum(self.positions**2))
        _, dist = self._separations()
        upper = np.triu_indices(self.N, 1)
        coulomb = float(np.sum(1.0 / dist[upper]))
        return -0.5 * kinetic + trap + coulomb

    # -- derivatives with respect to the variational parameters ------------
    def parameter_gradient(self):
        """d ln Psi_T / d(alpha, beta), the O_theta of Eq. (14.2).

        The beta derivative is elementary.  The alpha derivative acts through
        the frequency of the orbitals, and is assembled as
        Tr(D^{-1} dD/dalpha) for each spin block.
        """
        d_alpha = 0.0
        for index, Dinv in ((self.up, self.Dinv_up),
                            (self.down, self.Dinv_down)):
            n = len(index)
            dD = np.empty((n, n))
            for row, i in enumerate(index):
                for j in range(n):
                    dD[row, j] = float(
                        self.basis.d_value_d_omega(j, self.positions[i])[0])
            d_alpha += float(np.trace(Dinv @ dD)) * self.omega

        _, dist = self._separations()
        upper = np.triu_indices(self.N, 1)
        r = dist[upper]
        d = 1.0 / (1.0 + self.beta * r)
        d_beta = -float(np.sum(self.a[upper] * r * r * d * d))
        return np.array([d_alpha, d_beta])


# ---------------------------------------------------------------------------
#  The derivative of an orbital with respect to the trap frequency
# ---------------------------------------------------------------------------
def _d_value_d_omega(self, j, r):
    """d phi_j / d omega, needed only for the parameter gradient.

    phi = H_{nx}(s x) H_{ny}(s y) exp(-w (x^2+y^2)/2) with s = sqrt(w), so

        dphi/dw = [ (x/2s) H'_{nx} H_{ny} + (y/2s) H_{nx} H'_{ny}
                    - (x^2+y^2)/2 H_{nx} H_{ny} ] exp(...) .
    """
    nx, ny = self.quantum_numbers[j]
    w = self.omega
    s = np.sqrt(w)
    x, y = np.atleast_1d(r[..., 0]), np.atleast_1d(r[..., 1])
    Hx, Hy = self._hermite(nx, s * x), self._hermite(ny, s * y)
    dHx, dHy = self._hermite_d(nx, s * x), self._hermite_d(ny, s * y)
    g = np.exp(-0.5 * w * (x**2 + y**2))
    return (x / (2.0 * s) * dHx * Hy + y / (2.0 * s) * Hx * dHy
            - 0.5 * (x**2 + y**2) * Hx * Hy) * g


HarmonicOscillator2D.d_value_d_omega = _d_value_d_omega


# ---------------------------------------------------------------------------
#  Vectorised evaluation over an ensemble of walkers
# ---------------------------------------------------------------------------
#  The single-walker class above is the readable statement of the algorithm and
#  is what the verification uses.  For production runs the determinants are
#  small -- 3 x 3 for six electrons -- and numpy's batched LU over an array of
#  shape (M, n, n) beats a Python loop carrying Sherman-Morrison updates by a
#  wide margin.  The crossover is measured in the demo.
# ---------------------------------------------------------------------------
class Ensemble:
    """M walkers, each with N positions and N spins."""

    def __init__(self, n_walkers, n_particles, alpha, beta, omega=1.0,
                 rng=None):
        self.M, self.N = n_walkers, n_particles
        self.alpha, self.beta, self.omega = alpha, beta, omega
        self.n_orb = n_particles // 2
        self.basis = HarmonicOscillator2D(n_orb=self.n_orb,
                                          omega=alpha * omega)
        rng = np.random.default_rng(2024) if rng is None else rng
        self.positions = rng.normal(0.0, 1.0, (n_walkers, n_particles, 2))
        base = np.array([1] * self.n_orb + [-1] * self.n_orb)
        self.spins = np.tile(base, (n_walkers, 1))

    # -- bookkeeping --------------------------------------------------------
    def blocks(self, spins=None):
        """Row indices of the spin-up and spin-down electrons, per walker."""
        spins = self.spins if spins is None else spins
        up = np.argsort(-spins, kind="stable")[:, :self.n_orb]
        down = np.argsort(spins, kind="stable")[:, :self.n_orb]
        return up, down

    def _matrices(self, positions, index):
        """D[m, p, j] = phi_j(r_{index[m,p]}) for one spin block."""
        sub = np.take_along_axis(positions, index[:, :, None], axis=1)
        D = np.empty((positions.shape[0], self.n_orb, self.n_orb))
        for j in range(self.n_orb):
            D[:, :, j] = self.basis.value(j, sub)
        return D, sub

    # -- the wave function and its derivatives ------------------------------
    def log_psi(self, positions=None, spins=None):
        positions = self.positions if positions is None else positions
        spins = self.spins if spins is None else spins
        up, down = self.blocks(spins)
        total = 0.0
        for index in (up, down):
            D, _ = self._matrices(positions, index)
            total = total + np.linalg.slogdet(D)[1]
        return total + self._jastrow(positions, spins)

    def _pairs(self, positions):
        diff = positions[:, :, None, :] - positions[:, None, :, :]
        dist = np.linalg.norm(diff, axis=-1)
        eye = np.eye(self.N, dtype=bool)
        dist[:, eye] = np.inf
        return diff, dist

    def _cusp(self, spins):
        parallel = spins[:, :, None] == spins[:, None, :]
        return np.where(parallel, A_PARALLEL, A_ANTIPARALLEL)

    def _jastrow(self, positions, spins):
        _, dist = self._pairs(positions)
        a = self._cusp(spins)
        upper = np.triu_indices(self.N, 1)
        r = dist[:, upper[0], upper[1]]
        return np.sum(a[:, upper[0], upper[1]] * r / (1.0 + self.beta * r),
                      axis=1)

    def _determinant_derivatives(self, positions, spins):
        """grad_i |D|/|D| and lap_i |D|/|D| for every particle of every walker.

        Equations (2.49) and (2.50) applied blockwise.  Only the particles in a
        given block contribute to that block's determinant, so the two spin
        sectors are assembled independently and scattered back into an (M, N)
        array indexed by particle.
        """
        M, N = positions.shape[0], self.N
        grad = np.zeros((M, N, 2))
        lap = np.zeros((M, N))
        for index in self.blocks(spins):
            D, sub = self._matrices(positions, index)
            Dinv = safe_inverse(D)
            g = np.empty((M, self.n_orb, self.n_orb, 2))
            l = np.empty((M, self.n_orb, self.n_orb))
            for j in range(self.n_orb):
                g[:, :, j, :] = self.basis.gradient(j, sub)
                l[:, :, j] = self.basis.laplacian(j, sub)
            # sum_j (d phi_j / d r_p) Dinv[j, p]
            weights = np.transpose(Dinv, (0, 2, 1))       # weights[m, p, j]
            grad_block = np.einsum("mpjd,mpj->mpd", g, weights)
            lap_block = np.einsum("mpj,mpj->mp", l, weights)
            np.put_along_axis(grad, index[:, :, None], grad_block, axis=1)
            np.put_along_axis(lap, index, lap_block, axis=1)
        return grad, lap

    def _jastrow_derivatives(self, positions, spins):
        diff, dist = self._pairs(positions)
        a = self._cusp(spins)
        d = 1.0 / (1.0 + self.beta * dist)
        weight = a * d * d / dist                          # (M, N, N)
        grad = np.einsum("mij,mijd->mid", weight, diff)
        lap = np.sum(a * (d * d / dist - 2.0 * self.beta * d**3), axis=2)
        return grad, lap

    def quantum_force(self, positions=None, spins=None):
        positions = self.positions if positions is None else positions
        spins = self.spins if spins is None else spins
        g_det, _ = self._determinant_derivatives(positions, spins)
        g_jas, _ = self._jastrow_derivatives(positions, spins)
        return 2.0 * (g_det + g_jas)

    def local_energy(self, positions=None, spins=None):
        positions = self.positions if positions is None else positions
        spins = self.spins if spins is None else spins
        g_det, l_det = self._determinant_derivatives(positions, spins)
        g_jas, l_jas = self._jastrow_derivatives(positions, spins)
        kinetic = np.sum(l_det + l_jas + np.sum(g_jas * g_jas, axis=2)
                         + 2.0 * np.sum(g_det * g_jas, axis=2), axis=1)
        trap = 0.5 * self.omega**2 * np.sum(positions**2, axis=(1, 2))
        _, dist = self._pairs(positions)
        upper = np.triu_indices(self.N, 1)
        coulomb = np.sum(1.0 / dist[:, upper[0], upper[1]], axis=1)
        return -0.5 * kinetic + trap + coulomb

    def parameter_gradient(self, positions=None, spins=None):
        """O_theta = d ln Psi_T / d(alpha, beta) for every walker, (M, 2)."""
        positions = self.positions if positions is None else positions
        spins = self.spins if spins is None else spins
        M = positions.shape[0]
        d_alpha = np.zeros(M)
        for index in self.blocks(spins):
            D, sub = self._matrices(positions, index)
            Dinv = safe_inverse(D)
            dD = np.empty_like(D)
            for j in range(self.n_orb):
                dD[:, :, j] = self.basis.d_value_d_omega(j, sub)
            d_alpha += np.einsum("mij,mji->m", Dinv, dD) * self.omega
        _, dist = self._pairs(positions)
        a = self._cusp(spins)
        upper = np.triu_indices(self.N, 1)
        r = dist[:, upper[0], upper[1]]
        d = 1.0 / (1.0 + self.beta * r)
        d_beta = -np.sum(a[:, upper[0], upper[1]] * r * r * d * d, axis=1)
        return np.stack([d_alpha, d_beta], axis=1)


# ---------------------------------------------------------------------------
#  Moves
# ---------------------------------------------------------------------------
def limit_drift(force, time_step, a=1.0):
    """Umrigar's cutoff on the drift velocity.

    The quantum force F = 2 grad ln Psi_T diverges on a nodal surface of the
    determinant, where Psi_T vanishes.  A walker that strays close to a node is
    then thrown across the configuration space by a drift step of size
    D F dt, and the local energy it reports on the way is enormous.  For two
    electrons in Chapters 13 to 15 this never arose, because that trial
    function has no nodes.  Here it must be dealt with.

    The standard remedy replaces F by

        F -> F ( -1 + sqrt(1 + 2 a F^2 dt) ) / ( a F^2 dt ) ,

    which is the identity for small F and behaves as 1/(F dt) for large F, so
    the drift step saturates instead of diverging.  The important point is that
    the *same* modified drift is used both to propose the move and to evaluate
    the Green's function, so detailed balance is untouched and the sampled
    distribution is still exactly |Psi_T|^2.  This is the lesson of
    Section 13.6, where a deliberately wrong quantum force was shown to cost
    efficiency and not correctness -- used here on purpose.
    """
    v2 = np.sum(force * force, axis=-1, keepdims=True)
    scaled = a * v2 * time_step
    factor = np.where(scaled > 1e-12,
                      (-1.0 + np.sqrt(1.0 + 2.0 * scaled))
                      / np.where(scaled > 1e-12, scaled, 1.0),
                      1.0)
    return force * factor


def _greens_ratio(old, new, f_old, f_new, time_step, diffusion=0.5):
    """ln[G(x,y)/G(y,x)] for the drifting Gaussian, Eq. (13.20)."""
    return np.sum(0.5 * (f_old + f_new)
                  * (0.5 * diffusion * time_step * (f_old - f_new)
                     - new + old), axis=-1)


def _spatial_sweep(ens, positions, spins, log_old, force_old, time_step, rng,
                   drift_cutoff=True):
    """Move every particle once, with drift and Metropolis-Hastings."""
    M, N = positions.shape[0], ens.N
    root_dt, diffusion = math.sqrt(time_step), 0.5
    accepted = 0
    attempted = accepted_distance = 0.0

    cut = (lambda f: limit_drift(f, time_step)) if drift_cutoff else (lambda f: f)
    drift_old = cut(force_old)
    for i in range(N):
        trial = positions.copy()
        step = (diffusion * drift_old[:, i] * time_step
                + rng.normal(0.0, 1.0, (M, 2)) * root_dt)
        trial[:, i] = positions[:, i] + step
        log_new = ens.log_psi(trial, spins)
        force_new = ens.quantum_force(trial, spins)
        drift_new = cut(force_new)
        green = _greens_ratio(positions[:, i], trial[:, i],
                              drift_old[:, i], drift_new[:, i], time_step)
        take = (green + 2.0 * (log_new - log_old)
                > np.log(rng.random(M) + 1e-300))

        moved = np.sum(step * step, axis=1)
        attempted += float(moved.sum())
        accepted_distance += float(moved[take].sum())
        accepted += int(take.sum())

        positions[take] = trial[take]
        log_old = np.where(take, log_new, log_old)
        force_old[take] = force_new[take]
        drift_old[take] = drift_new[take]

    return (positions, log_old, force_old, accepted / (M * N),
            attempted, accepted_distance)


def _spin_exchange(ens, positions, spins, log_old, rng):
    """Propose swapping the spins of one unlike pair, per walker.

    The proposal is symmetric -- an up electron and a down electron are drawn
    uniformly -- so the ordinary Metropolis test applies.  S_z is conserved
    exactly.  The move changes the wave function through both determinants (the
    two electrons change blocks) and through the Jastrow factor (their cusp
    coefficients with every other electron change), which is why it has a
    non-trivial acceptance rate at all.
    """
    M, n = positions.shape[0], ens.n_orb
    up, down = ens.blocks(spins)
    pick_up = up[np.arange(M), rng.integers(n, size=M)]
    pick_down = down[np.arange(M), rng.integers(n, size=M)]

    trial = spins.copy()
    trial[np.arange(M), pick_up] = -1
    trial[np.arange(M), pick_down] = +1
    log_new = ens.log_psi(positions, trial)
    take = 2.0 * (log_new - log_old) > np.log(rng.random(M) + 1e-300)

    spins[take] = trial[take]
    log_old = np.where(take, log_new, log_old)
    return spins, log_old, float(take.mean())


# ---------------------------------------------------------------------------
#  Variational Monte Carlo
# ---------------------------------------------------------------------------
def vmc(n_particles=6, alpha=1.0, beta=0.4, omega=1.0, n_walkers=400,
        n_steps=2000, time_step=0.05, burn_in=300, rng=None,
        sample_spins=True, want_gradient=False, drift_cutoff=True):
    """Importance-sampled VMC, optionally sampling the spin assignment too."""
    rng = np.random.default_rng(2024) if rng is None else rng
    ens = Ensemble(n_walkers, n_particles, alpha, beta, omega, rng)
    positions, spins = ens.positions, ens.spins
    log_old = ens.log_psi(positions, spins)
    force_old = ens.quantum_force(positions, spins)

    series = np.empty(n_steps)
    gradient_sum = np.zeros(2)
    gradient_energy_sum = np.zeros(2)
    energy_sum = 0.0
    acceptance = spin_acceptance = 0.0

    for step in range(n_steps + burn_in):
        (positions, log_old, force_old, rate,
         _, _) = _spatial_sweep(ens, positions, spins, log_old, force_old,
                                time_step, rng, drift_cutoff)
        acceptance += rate
        if sample_spins:
            spins, log_old, srate = _spin_exchange(ens, positions, spins,
                                                   log_old, rng)
            spin_acceptance += srate
            force_old = ens.quantum_force(positions, spins)
        if step >= burn_in:
            energies = ens.local_energy(positions, spins)
            series[step - burn_in] = energies.mean()
            if want_gradient:
                O = ens.parameter_gradient(positions, spins)
                gradient_sum += O.mean(axis=0)
                gradient_energy_sum += (O * energies[:, None]).mean(axis=0)
                energy_sum += energies.mean()

    steps = n_steps + burn_in
    out = dict(series=series, energy=float(series.mean()),
               acceptance=acceptance / steps,
               spin_acceptance=spin_acceptance / steps if sample_spins else 0.0,
               variance=float(ens.local_energy(positions, spins).var(ddof=1)),
               alpha=alpha, beta=beta, omega=omega, n_particles=n_particles,
               sample_spins=sample_spins)
    if want_gradient:
        mean_energy = energy_sum / n_steps
        out["gradient"] = 2.0 * (gradient_energy_sum / n_steps
                                 - gradient_sum / n_steps * mean_energy)
    return out


def optimise(n_particles=6, alpha=1.0, beta=0.4, omega=1.0, learning_rate=0.02,
             max_iter=25, tol=2e-3, n_walkers=200, n_steps=200, rng=None,
             sample_spins=True, verbose=False):
    """Gradient descent on (alpha, beta) using Eq. (14.2).

    Short, noisy runs are enough: we are locating a minimum, not measuring an
    energy.  That was the lesson of chapter 14 and it is unchanged here.
    """
    rng = np.random.default_rng(7) if rng is None else rng
    theta = np.array([alpha, beta], dtype=float)
    history = []
    for iteration in range(max_iter):
        run = vmc(n_particles, theta[0], theta[1], omega,
                  n_walkers=n_walkers, n_steps=n_steps, burn_in=100,
                  time_step=0.05, rng=rng, sample_spins=sample_spins,
                  want_gradient=True)
        gradient = run["gradient"]
        history.append((theta.copy(), run["energy"], gradient.copy()))
        if verbose:
            print(f"   {iteration:3d}  alpha {theta[0]:.4f}  "
                  f"beta {theta[1]:.4f}  E {run['energy']:12.6f}  "
                  f"|grad| {np.linalg.norm(gradient):.4f}")
        if np.linalg.norm(gradient) < tol:
            break
        theta = theta - learning_rate * gradient
        theta[0] = min(max(theta[0], 0.3), 2.0)
        theta[1] = min(max(theta[1], 0.02), 3.0)
    return theta, history


# ---------------------------------------------------------------------------
#  Diffusion Monte Carlo
# ---------------------------------------------------------------------------
def dmc(n_particles=6, alpha=1.0, beta=0.4, omega=1.0, n_walkers=400,
        n_steps=3000, time_step=0.01, burn_in=600, rng=None,
        feedback=0.05, sample_spins=True, effective_time_step=True,
        drift_cutoff=True):
    """Importance-sampled DMC, following chapter 15 with a Slater-Jastrow
    guiding function and, optionally, spin-exchange moves."""
    rng = np.random.default_rng(2024) if rng is None else rng
    ens = Ensemble(n_walkers, n_particles, alpha, beta, omega, rng)
    positions, spins = ens.positions, ens.spins

    # equilibrate in the variational distribution before branching starts
    log_old = ens.log_psi(positions, spins)
    force_old = ens.quantum_force(positions, spins)
    for _ in range(200):
        positions, log_old, force_old, _, _, _ = _spatial_sweep(
            ens, positions, spins, log_old, force_old, 0.05, rng)
        if sample_spins:
            spins, log_old, _ = _spin_exchange(ens, positions, spins,
                                               log_old, rng)
            force_old = ens.quantum_force(positions, spins)

    target = n_walkers
    energy_old = ens.local_energy(positions, spins)
    trial_energy = float(energy_old.mean())
    running = trial_energy
    series = np.empty(n_steps)
    growth = np.empty(n_steps)
    population = np.empty(n_steps + burn_in, dtype=int)
    acceptance = spin_acceptance = 0.0

    for step in range(n_steps + burn_in):
        (positions, log_old, force_old, rate, attempted,
         accepted_distance) = _spatial_sweep(ens, positions, spins, log_old,
                                             force_old, time_step, rng,
                                             drift_cutoff)
        acceptance += rate
        if sample_spins:
            spins, log_old, srate = _spin_exchange(ens, positions, spins,
                                                   log_old, rng)
            spin_acceptance += srate
            force_old = ens.quantum_force(positions, spins)
        energy_new = ens.local_energy(positions, spins)

        dt_branch = time_step
        if effective_time_step and attempted > 0.0:
            dt_branch = time_step * accepted_distance / attempted

        weight = np.exp(-dt_branch * (0.5 * (energy_old + energy_new)
                                      - trial_energy))
        copies = np.clip(np.floor(weight + rng.random(len(weight))),
                         0, 3).astype(int)
        if copies.sum() == 0:
            copies[rng.integers(len(copies))] = 1

        positions = np.repeat(positions, copies, axis=0)
        spins = np.repeat(spins, copies, axis=0)
        log_old = np.repeat(log_old, copies)
        force_old = np.repeat(force_old, copies, axis=0)
        energy_old = np.repeat(energy_new, copies)
        if len(positions) > 10 * target:
            keep = rng.choice(len(positions), 10 * target, replace=False)
            positions, spins = positions[keep], spins[keep]
            log_old, force_old = log_old[keep], force_old[keep]
            energy_old = energy_old[keep]
        ens.M = len(positions)
        population[step] = len(positions)

        mixed = float(energy_old.mean())
        running += 0.01 * (mixed - running)
        trial_energy = running - (feedback / time_step) * math.log(
            len(positions) / target)
        if step >= burn_in:
            series[step - burn_in] = mixed
            growth[step - burn_in] = trial_energy

    steps = n_steps + burn_in
    return dict(series=series, energy=float(series.mean()),
                growth=float(growth.mean()),
                acceptance=acceptance / steps,
                spin_acceptance=spin_acceptance / steps if sample_spins else 0.0,
                population=population,
                mean_population=float(population[burn_in:].mean()),
                time_step=time_step, alpha=alpha, beta=beta, omega=omega,
                n_particles=n_particles, sample_spins=sample_spins)
