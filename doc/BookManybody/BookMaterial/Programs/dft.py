"""
Density functional theory: densities, functionals and the Kohn-Sham equations.

Companion code to chapter 8 of *Quantum mechanics for Many-particle Systems*.

Three things are done here.

First, the local density approximation is *derived* rather than quoted.  For
the three-dimensional electron gas the exchange energy per particle is
obtained from the one-body density matrix of the filled Fermi sea and checked
against the analytic 0.916/r_s of chapter 6.  The same construction is then
carried out numerically for a one-dimensional gas interacting through the
softened Coulomb potential used in chapters 2 and 6, which produces a
tabulated exchange functional for exactly that system.

Second, the Kohn-Sham equations are solved self-consistently in the same
oscillator basis as the Hartree-Fock calculation of chapter 6, so that
Hartree-Fock, Kohn-Sham and (for two particles) the exact answer can be
compared for one and the same Hamiltonian.

Third, the two characteristic errors of the local approximation are exhibited:
the exchange hole is too shallow and too short-ranged, and the self-interaction
that Hartree-Fock cancels exactly survives in the Kohn-Sham energy.

    ExchangeLDA3D   -- eps_x[n] for the three-dimensional electron gas
    ExchangeLDA1D   -- the same, built numerically for the model interaction
    ThomasFermi     -- the kinetic-energy functional and its failure
    KohnSham1D      -- a self-consistent Kohn-Sham solver
    exchange_hole   -- the pair correlation function of the uniform gas

Author: Morten Hjorth-Jensen
"""

from math import factorial, pi

import numpy as np
from numpy.polynomial.hermite import hermval
from scipy.linalg import eigh


# ---------------------------------------------------------------------------
#  The three-dimensional electron gas
# ---------------------------------------------------------------------------
class ExchangeLDA3D:
    """Exchange energy per particle of the homogeneous electron gas.

    The one-body density matrix of the filled Fermi sea is

        rho(s) = 3 n [sin(k_F s) - k_F s cos(k_F s)] / (k_F s)^3 ,

    and the exchange energy is

        E_x / V = -(1/4) int d^3s  v(s) |rho(s)|^2 ,

    with v(s) = 1/s in Hartree atomic units.  The integral can be done in
    closed form and gives the familiar

        eps_x(n) = -(3/4) (3/pi)^(1/3) n^(1/3) .
    """

    @staticmethod
    def eps_x(n):
        """Exchange energy per particle, in Hartree."""
        n = np.asarray(n, dtype=float)
        return -0.75 * (3.0 / pi) ** (1.0 / 3.0) * n ** (1.0 / 3.0)

    @staticmethod
    def v_x(n):
        """Exchange potential  d(n eps_x)/dn = (4/3) eps_x."""
        return (4.0 / 3.0) * ExchangeLDA3D.eps_x(n)

    @staticmethod
    def density_from_rs(rs):
        """n from the dimensionless radius r_s (in units of the Bohr radius)."""
        return 3.0 / (4.0 * pi * np.asarray(rs, dtype=float) ** 3)

    @staticmethod
    def coefficient_rydberg():
        """The 0.916 of chapter 6, obtained from the LDA expression."""
        return 0.75 * (3.0 / pi) ** (1.0 / 3.0) \
            * (3.0 / (4.0 * pi)) ** (1.0 / 3.0) * 2.0

    @staticmethod
    def kinetic_density(n):
        """Thomas-Fermi kinetic energy density, in Hartree per unit volume."""
        n = np.asarray(n, dtype=float)
        return 0.3 * (3.0 * pi ** 2) ** (2.0 / 3.0) * n ** (5.0 / 3.0)


def exchange_hole(y):
    """The exchange part of the pair correlation function of the 3D gas.

    With y = k_F s,

        g(s) - 1 = -(9/2) [ (sin y - y cos y) / y^3 ]^2 ,

    for a spin-unpolarised gas, so that g(0) = 1/2: two electrons of the same
    spin never coincide, two of opposite spin are uncorrelated at this level.
    """
    y = np.asarray(y, dtype=float)
    out = np.empty_like(y)
    small = np.abs(y) < 1e-8
    out[small] = 1.0 / 3.0
    ys = y[~small]
    out[~small] = (np.sin(ys) - ys * np.cos(ys)) / ys ** 3
    return 1.0 - 4.5 * out ** 2


# ---------------------------------------------------------------------------
#  A one-dimensional gas with the interaction of chapters 2 and 6
# ---------------------------------------------------------------------------
class ExchangeLDA1D:
    """Exchange functional of a 1D gas of spinless fermions.

    The interaction is the softened Coulomb potential of chapters 2 and 6,

        v(s) = strength / sqrt(s^2 + softening^2) ,

    for which no closed form exists.  We therefore do what is done in
    practice: solve the uniform gas once, tabulate eps_x(n), and interpolate.

    At density n the Fermi momentum of a spinless one-dimensional gas is
    k_F = pi n and the one-body density matrix is

        rho(s) = sin(k_F s) / (pi s) ,

    so that

        eps_x(n) = -(1/n) int_0^inf ds  v(s) [sin(k_F s) / (pi s)]^2 ,

    the factor two from the two half-lines cancelling the one half of the
    exchange integral.
    """

    def __init__(self, strength=1.0, softening=0.5,
                 n_max=3.0, n_points=241, cutoff=400.0, n_quad=200001):
        self.strength = strength
        self.softening = softening
        self.grid = np.linspace(1e-4, n_max, n_points)
        self.table = np.array([self._integrate(n) for n in self.grid])
        # a smooth interpolant, plus its derivative for the potential
        self._spline = None

    # ------------------------------------------------------------------
    def interaction(self, s):
        return self.strength / np.sqrt(s ** 2 + self.softening ** 2)

    def _integrate(self, n, cutoff=400.0, n_quad=200001):
        """eps_x(n) by Simpson integration on a fine grid."""
        k_f = pi * n
        s = np.linspace(0.0, cutoff, n_quad)
        ds = s[1] - s[0]
        rho = np.empty_like(s)
        rho[0] = k_f / pi
        rho[1:] = np.sin(k_f * s[1:]) / (pi * s[1:])
        f = self.interaction(s) * rho ** 2
        # composite Simpson
        integral = ds / 3.0 * (f[0] + f[-1]
                               + 4.0 * f[1:-1:2].sum()
                               + 2.0 * f[2:-2:2].sum())
        return -integral / n

    # ------------------------------------------------------------------
    def _build_spline(self):
        from scipy.interpolate import CubicSpline
        self._spline = CubicSpline(self.grid, self.grid * self.table)
        return self._spline

    def energy_density(self, n):
        """n * eps_x(n), the exchange energy per unit length."""
        if self._spline is None:
            self._build_spline()
        n = np.clip(np.asarray(n, dtype=float), 0.0, self.grid[-1])
        return np.where(n < self.grid[0], 0.0, self._spline(n))

    def eps_x(self, n):
        n = np.asarray(n, dtype=float)
        out = np.zeros_like(n)
        mask = n > self.grid[0]
        out[mask] = self.energy_density(n[mask]) / n[mask]
        return out

    def v_x(self, n):
        """d(n eps_x)/dn, the exchange potential."""
        if self._spline is None:
            self._build_spline()
        n = np.clip(np.asarray(n, dtype=float), 0.0, self.grid[-1])
        return np.where(n < self.grid[0], 0.0, self._spline(n, 1))


# ---------------------------------------------------------------------------
#  Thomas-Fermi
# ---------------------------------------------------------------------------
class ThomasFermi:
    """The Thomas-Fermi model for the one-dimensional trap.

    The kinetic functional of a spinless one-dimensional gas follows from the
    same Fermi-sea average that gives Thomas-Fermi in three dimensions,

        t[n] = (pi^2 / 6) n^3 ,

    since  T/L = int_{-k_F}^{k_F} dk/(2 pi) k^2/2 = k_F^3/(6 pi)
    and  k_F = pi n.  Minimising

        E[n] = int dx { t[n] + v_ext n } + E_H[n] + E_x[n]

    subject to a fixed particle number gives an algebraic equation for n(x),
    with no orbitals anywhere -- which is exactly why it has no shell
    structure.
    """

    @staticmethod
    def kinetic_density(n):
        return (pi ** 2 / 6.0) * np.asarray(n, dtype=float) ** 3

    @staticmethod
    def kinetic_potential(n):
        return 0.5 * pi ** 2 * np.asarray(n, dtype=float) ** 2


# ---------------------------------------------------------------------------
#  The trap: the same system as chapters 2 and 6
# ---------------------------------------------------------------------------
def oscillator_basis(n_orbitals=8, grid=None):
    """Harmonic-oscillator eigenfunctions on a grid, and the grid itself."""
    if grid is None:
        grid = np.linspace(-8.0, 8.0, 401)
    phi = np.zeros((n_orbitals, grid.size))
    for k in range(n_orbitals):
        coef = np.zeros(k + 1)
        coef[k] = 1.0
        norm = 1.0 / np.sqrt(2.0 ** k * factorial(k) * np.sqrt(pi))
        phi[k] = norm * hermval(grid, coef) * np.exp(-0.5 * grid ** 2)
    return grid, phi


class KohnSham1D:
    """Self-consistent Kohn-Sham for N spinless fermions in a 1D trap.

    The external potential is the harmonic trap, so the one-body part is
    diagonal in the oscillator basis with eigenvalues k + 1/2.  The
    interaction is the softened Coulomb potential of chapters 2 and 6.  At
    each iteration we build

        n(x)   = sum_i |psi_i(x)|^2 ,
        v_H(x) = int dy v(x, y) n(y) ,
        v_x(x) = v_x^unif( n(x) )        (the local approximation)

    project v_H + v_x onto the oscillator basis, diagonalise, and repeat.
    """

    def __init__(self, n_particles, n_orbitals=8, strength=1.0,
                 softening=0.5, lda=None, grid=None):
        self.N = n_particles
        self.n_orbitals = n_orbitals
        self.strength = strength
        self.softening = softening
        self.grid, self.phi = oscillator_basis(n_orbitals, grid)
        self.dx = self.grid[1] - self.grid[0]
        self.h0 = np.diag(np.arange(n_orbitals) + 0.5)
        self.interaction = strength / np.sqrt(
            (self.grid[:, None] - self.grid[None, :]) ** 2 + softening ** 2)
        self.lda = lda if lda is not None else ExchangeLDA1D(strength, softening)

    # ------------------------------------------------------------------
    def density(self, C):
        occupied = self.phi.T @ C[:, :self.N]        # (grid, N)
        return (occupied ** 2).sum(axis=1)

    def hartree_potential(self, n):
        return self.interaction @ n * self.dx

    def project(self, potential):
        """<a| potential |b> on the grid."""
        return (self.phi * potential) @ self.phi.T * self.dx

    # ------------------------------------------------------------------
    def hartree_energy(self, n):
        return 0.5 * float(n @ (self.interaction @ n)) * self.dx ** 2

    def exchange_energy(self, n):
        return float(np.sum(self.lda.energy_density(n)) * self.dx)

    def total_energy(self, eps, n):
        """E = sum_i eps_i - E_H + E_x - int v_x n."""
        band = float(eps[:self.N].sum())
        v_x = self.lda.v_x(n)
        return (band - self.hartree_energy(n) + self.exchange_energy(n)
                - float(n @ v_x) * self.dx)

    # ------------------------------------------------------------------
    def run(self, max_iter=400, tol=1e-10, mixing=0.3, verbose=False):
        C = np.eye(self.n_orbitals)
        n = self.density(C)
        history = []
        energy = np.nan
        for iteration in range(1, max_iter + 1):
            v_h = self.hartree_potential(n)
            v_x = self.lda.v_x(n)
            F = self.h0 + self.project(v_h + v_x)
            eps, C = eigh(F)
            n_new = self.density(C)
            drift = float(np.abs(n_new - n).max())
            n = (1.0 - mixing) * n + mixing * n_new
            energy = self.total_energy(eps, n)
            history.append((iteration, energy, drift))
            if verbose:
                print(f"   {iteration:4d}  E = {energy:14.10f}   "
                      f"max |dn| = {drift:.3e}")
            if drift < tol:
                break
        self.C, self.eps, self.n, self.E = C, eps, n, energy
        self.iterations = iteration
        self.history = history
        return self.E

    # ------------------------------------------------------------------
    def components(self):
        """The pieces of the total energy, for inspection."""
        occupied = self.phi.T @ self.C[:, :self.N]
        kinetic = float(np.einsum("ai,ab,bi->", self.C[:, :self.N],
                                  self.h0, self.C[:, :self.N]))
        return dict(band=float(self.eps[:self.N].sum()),
                    one_body=kinetic,
                    hartree=self.hartree_energy(self.n),
                    exchange=self.exchange_energy(self.n),
                    total=self.E,
                    particles=float(self.n.sum() * self.dx))


# ---------------------------------------------------------------------------
#  Hartree-Fock in the same basis, for comparison
# ---------------------------------------------------------------------------
def trap_two_body(n_orbitals=8, strength=1.0, softening=0.5, grid=None):
    """Antisymmetrised matrix elements of the softened Coulomb interaction."""
    grid, phi = oscillator_basis(n_orbitals, grid)
    dx = grid[1] - grid[0]
    interaction = strength / np.sqrt(
        (grid[:, None] - grid[None, :]) ** 2 + softening ** 2)
    direct = np.einsum("ax,by,xy,cx,dy->abcd",
                       phi, phi, interaction * dx * dx, phi, phi,
                       optimize=True)
    return direct - direct.transpose(0, 1, 3, 2)


def hartree_fock(n_particles, n_orbitals=8, strength=1.0, softening=0.5,
                 max_iter=400, tol=1e-12):
    """The Hartree-Fock calculation of chapter 6, repeated here."""
    h0 = np.diag(np.arange(n_orbitals) + 0.5)
    vbar = trap_two_body(n_orbitals, strength, softening)
    C = np.eye(n_orbitals)
    previous = None
    for iteration in range(1, max_iter + 1):
        rho = C[:, :n_particles] @ C[:, :n_particles].T
        F = h0 + np.einsum("acbd,dc->ab", vbar, rho)
        eps, C = eigh(F)
        if previous is not None and abs(eps.sum() - previous) < tol:
            break
        previous = eps.sum()
    occ = C[:, :n_particles]
    one = np.einsum("ab,ab->", h0, occ @ occ.T)
    two = 0.5 * np.einsum("ai,bj,abcd,ci,dj->", occ, occ, vbar, occ, occ)
    return dict(energy=one + two, eps=eps, C=C, iterations=iteration,
                vbar=vbar)


def exact_two_particle(n_orbitals=8, strength=1.0, softening=0.5):
    """Full configuration interaction for two spinless fermions."""
    from itertools import combinations
    vbar = trap_two_body(n_orbitals, strength, softening)
    sp = np.arange(n_orbitals) + 0.5
    pairs = list(combinations(range(n_orbitals), 2))
    dim = len(pairs)
    H = np.zeros((dim, dim))
    for a, (i, j) in enumerate(pairs):
        H[a, a] += sp[i] + sp[j]
        for b, (k, l) in enumerate(pairs):
            H[a, b] += vbar[i, j, k, l]
    return float(np.linalg.eigvalsh(H)[0]), dim


# ---------------------------------------------------------------------------
#  Demonstrations
# ---------------------------------------------------------------------------
def demo_lda_from_the_gas():
    print("=" * 74)
    print("1. Where the local density approximation comes from")
    print("=" * 74)
    print("The exchange energy per particle of the three-dimensional gas,")
    print("   eps_x(n) = -(3/4)(3/pi)^(1/3) n^(1/3),")
    print("is the whole content of the local density approximation for")
    print("exchange.  Evaluated on the uniform gas it must reproduce the")
    print("0.916/r_s of chapter 6, and it does:")
    print()
    print(f"{'r_s':>8s} {'eps_x [Ha]':>14s} {'eps_x [Ry]':>14s} "
          f"{'-0.916/r_s':>14s}")
    for rs in (1.0, 2.0, 4.0, 4.8253, 6.0):
        n = ExchangeLDA3D.density_from_rs(rs)
        e = float(ExchangeLDA3D.eps_x(n))
        print(f"{rs:8.4f} {e:14.6f} {2 * e:14.6f} {-0.916 / rs:14.6f}")
    print()
    print(f"the coefficient itself: {ExchangeLDA3D.coefficient_rydberg():.6f}")
    print("against the 0.9163 exchange coefficient of chapter 6.  This is not a")
    print("coincidence but a tautology: the LDA is *defined* so that it is")
    print("exact for the uniform gas.  What is not obvious, and what makes it")
    print("useful, is that it remains good when the density varies.")
    print()
    print("The Thomas-Fermi kinetic energy density comes from the same")
    print("calculation, and equals the kinetic term of chapter 6 divided by")
    print("the volume:")
    print(f"{'r_s':>8s} {'t_TF [Ha/a0^3]':>18s} {'3/5 eps_F [Ha]':>18s}")
    for rs in (1.0, 2.0, 4.0):
        n = ExchangeLDA3D.density_from_rs(rs)
        t = float(ExchangeLDA3D.kinetic_density(n))
        k_f = (3.0 * pi ** 2 * n) ** (1.0 / 3.0)
        print(f"{rs:8.4f} {t:18.8f} {0.6 * 0.5 * k_f ** 2 * n:18.8f}")


def demo_exchange_hole():
    print("=" * 74)
    print("2. The exchange hole")
    print("=" * 74)
    print("Exchange is not a mysterious correction: it is the statement that")
    print("each electron carries around it a region depleted of electrons of")
    print("the same spin.  For the uniform gas the pair correlation function")
    print("is known exactly at the Hartree-Fock level,")
    print("   g(s) = 1 - (9/2)[(sin y - y cos y)/y^3]^2,   y = k_F s.")
    print()
    print(f"{'k_F s':>8s} {'g(s)':>12s}   interpretation")
    for y in (0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0):
        g = float(exchange_hole(np.array([y]))[0])
        note = ""
        if y == 0.0:
            note = "g(0) = 1/2: same spins never coincide"
        elif y > 10:
            note = "the hole has healed"
        print(f"{y:8.2f} {g:12.6f}   {note}")
    print()
    print("The hole integrates to exactly one missing electron,")
    print("   n int d^3s [g(s) - 1] = (4 / 3 pi) int_0^inf dy y^2 [g(y) - 1]")
    print("                         = -1,")
    print("which is the sum rule that any decent approximate functional must")
    print("respect.  Checking it numerically:")
    y = np.linspace(1e-8, 400.0, 4000001)
    integrand = y ** 2 * (exchange_hole(y) - 1.0)
    total = np.trapezoid(integrand, y) * 4.0 / (3.0 * pi)
    print(f"   (4 / 3 pi) int dy y^2 [g(y) - 1] = {total:.6f}")
    print()
    print("The exact value of the integral is -pi/6 times 9/2, and the")
    print("prefactor turns it into exactly -1.  Local functionals are built")
    print("to satisfy this sum rule even when they get the shape of the hole")
    print("wrong, which is a large part of why they work at all.")


def demo_lda_1d():
    print("=" * 74)
    print("3. Building a local functional for the model interaction")
    print("=" * 74)
    print("The trap of chapters 2 and 6 has a softened Coulomb interaction in")
    print("one dimension, for which no closed-form exchange functional")
    print("exists.  So we do what is done in practice: solve the uniform gas")
    print("once, tabulate, and interpolate.  For a spinless one-dimensional")
    print("gas k_F = pi n, the density matrix is sin(k_F s)/(pi s), and")
    print("   eps_x(n) = -(1/n) int_0^inf ds v(s) [sin(k_F s)/(pi s)]^2.")
    print()
    lda = ExchangeLDA1D()
    print(f"{'n':>8s} {'eps_x(n)':>14s} {'v_x(n)':>14s}")
    for n in (0.05, 0.1, 0.2, 0.4, 0.8, 1.5, 2.5):
        arr = np.array([n])
        print(f"{n:8.3f} {lda.eps_x(arr)[0]:14.8f} "
              f"{lda.v_x(arr)[0]:14.8f}")
    print()
    print("Two checks.  The interpolated energy density must reproduce the")
    print("integrated values it was built from, and the potential must be")
    print("the derivative of the energy density:")
    test = np.array([0.15, 0.35, 0.9, 1.8])
    direct = np.array([lda._integrate(float(x)) for x in test])
    print(f"   max |interpolated - integrated| = "
          f"{np.abs(lda.eps_x(test) - direct).max():.2e}")
    h = 1e-5
    numeric = (lda.energy_density(test + h)
               - lda.energy_density(test - h)) / (2 * h)
    print(f"   max |v_x - d(n eps_x)/dn|       = "
          f"{np.abs(lda.v_x(test) - numeric).max():.2e}")


def demo_kohn_sham():
    print("=" * 74)
    print("4. Solving the Kohn-Sham equations")
    print("=" * 74)
    print("The same trap, the same basis of eight oscillator orbitals and the")
    print("same interaction as the Hartree-Fock calculation of chapter 6, so")
    print("that the two can be compared directly.  For two particles the")
    print("exact answer is available as well.")
    print()
    lda = ExchangeLDA1D()
    print(f"{'N':>3s} {'E(Kohn-Sham)':>15s} {'E(Hartree-Fock)':>17s} "
          f"{'E(exact)':>12s} {'KS iter':>9s} {'HF iter':>9s}")
    for N in (1, 2, 3, 4):
        ks = KohnSham1D(N, lda=lda)
        ks.run()
        hf = hartree_fock(N)
        if N == 2:
            exact, _ = exact_two_particle()
            exact_str = f"{exact:12.8f}"
        elif N == 1:
            exact_str = f"{0.5:12.8f}"
        else:
            exact_str = f"{'--':>12s}"
        print(f"{N:3d} {ks.E:15.8f} {hf['energy']:17.8f} {exact_str} "
              f"{ks.iterations:9d} {hf['iterations']:9d}")
    print()
    print("Kohn-Sham lies close to Hartree-Fock but not on top of it, and the")
    print("difference is not correlation -- we have put no correlation")
    print("functional in.  It is the error of treating exchange locally, plus")
    print("the self-interaction discussed below.")


def demo_self_interaction():
    print("=" * 74)
    print("5. The self-interaction error")
    print("=" * 74)
    print("For one particle there is nothing to interact with, so the exact")
    print("energy is the oscillator ground state, 0.5.  Hartree-Fock gets")
    print("this right identically: the direct and exchange terms cancel term")
    print("by term.  A local exchange functional cannot, because it knows")
    print("only the density, and the density of one electron looks exactly")
    print("like the density of half of two.")
    print()
    lda = ExchangeLDA1D()
    ks = KohnSham1D(1, lda=lda)
    ks.run()
    parts = ks.components()
    hf = hartree_fock(1)
    print(f"   exact                     0.50000000")
    print(f"   Hartree-Fock         {hf['energy']:15.8f}")
    print(f"   Kohn-Sham (LDA x)    {ks.E:15.8f}")
    print()
    print("   the Kohn-Sham pieces:")
    print(f"      Hartree energy    {parts['hartree']:15.8f}")
    print(f"      exchange energy   {parts['exchange']:15.8f}")
    print(f"      sum               {parts['hartree'] + parts['exchange']:15.8f}"
          "   (should vanish)")
    print(f"      particle number   {parts['particles']:15.8f}")
    print()
    print("The residue is the self-interaction error.  It is the reason plain")
    print("LDA overbinds anions, underestimates barriers and gives band gaps")
    print("that are too small, and it is what self-interaction corrections and")
    print("hybrid functionals -- which put back a fraction of exact exchange --")
    print("are designed to remove.")


def demo_local_approximation():
    print("=" * 74)
    print("6. How good is the local approximation?")
    print("=" * 74)
    print("We can isolate the error of locality by comparing, on the *same*")
    print("density, the exact exchange energy of a determinant with the local")
    print("estimate.  Take the Kohn-Sham orbitals, compute both, and look at")
    print("the ratio.")
    print()
    lda = ExchangeLDA1D()
    print(f"{'N':>3s} {'E_x (exact)':>15s} {'E_x (LDA)':>15s} {'ratio':>9s}")
    for N in (1, 2, 3, 4):
        ks = KohnSham1D(N, lda=lda)
        ks.run()
        vbar = trap_two_body(ks.n_orbitals)
        occ = ks.C[:, :N]
        direct_plus_exchange = 0.5 * np.einsum(
            "ai,bj,abcd,ci,dj->", occ, occ, vbar, occ, occ)
        hartree = ks.hartree_energy(ks.n)
        exact_x = direct_plus_exchange - hartree
        local_x = ks.exchange_energy(ks.n)
        print(f"{N:3d} {exact_x:15.8f} {local_x:15.8f} "
              f"{local_x / exact_x:9.4f}")
    print()
    print("The local approximation recovers most of the exchange energy, and")
    print("it recovers more of it as the particle number grows: the ratio")
    print("climbs from 0.92 to 0.98 between one and four particles.  That is")
    print("exactly what the construction promises -- the density becomes")
    print("smoother on the scale of the local Fermi wavelength, so the")
    print("uniform-gas estimate becomes more appropriate.  For one particle")
    print("it is worst, because there the exchange energy is nothing but the")
    print("cancellation of a spurious self-interaction.  In three dimensions")
    print("the corresponding ratio is close to 0.9 for atoms, which is why")
    print("gradient corrections were introduced.")


def demo_thomas_fermi():
    print("=" * 74)
    print("7. Thomas-Fermi, and why orbitals came back")
    print("=" * 74)
    print("The Thomas-Fermi model applies the local approximation to the")
    print("kinetic energy as well.  For the one-dimensional trap the kinetic")
    print("density is t[n] = (pi^2/6) n^3, and minimising")
    print("   E[n] = int dx {t[n] + x^2 n / 2} + E_H[n] + E_x[n]")
    print("gives an algebraic equation for n(x) with no orbitals at all.")
    print()
    print("The trouble is visible already for non-interacting particles,")
    print("where the exact kinetic energy is known:")
    print(f"{'N':>3s} {'T (exact)':>12s} {'T (Thomas-Fermi)':>18s} "
          f"{'ratio':>9s}")
    grid, phi = oscillator_basis(8)
    dx = grid[1] - grid[0]
    for N in (1, 2, 3, 4, 6, 8):
        n = (phi[:N] ** 2).sum(axis=0)
        exact = sum(k + 0.5 for k in range(N)) \
            - 0.5 * float(n @ grid ** 2) * dx
        tf = float(ThomasFermi.kinetic_density(n).sum()) * dx
        print(f"{N:3d} {exact:12.6f} {tf:18.6f} {tf / exact:9.4f}")
    print()
    print("Thomas-Fermi overshoots, badly for one particle and by less than")
    print("one per cent by eight, approaching the exact value from above as")
    print("the density smooths out.  But the error in magnitude is not the")
    print("real problem.  What the model misses is structure: the density it")
    print("produces is smooth and featureless, with no shells whatever, it")
    print("has no discrete spectrum to speak of, and it famously cannot bind")
    print("molecules -- Teller proved this.  The Kohn-Sham construction exists")
    print("precisely to avoid this.  By reintroducing orbitals for the")
    print("kinetic energy alone -- the one term we cannot approximate")
    print("locally -- it keeps shell structure exactly while leaving only the")
    print("small remainder to a local functional.")


def _demo():
    for f in (demo_lda_from_the_gas, demo_exchange_hole, demo_lda_1d,
              demo_kohn_sham, demo_self_interaction,
              demo_local_approximation, demo_thomas_fermi):
        f()
        print()


if __name__ == "__main__":
    _demo()
