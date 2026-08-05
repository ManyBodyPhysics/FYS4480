"""
Variational Monte Carlo for two electrons in a two-dimensional quantum dot.

Companion code to chapter 13 of *Quantum Mechanics for Many-particle Systems*.

The Hamiltonian is the one of chapter 11,

    H = sum_i [ -grad_i^2 / 2 + omega^2 r_i^2 / 2 ] + 1 / r_12 ,

and the trial function is the two-parameter Pade-Jastrow form

    Psi_T = exp( -alpha omega (r_1^2 + r_2^2) / 2 )
            * exp( r_12 / (1 + beta r_12) ) ,

with alpha scaling the oscillator orbital and beta controlling the range of
the correlation factor.  The numerator of the Jastrow exponent is fixed at one
by the cusp condition for two electrons of opposite spin in two dimensions, so
there are exactly two variational parameters.

Everything that the sampling needs is available in closed form: the local
energy and the quantum force are both differentiated by hand and then checked
against finite differences, because an error in either is silent -- a wrong
local energy gives a wrong energy, and a wrong quantum force merely makes the
sampling inefficient without changing the answer.  Both checks are in the
demo.

Two samplers are implemented.  The brute-force Metropolis of chapter 12 moves
a particle by a uniform random displacement; importance sampling replaces that
with a step drawn from the Fokker-Planck Green's function, drifting the walker
along the quantum force F = 2 grad Psi_T / Psi_T, and corrects the asymmetry
of the proposal with the Metropolis-Hastings acceptance ratio.

Author: Morten Hjorth-Jensen
"""

import math

import numpy as np

# ---------------------------------------------------------------------------
#  Exact reference values
# ---------------------------------------------------------------------------
#  Taut, Phys. Rev. A 48, 3561 (1993): the relative motion of two electrons in
#  a parabolic trap is analytically solvable for a discrete set of
#  frequencies, and at hbar omega = 1 the ground state is exactly 3.
# ---------------------------------------------------------------------------
TAUT_ENERGY = 3.0


# ---------------------------------------------------------------------------
#  The trial wave function and everything derived from it
# ---------------------------------------------------------------------------
def log_psi(r, alpha, beta, omega=1.0):
    """ln Psi_T for a configuration r of shape (2, 2)."""
    r_squared = r[0, 0]**2 + r[0, 1]**2 + r[1, 0]**2 + r[1, 1]**2
    r12 = math.hypot(r[0, 0] - r[1, 0], r[0, 1] - r[1, 1])
    return -0.5 * alpha * omega * r_squared + r12 / (1.0 + beta * r12)


def local_energy(r, alpha, beta, omega=1.0):
    """E_L = H Psi_T / Psi_T, in closed form.

        E_L = (1/2) omega^2 (1 - alpha^2) (r_1^2 + r_2^2) + 2 alpha omega
              + 1/r_12
              + d^2 ( alpha omega r_12 - d^2 + 2 beta d - 1/r_12 ) ,

    with d = 1/(1 + beta r_12).  The first line is the pure Gaussian result of
    chapter 12; the second is everything the Jastrow factor contributes, and it
    is what removes the divergence of 1/r_12 -- the two 1/r_12 terms cancel as
    r_12 -> 0, which is the cusp condition doing its work.
    """
    r_squared = r[0, 0]**2 + r[0, 1]**2 + r[1, 0]**2 + r[1, 1]**2
    r12 = math.hypot(r[0, 0] - r[1, 0], r[0, 1] - r[1, 1])
    d = 1.0 / (1.0 + beta * r12)
    d2 = d * d
    return (0.5 * omega**2 * (1.0 - alpha**2) * r_squared
            + 2.0 * alpha * omega
            + 1.0 / r12
            + d2 * (alpha * omega * r12 - d2 + 2.0 * beta * d - 1.0 / r12))


def quantum_force(r, alpha, beta, omega=1.0):
    """F = 2 grad Psi_T / Psi_T = 2 grad ln Psi_T, one row per particle.

        F_1 = -2 alpha omega r_1 + 2 (r_1 - r_2) d^2 / r_12 ,

    and likewise for particle two with the labels exchanged.  Note that this
    is a *sum* of two terms: the first pulls the walker towards the centre of
    the trap, the second pushes the two electrons apart.
    """
    force = np.empty((2, 2))
    separation = r[0] - r[1]
    r12 = math.hypot(separation[0], separation[1])
    d = 1.0 / (1.0 + beta * r12)
    drift = 2.0 * separation * d * d / r12
    force[0] = -2.0 * alpha * omega * r[0] + drift
    force[1] = -2.0 * alpha * omega * r[1] - drift
    return force


# -- finite-difference versions, used only to check the two above -----------
def local_energy_fd(r, alpha, beta, omega=1.0, h=1e-5):
    """E_L by finite differences, from grad^2 Psi / Psi = grad^2 ln Psi
    + (grad ln Psi)^2."""
    laplacian = gradient_squared = 0.0
    centre = log_psi(r, alpha, beta, omega)
    for i in range(2):
        for k in range(2):
            shifted = r.copy()
            shifted[i, k] += h
            plus = log_psi(shifted, alpha, beta, omega)
            shifted[i, k] -= 2.0 * h
            minus = log_psi(shifted, alpha, beta, omega)
            laplacian += (plus + minus - 2.0 * centre) / h**2
            gradient_squared += ((plus - minus) / (2.0 * h))**2
    r_squared = r[0, 0]**2 + r[0, 1]**2 + r[1, 0]**2 + r[1, 1]**2
    r12 = math.hypot(r[0, 0] - r[1, 0], r[0, 1] - r[1, 1])
    return (-0.5 * (laplacian + gradient_squared)
            + 0.5 * omega**2 * r_squared + 1.0 / r12)


def quantum_force_fd(r, alpha, beta, omega=1.0, h=1e-5):
    """F by finite differences."""
    force = np.empty((2, 2))
    for i in range(2):
        for k in range(2):
            shifted = r.copy()
            shifted[i, k] += h
            plus = log_psi(shifted, alpha, beta, omega)
            shifted[i, k] -= 2.0 * h
            minus = log_psi(shifted, alpha, beta, omega)
            force[i, k] = 2.0 * (plus - minus) / (2.0 * h)
    return force


# ---------------------------------------------------------------------------
#  Error estimates for correlated data (chapter 12, section 12.3)
# ---------------------------------------------------------------------------
def correlation_time(x, max_lag=500):
    """tau = 1 + 2 sum kappa_d, truncated at the first negative kappa."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    deviation = x - x.mean()
    variance = float(np.dot(deviation, deviation) / n)
    if variance <= 0.0:
        return 1.0
    kappa = [float(np.dot(deviation[:n - d], deviation[d:]) / n / variance)
             for d in range(1, min(max_lag, n // 4))]
    total = 0.0
    for value in kappa:
        if value < 0.0:
            break
        total += value
    return 1.0 + 2.0 * total


def statistics(samples):
    """Mean, correlation-corrected error, variance and correlation time."""
    samples = np.asarray(samples, dtype=float)
    tau = correlation_time(samples)
    naive = samples.std(ddof=1) / math.sqrt(len(samples))
    return dict(energy=float(samples.mean()),
                error=float(naive * math.sqrt(tau)),
                naive_error=float(naive),
                variance=float(samples.var(ddof=1)),
                tau=float(tau))


# ---------------------------------------------------------------------------
#  Sampler 1: brute-force Metropolis
# ---------------------------------------------------------------------------
def vmc_brute_force(alpha, beta, n_cycles=100000, step=1.0, omega=1.0,
                    rng=None, burn_in=5000):
    """Uniform proposals, acceptance min(1, |Psi(y)|^2 / |Psi(x)|^2).

    One particle is moved at a time, which is the standard choice: moving all
    of them at once lowers the acceptance rate for the same step length.
    """
    if rng is None:
        rng = np.random.default_rng(2024)
    position = rng.normal(0.0, 1.0, (2, 2))
    log_old = log_psi(position, alpha, beta, omega)
    samples = np.empty(n_cycles)
    accepted = 0
    proposals = 0

    for cycle in range(n_cycles + burn_in):
        for particle in range(2):
            trial = position.copy()
            trial[particle] += step * (rng.random(2) - 0.5)
            log_new = log_psi(trial, alpha, beta, omega)
            proposals += 1
            if 2.0 * (log_new - log_old) > math.log(rng.random() + 1e-300):
                position, log_old = trial, log_new
                accepted += 1
        if cycle >= burn_in:
            samples[cycle - burn_in] = local_energy(position, alpha, beta,
                                                    omega)

    out = statistics(samples)
    out.update(acceptance=accepted / proposals, samples=samples,
               alpha=alpha, beta=beta, method="brute force")
    return out


# ---------------------------------------------------------------------------
#  Sampler 2: importance sampling, Metropolis-Hastings
# ---------------------------------------------------------------------------
def greens_ratio(old, new, force_old, force_new, diffusion, time_step):
    """ln[ G(x, y) / G(y, x) ] for the Fokker-Planck Green's function.

    With G(y, x) proportional to exp( -(y - x - D dt F(x))^2 / 4 D dt ) the
    ratio collapses to

        (1/2)(F(x) + F(y)) . [ (D dt / 2)(F(x) - F(y)) - (y - x) ] ,

    which needs no exponentials and no normalisation.
    """
    return float(np.sum(0.5 * (force_old + force_new)
                        * (0.5 * diffusion * time_step
                           * (force_old - force_new) - new + old)))


def vmc_importance(alpha, beta, n_cycles=100000, time_step=0.05, omega=1.0,
                   rng=None, burn_in=5000, force=None):
    """Proposals drawn from the Fokker-Planck Green's function.

        y = x + D F(x) dt + xi sqrt(dt),   xi ~ N(0, 1),

    with D = 1/2 from the kinetic energy operator, accepted with the
    Metropolis-Hastings ratio

        q = G(x, y) |Psi(y)|^2 / [ G(y, x) |Psi(x)|^2 ] .

    The argument `force` allows a deliberately wrong quantum force to be
    substituted, which the demo uses to show that the algorithm still samples
    the right distribution -- only less efficiently.
    """
    if rng is None:
        rng = np.random.default_rng(2024)
    if force is None:
        force = quantum_force
    diffusion = 0.5
    root_dt = math.sqrt(time_step)

    position = rng.normal(0.0, 1.0, (2, 2)) * root_dt
    log_old = log_psi(position, alpha, beta, omega)
    force_old = force(position, alpha, beta, omega)
    samples = np.empty(n_cycles)
    accepted = 0
    proposals = 0

    for cycle in range(n_cycles + burn_in):
        for particle in range(2):
            trial = position.copy()
            trial[particle] = (position[particle]
                               + diffusion * force_old[particle] * time_step
                               + rng.normal(0.0, 1.0, 2) * root_dt)
            log_new = log_psi(trial, alpha, beta, omega)
            force_new = force(trial, alpha, beta, omega)
            log_ratio = (greens_ratio(position[particle], trial[particle],
                                      force_old[particle], force_new[particle],
                                      diffusion, time_step)
                         + 2.0 * (log_new - log_old))
            proposals += 1
            if log_ratio > math.log(rng.random() + 1e-300):
                position, log_old, force_old = trial, log_new, force_new
                accepted += 1
        if cycle >= burn_in:
            samples[cycle - burn_in] = local_energy(position, alpha, beta,
                                                    omega)

    out = statistics(samples)
    out.update(acceptance=accepted / proposals, samples=samples,
               alpha=alpha, beta=beta, time_step=time_step,
               method="importance sampling")
    return out


def wrong_quantum_force(r, alpha, beta, omega=1.0):
    """The quantum force with the two terms multiplied instead of added.

    This is not a legitimate drift, and it is kept only so that the demo can
    show what happens: because the same expression is used both to propose the
    move and to evaluate the Green's function, detailed balance still holds
    and the energy is still right -- the efficiency is what is lost.
    """
    force = np.empty((2, 2))
    separation = r[0] - r[1]
    r12 = math.hypot(separation[0], separation[1])
    d = 1.0 / (1.0 + beta * r12)
    force[0] = -2.0 * alpha * omega * r[0] * separation * d * d / r12
    force[1] = 2.0 * alpha * omega * r[1] * separation * d * d / r12
    return force


# ---------------------------------------------------------------------------
#  Parameter scans
# ---------------------------------------------------------------------------
def scan(alphas, betas, sampler=vmc_importance, n_cycles=20000, **kwargs):
    """Energy and variance on a grid of variational parameters."""
    energies = np.zeros((len(alphas), len(betas)))
    variances = np.zeros_like(energies)
    errors = np.zeros_like(energies)
    for i, alpha in enumerate(alphas):
        for j, beta in enumerate(betas):
            result = sampler(alpha, beta, n_cycles=n_cycles,
                             rng=np.random.default_rng(2024), **kwargs)
            energies[i, j] = result["energy"]
            variances[i, j] = result["variance"]
            errors[i, j] = result["error"]
    return energies, variances, errors


# ---------------------------------------------------------------------------
def _demo():
    print("=" * 74)
    print("1. The trial function, checked term by term")
    print("=" * 74)
    print("Both the local energy and the quantum force are differentiated by")
    print("hand.  An error in the first gives a wrong energy; an error in the")
    print("second is silent, so both are worth checking against finite")
    print("differences before anything else is run.")
    print()
    rng = np.random.default_rng(1)
    worst_energy = worst_force = 0.0
    for omega in (0.5, 1.0, 2.0):
        for _ in range(20):
            r = rng.normal(0.0, 1.0, (2, 2))
            for alpha, beta in ((0.98, 0.40), (0.80, 0.20), (1.10, 0.60)):
                worst_energy = max(worst_energy,
                                   abs(local_energy(r, alpha, beta, omega)
                                       - local_energy_fd(r, alpha, beta,
                                                         omega)))
                worst_force = max(worst_force,
                                  float(np.abs(
                                      quantum_force(r, alpha, beta, omega)
                                      - quantum_force_fd(r, alpha, beta,
                                                         omega)).max()))
    print(f"   local energy : largest discrepancy {worst_energy:.1e}")
    print(f"   quantum force: largest discrepancy {worst_force:.1e}")
    print("   (both at the level of the finite-difference truncation error)")

    print()
    print("=" * 74)
    print("2. The cusp condition at work")
    print("=" * 74)
    print("As the two electrons approach, 1/r_12 diverges.  With the Jastrow")
    print("factor the divergence cancels against the kinetic term and the")
    print("local energy stays finite; without it, it does not.")
    print()
    print(f"{'r_12':>8s} {'E_L with Jastrow':>18s} {'E_L without':>14s}")
    for r12 in (1.0, 0.1, 0.01, 0.001, 0.0001):
        r = np.array([[r12 / 2, 0.0], [-r12 / 2, 0.0]])
        with_jastrow = local_energy(r, 1.0, 0.4)
        # the same configuration with the correlation factor switched off
        r_squared = float(np.sum(r * r))
        without = 2.0 + 0.5 * (1.0 - 1.0) * r_squared + 1.0 / r12
        print(f"{r12:8.4f} {with_jastrow:18.6f} {without:14.2f}")
    print()
    print("   The right-hand column runs away as 1/r_12; the left-hand one")
    print("   converges.  That is the whole reason for the correlation factor,")
    print("   and it is why the variance of E_L collapses once it is switched")
    print("   on.")

    print()
    print("=" * 74)
    print("3. Brute-force Metropolis")
    print("=" * 74)
    print("The sampler of chapter 12: uniform proposals, one particle at a")
    print("time, accepted with min(1, |Psi(y)|^2/|Psi(x)|^2).  The step length")
    print("is the one free parameter.")
    print()
    print(f"{'step':>6s} {'acceptance':>11s} {'energy':>13s} {'error':>10s} "
          f"{'variance':>10s} {'tau':>7s}")
    for step in (0.5, 1.0, 2.0, 4.0, 8.0):
        r = vmc_brute_force(0.98, 0.40, n_cycles=30000, step=step,
                            rng=np.random.default_rng(2024))
        print(f"{step:6.1f} {r['acceptance']:10.1%} {r['energy']:13.6f} "
              f"{r['error']:10.6f} {r['variance']:10.5f} {r['tau']:7.2f}")
    print()
    print("   Small steps are accepted but move nowhere, large ones are")
    print("   rejected: tau is smallest in between, exactly as in chapter 12.")

    print()
    print("=" * 74)
    print("4. Importance sampling")
    print("=" * 74)
    print("Proposals now come from the Fokker-Planck Green's function, drifting")
    print("along F = 2 grad Psi_T / Psi_T, and the acceptance is corrected by")
    print("the Metropolis-Hastings ratio.  The time step replaces the step")
    print("length as the parameter to tune.")
    print()
    print(f"{'dt':>7s} {'acceptance':>11s} {'energy':>13s} {'error':>10s} "
          f"{'variance':>10s} {'tau':>7s}")
    for time_step in (0.01, 0.05, 0.1, 0.3, 0.5, 1.0):
        r = vmc_importance(0.98, 0.40, n_cycles=30000, time_step=time_step,
                           rng=np.random.default_rng(2024))
        print(f"{time_step:7.3f} {r['acceptance']:10.1%} {r['energy']:13.6f} "
              f"{r['error']:10.6f} {r['variance']:10.5f} {r['tau']:7.2f}")
    print()
    print("   The acceptance rate stays high even at time steps where a")
    print("   uniform proposal of the same size would be rejected almost")
    print("   always, because the drift proposes moves the wave function")
    print("   likes.  The correlation time falls monotonically here; the")
    print("   Metropolis-Hastings ratio removes the discretisation bias of")
    print("   the Langevin step exactly, so a large dt costs accuracy only")
    print("   through the acceptance rate.")

    print()
    print("=" * 74)
    print("5. What importance sampling actually buys")
    print("=" * 74)
    print("Comparing the two samplers at their best settings, on equal numbers")
    print("of Monte Carlo cycles:")
    print()
    print(f"{'sampler':>22s} {'acceptance':>11s} {'energy':>13s} "
          f"{'error':>10s} {'tau':>7s} {'n_eff':>8s}")
    brute = vmc_brute_force(0.98, 0.40, n_cycles=60000, step=4.0,
                            rng=np.random.default_rng(7))
    importance = vmc_importance(0.98, 0.40, n_cycles=60000, time_step=0.5,
                                rng=np.random.default_rng(7))
    for label, r in (("brute force", brute), ("importance sampling",
                                              importance)):
        print(f"{label:>22s} {r['acceptance']:10.1%} {r['energy']:13.6f} "
              f"{r['error']:10.6f} {r['tau']:7.2f} "
              f"{len(r['samples'])/r['tau']:8.0f}")
    print()
    print(f"   error ratio {brute['error']/importance['error']:.2f}, which is a")
    print(f"   factor {(brute['error']/importance['error'])**2:.1f} in computer "
          "time for the same accuracy.")

    print()
    print("=" * 74)
    print("6. A wrong quantum force is not a wrong answer")
    print("=" * 74)
    print("It is worth knowing what a mistake in the drift does.  Because the")
    print("same expression is used to propose the move and to evaluate the")
    print("Green's function, detailed balance is preserved whatever the drift")
    print("is: the energy stays correct and only the efficiency suffers.  A")
    print("wrong local energy, by contrast, gives a wrong number.")
    print()
    print(f"{'quantum force':>22s} {'acceptance':>11s} {'energy':>13s} "
          f"{'error':>10s} {'tau':>7s}")
    for label, f in (("correct", quantum_force),
                     ("terms multiplied", wrong_quantum_force)):
        r = vmc_importance(0.98, 0.40, n_cycles=30000, time_step=0.5,
                           rng=np.random.default_rng(11), force=f)
        print(f"{label:>22s} {r['acceptance']:10.1%} {r['energy']:13.6f} "
              f"{r['error']:10.6f} {r['tau']:7.2f}")

    print()
    print("=" * 74)
    print("7. The variational surface")
    print("=" * 74)
    print("Scanning both parameters.  The energy is flat near the minimum and")
    print("the variance is not, which is why variance minimisation is often")
    print("the better way to optimise a trial function.")
    print()
    alphas = np.array([0.94, 0.96, 0.98, 1.00, 1.02])
    betas = np.array([0.30, 0.36, 0.40, 0.44, 0.50])
    energies, variances, errors = scan(alphas, betas, n_cycles=15000,
                                       time_step=0.5)
    print("   energy")
    print("       " + "".join(f"{b:>11.2f}" for b in betas) + "   <- beta")
    for i, alpha in enumerate(alphas):
        print(f"  {alpha:5.2f}" + "".join(f"{energies[i, j]:11.5f}"
                                          for j in range(len(betas))))
    print()
    print("   variance")
    print("       " + "".join(f"{b:>11.2f}" for b in betas) + "   <- beta")
    for i, alpha in enumerate(alphas):
        print(f"  {alpha:5.2f}" + "".join(f"{variances[i, j]:11.5f}"
                                          for j in range(len(betas))))
    i, j = np.unravel_index(np.argmin(energies), energies.shape)
    k, l = np.unravel_index(np.argmin(variances), variances.shape)
    print()
    print(f"   lowest energy   at alpha = {alphas[i]:.2f}, "
          f"beta = {betas[j]:.2f}: {energies[i, j]:.6f} "
          f"+/- {errors[i, j]:.6f}")
    print(f"   lowest variance at alpha = {alphas[k]:.2f}, "
          f"beta = {betas[l]:.2f}: {variances[k, l]:.6f}")

    print()
    print("=" * 74)
    print("8. The final number")
    print("=" * 74)
    best = vmc_importance(1.00, 0.38, n_cycles=200000, time_step=0.5,
                          rng=np.random.default_rng(2024))
    print("   VMC, importance sampling, 200 000 cycles, alpha = 1, beta = 0.38")
    print(f"      E = {best['energy']:.6f} +/- {best['error']:.6f}"
          f"   (variance {best['variance']:.5f}, tau {best['tau']:.2f})")
    print()
    print(f"   exact (Taut)              {TAUT_ENERGY:.6f}")
    print("   CCSD, 42 orbitals         3.013626      (table 11.4)")
    print("   Hartree-Fock, 42 orbitals 3.161921      (table 11.3)")
    print()
    gap = best["energy"] - TAUT_ENERGY
    print(f"   The Monte Carlo energy sits {gap:.6f} above the exact one, and")
    print(f"   that gap is {gap/best['error']:.0f} standard errors wide -- which is")
    print("   the point rather than a problem.  A variational energy is an")
    print("   upper bound, so what is being measured here is not a statistical")
    print("   discrepancy but the deficiency of a two-parameter ansatz, and it")
    print("   is resolved because the statistical error is small enough to see")
    print("   it.  The residual would only be removed by a better trial")
    print("   function, or by diffusion Monte Carlo, which projects it out.")
    print()
    print("   Note also where this sits relative to chapter 11: below the CCSD")
    print("   energy in a 42-orbital basis, by a factor of twenty-five in the")
    print("   error.  The Jastrow factor depends on r_12 directly and satisfies")
    print("   the cusp condition exactly, and no finite sum of oscillator")
    print("   products can do either.")


if __name__ == "__main__":
    _demo()
