"""
Diffusion Monte Carlo for two electrons in a two-dimensional quantum dot.

Companion code to chapter 15 of *Quantum Mechanics for Many-particle Systems*.

Chapters 13 and 14 pushed a two-parameter trial function as far as it goes:
optimised properly, with a trustworthy error bar, it gives 3.000549(49)
against Taut's exact 3.  The residual is not statistical.  It is what the
ansatz cannot represent, and no amount of sampling removes it.

Diffusion Monte Carlo removes it by projecting instead of guessing.  The
imaginary-time Schroedinger equation

    -d Phi / d tau = ( H - E_T ) Phi

has the formal solution Phi(tau) = exp[-(H - E_T) tau] Phi(0), and expanding
Phi(0) in eigenstates shows that every excited component decays relative to the
ground state like exp[-(E_n - E_0) tau].  Whatever we start from, we end at
Phi_0 -- provided the starting state is not orthogonal to it.

The equation for Phi alone is a diffusion equation with a position-dependent
source term, and the source is the potential, which diverges as 1/r_12.  That
makes the naive algorithm useless.  Importance sampling fixes it: propagate

    f(R, tau) = Psi_T(R) Phi(R, tau)

instead, where Psi_T is the optimised trial function of chapter 14.  Then

    -d f / d tau = -D grad^2 f + D grad . ( f F ) + ( E_L(R) - E_T ) f ,

with D = 1/2, the quantum force F = 2 grad Psi_T / Psi_T of chapter 13, and the
local energy E_L in the source term rather than the potential.  Three things
have been bought at once: the drift F steers walkers into the region where
Psi_T is large, the source term is bounded because the Jastrow cusp makes E_L
finite, and the ground-state energy comes out of the *mixed* estimator

    E_0 = < E_L >_f ,

which is the same average of the same local energy as in variational Monte
Carlo -- only over a different distribution.

The algorithm is therefore the importance-sampled Metropolis-Hastings step of
chapter 13, followed by a birth-death step with weight

    w = exp[ -dt ( (E_L(R) + E_L(R')) / 2 - E_T ) ] .

Everything else -- trial function, local energy, quantum force -- is imported
unchanged from `vmc.py`, and the parameters come from `vmcoptimise.py`.

Author: Morten Hjorth-Jensen
"""

import math

import numpy as np

from vmc import TAUT_ENERGY, correlation_time, local_energy, log_psi, \
    quantum_force
from vmcoptimise import (_local_energy_many, _log_psi_many,
                         _quantum_force_many, blocking, bootstrap)

#  The variational optimum of chapter 14.
ALPHA_OPT = 1.00
BETA_OPT = 0.40
E_VMC = 3.000549


# ---------------------------------------------------------------------------
#  A drift-diffusion Metropolis-Hastings sweep over an ensemble of walkers
# ---------------------------------------------------------------------------
def _sweep(r, log_old, force_old, alpha, beta, omega, time_step, rng):
    """Move every particle of every walker once, with accept/reject.

    This is exactly the importance-sampling step of chapter 13, applied to a
    whole ensemble at once.  It returns the updated arrays together with the
    numbers needed for the effective time step: the sum of the squared
    displacements attempted and the sum of those actually accepted.
    """
    n_walkers = r.shape[0]
    diffusion = 0.5
    root_dt = math.sqrt(time_step)
    attempted = accepted_distance = 0.0
    accepted = 0

    for particle in range(2):
        trial = r.copy()
        drift = diffusion * force_old[:, particle] * time_step
        noise = rng.normal(0.0, 1.0, (n_walkers, 2)) * root_dt
        trial[:, particle] = r[:, particle] + drift + noise
        log_new = _log_psi_many(trial, alpha, beta, omega)
        force_new = _quantum_force_many(trial, alpha, beta, omega)

        green = np.sum(
            0.5 * (force_old[:, particle] + force_new[:, particle])
            * (0.5 * diffusion * time_step
               * (force_old[:, particle] - force_new[:, particle])
               - trial[:, particle] + r[:, particle]), axis=1)
        take = (green + 2.0 * (log_new - log_old)
                > np.log(rng.random(n_walkers) + 1e-300))

        displacement = np.sum((trial[:, particle] - r[:, particle])**2, axis=1)
        attempted += float(displacement.sum())
        accepted_distance += float(displacement[take].sum())
        accepted += int(take.sum())

        r[take] = trial[take]
        log_old[take] = log_new[take]
        force_old[take] = force_new[take]

    return r, log_old, force_old, attempted, accepted_distance, \
        accepted / (2 * n_walkers)


# ---------------------------------------------------------------------------
#  Equilibration: start the walkers from the variational distribution
# ---------------------------------------------------------------------------
def vmc_ensemble(alpha, beta, n_walkers=500, n_steps=300, time_step=0.1,
                 omega=1.0, rng=None):
    """An ensemble distributed as |Psi_T|^2, for use as the DMC start.

    Starting a diffusion Monte Carlo run from the variational distribution is
    not merely convenient: |Psi_T|^2 is already close to f = Psi_T Phi_0 when
    Psi_T is good, so the projection has little to do and the equilibration
    phase is short.
    """
    if rng is None:
        rng = np.random.default_rng(2024)
    r = rng.normal(0.0, 1.0, (n_walkers, 2, 2))
    log_old = _log_psi_many(r, alpha, beta, omega)
    force_old = _quantum_force_many(r, alpha, beta, omega)
    for _ in range(n_steps):
        r, log_old, force_old, _, _, _ = _sweep(
            r, log_old, force_old, alpha, beta, omega, time_step, rng)
    return r


# ---------------------------------------------------------------------------
#  The diffusion Monte Carlo algorithm
# ---------------------------------------------------------------------------
def dmc(alpha=ALPHA_OPT, beta=BETA_OPT, n_walkers=500, n_steps=4000,
        time_step=0.01, omega=1.0, rng=None, burn_in=1000,
        feedback=0.05, effective_time_step=True, max_walkers=None,
        start=None):
    """Importance-sampled diffusion Monte Carlo.

    One step consists of

      1. a drift-diffusion move of every particle, accepted or rejected with
         the Metropolis-Hastings ratio of chapter 13.  The accept/reject step
         is not part of the original algorithm and it is not optional here:
         without it the walkers sample the Green's function only to O(dt) and
         the time-step error is an order of magnitude larger;

      2. a branching step, in which the walker is replaced by

             m = floor( w + u ),   u ~ U(0, 1),

         copies of itself, with the branching weight

             w = exp[ -dt_eff ( (E_L(R) + E_L(R')) / 2 - E_T ) ] .

         The average of E_L over the endpoints is the trapezoidal rule for the
         integral of the source term along the step, and it is one order more
         accurate than using either endpoint alone;

      3. an update of the trial energy E_T, which is the only thing holding the
         population at its target value.

    `effective_time_step` replaces dt in the branching factor by

        dt_eff = dt * ( sum of accepted squared displacements )
                    / ( sum of attempted squared displacements ) ,

    which corrects the branching for the moves that the accept/reject step
    threw away.  Rejected moves mean the walker diffused less than dt would
    suggest, and it should therefore be reweighted by less.  Umrigar,
    Nightingale and Runge (J. Chem. Phys. 99, 2865 (1993)) introduced it, and
    it removes most of the remaining time-step error, as the demo shows.

    Returns the mixed-estimator energy and its blocking error, the growth
    estimator, the population history and the acceptance rate.
    """
    if rng is None:
        rng = np.random.default_rng(2024)
    if max_walkers is None:
        max_walkers = 10 * n_walkers

    r = vmc_ensemble(alpha, beta, n_walkers, 300, max(time_step, 0.05),
                     omega, rng) if start is None else start.copy()
    target = r.shape[0]
    log_old = _log_psi_many(r, alpha, beta, omega)
    force_old = _quantum_force_many(r, alpha, beta, omega)
    energy_old = _local_energy_many(r, alpha, beta, omega)

    trial_energy = float(energy_old.mean())
    running = trial_energy
    series = np.empty(n_steps)
    growth = np.empty(n_steps)
    spread = np.empty(n_steps)
    population = np.empty(n_steps + burn_in, dtype=int)
    acceptance = 0.0
    dt_eff_sum = 0.0

    for step in range(n_steps + burn_in):
        r, log_new, force_new, attempted, accepted_distance, rate = _sweep(
            r, log_old.copy(), force_old.copy(), alpha, beta, omega,
            time_step, rng)
        acceptance += rate
        energy_new = _local_energy_many(r, alpha, beta, omega)

        dt_branch = time_step
        if effective_time_step and attempted > 0.0:
            dt_branch = time_step * accepted_distance / attempted
        dt_eff_sum += dt_branch

        # -- birth and death ------------------------------------------------
        weight = np.exp(-dt_branch * (0.5 * (energy_old + energy_new)
                                      - trial_energy))
        copies = np.floor(weight + rng.random(r.shape[0])).astype(int)
        copies = np.clip(copies, 0, 3)
        if copies.sum() == 0:                      # pathological, restart
            copies[rng.integers(len(copies))] = 1

        r = np.repeat(r, copies, axis=0)
        log_new = np.repeat(log_new, copies)
        force_new = np.repeat(force_new, copies, axis=0)
        energy_new = np.repeat(energy_new, copies)

        if len(r) > max_walkers:                   # hard ceiling, rarely hit
            keep = rng.choice(len(r), max_walkers, replace=False)
            r, log_new = r[keep], log_new[keep]
            force_new, energy_new = force_new[keep], energy_new[keep]

        r, log_old, force_old, energy_old = r.copy(), log_new, force_new, \
            energy_new
        population[step] = len(r)

        # -- population control ---------------------------------------------
        #  E_T is pushed down when the population grows and up when it
        #  shrinks.  The feedback must be weak: a strong one imposes a
        #  fluctuation on E_T that feeds straight back into the branching.
        mixed = float(energy_old.mean())
        running += 0.01 * (mixed - running)
        trial_energy = running - (feedback / time_step) * math.log(
            len(r) / target)

        if step >= burn_in:
            series[step - burn_in] = mixed
            growth[step - burn_in] = trial_energy
            spread[step - burn_in] = float(energy_old.var(ddof=1))

    energy, error, _ = blocking(series)
    return dict(energy=energy, error=error,
                growth=float(growth.mean()),
                growth_error=blocking(growth)[1],
                variance=float(spread.mean()),
                tau=correlation_time(series),
                acceptance=acceptance / (n_steps + burn_in),
                time_step=time_step,
                effective_time_step=dt_eff_sum / (n_steps + burn_in),
                population=population, series=series,
                mean_population=float(population[burn_in:].mean()),
                alpha=alpha, beta=beta)


# ---------------------------------------------------------------------------
#  Time-step extrapolation
# ---------------------------------------------------------------------------
def time_step_scan(time_steps, n_walkers=500, n_steps=4000, seed=2024,
                   **kwargs):
    """Run DMC at a sequence of time steps, each with its own generator."""
    results = []
    for time_step in time_steps:
        results.append(dmc(n_walkers=n_walkers, n_steps=n_steps,
                           time_step=time_step,
                           rng=np.random.default_rng(seed), **kwargs))
    return results


def extrapolate(time_steps, energies, errors=None, degree=1):
    """Weighted least-squares fit E(dt) = E(0) + c dt (+ c2 dt^2).

    Returns the intercept, its standard error, and the fitted coefficients.
    The intercept is the number to quote: it is the dt -> 0 limit, which is
    where the short-time approximation to the Green's function becomes exact.
    """
    time_steps = np.asarray(time_steps, dtype=float)
    energies = np.asarray(energies, dtype=float)
    design = np.vander(time_steps, degree + 1, increasing=True)
    if errors is None:
        weights = np.ones_like(energies)
    else:
        weights = 1.0 / np.asarray(errors, dtype=float)
    weighted_design = design * weights[:, None]
    coefficients, *_ = np.linalg.lstsq(weighted_design, energies * weights,
                                       rcond=None)
    covariance = np.linalg.inv(weighted_design.T @ weighted_design)
    residual = energies - design @ coefficients
    chi_square = float(np.sum((residual * weights)**2))
    dof = max(len(energies) - degree - 1, 1)
    scale = max(chi_square / dof, 1.0)
    return dict(intercept=float(coefficients[0]),
                intercept_error=float(math.sqrt(covariance[0, 0] * scale)),
                coefficients=coefficients,
                chi_square=chi_square, dof=dof)


# ---------------------------------------------------------------------------
#  A one-dimensional illustration: DMC without importance sampling
# ---------------------------------------------------------------------------
def naive_dmc_harmonic(n_walkers=2000, n_steps=4000, time_step=0.01,
                       burn_in=1000, rng=None):
    """Pure branching random walk for the one-dimensional oscillator.

    No trial function, no drift, no accept/reject: walkers diffuse freely and
    are branched with weight exp[-dt (V(x) - E_T)], V = x^2/2.  The exact
    ground-state energy is 1/2.  This is the algorithm in its original form,
    and it works here only because the potential is bounded below and smooth.
    For the Coulomb problem the same code would be useless, which is the whole
    argument for importance sampling.
    """
    if rng is None:
        rng = np.random.default_rng(2024)
    x = rng.normal(0.0, 1.0, n_walkers)
    target = n_walkers
    trial_energy = 0.5 * float((x * x).mean())
    running = trial_energy
    series = np.empty(n_steps)
    root_dt = math.sqrt(time_step)

    for step in range(n_steps + burn_in):
        potential_old = 0.5 * x * x
        x = x + rng.normal(0.0, 1.0, len(x)) * root_dt
        potential_new = 0.5 * x * x
        weight = np.exp(-time_step * (0.5 * (potential_old + potential_new)
                                      - trial_energy))
        copies = np.clip(np.floor(weight + rng.random(len(x))), 0,
                         3).astype(int)
        x = np.repeat(x, copies)
        if len(x) == 0:
            x = rng.normal(0.0, 1.0, target)
        potential = 0.5 * float((x * x).mean())
        running += 0.01 * (potential - running)
        trial_energy = running - 0.05 / time_step * math.log(len(x) / target)
        if step >= burn_in:
            series[step - burn_in] = trial_energy
    mean, error, _ = blocking(series)
    return dict(energy=mean, error=error, exact=0.5)


# ---------------------------------------------------------------------------
def _demo():
    print("=" * 74)
    print("1. Diffusion Monte Carlo in its original form")
    print("=" * 74)
    print("Free diffusion plus branching on exp[-dt(V - E_T)], no trial")
    print("function at all, for the one-dimensional harmonic oscillator.")
    print("The energy comes from the growth of the population alone.")
    print()
    naive = naive_dmc_harmonic(n_steps=3000)
    print(f"   E = {naive['energy']:.6f} +/- {naive['error']:.6f}"
          f"   (exact {naive['exact']:.6f})")
    print()
    print("   That works because x^2/2 is smooth and bounded below.  Put")
    print("   1/r_12 in the exponent instead and the weights diverge as two")
    print("   walkers approach: the population is then dominated by a handful")
    print("   of configurations and the variance is unbounded.  Importance")
    print("   sampling replaces V by E_L, which the cusp condition keeps")
    print("   finite, and that is what makes the method usable.")

    print()
    print("=" * 74)
    print("2. The quantum dot, from the variational starting point")
    print("=" * 74)
    print(f"   trial function: alpha = {ALPHA_OPT:.2f}, beta = {BETA_OPT:.2f}"
          "  (chapter 14)")
    print(f"   VMC energy    : {E_VMC:.6f}")
    print(f"   exact (Taut)  : {TAUT_ENERGY:.6f}")
    print()
    result = dmc(n_walkers=500, n_steps=4000, time_step=0.05,
                 rng=np.random.default_rng(2024))
    print(f"   DMC, dt = {result['time_step']:.3f}, "
          f"{result['mean_population']:.0f} walkers, "
          f"{result['acceptance']:.1%} acceptance")
    print(f"      mixed estimator  E = {result['energy']:.6f} "
          f"+/- {result['error']:.6f}")
    print(f"      growth estimator E = {result['growth']:.6f} "
          f"+/- {result['growth_error']:.6f}")
    print(f"      variance {result['variance']:.5f}, tau {result['tau']:.2f}")
    print()
    print("   The two agree, as they must: the mixed estimator measures the")
    print("   eigenvalue directly and the growth estimator measures the rate")
    print("   at which the population would grow if E_T were held fixed.  The")
    print("   mixed one has by far the smaller variance, which is why it is")
    print("   the one quoted.")

    print()
    print("=" * 74)
    print("3. The time step is the systematic error")
    print("=" * 74)
    print("The short-time Green's function is exact only as dt -> 0.  Running")
    print("at several time steps and extrapolating is not optional -- it is")
    print("how the answer is obtained.")
    print()
    time_steps = [0.5, 0.3, 0.2, 0.1, 0.05, 0.02]
    for label, use_eff in (("plain dt", False), ("with dt_eff", True)):
        print(f"   {label}")
        print(f"   {'dt':>8s} {'energy':>12s} {'error':>10s} "
              f"{'accept':>8s} {'dt_eff/dt':>10s} {'tau':>7s}")
        energies, errors = [], []
        for time_step in time_steps:
            out = dmc(n_walkers=500, n_steps=4000, time_step=time_step,
                      effective_time_step=use_eff,
                      rng=np.random.default_rng(2024))
            energies.append(out["energy"])
            errors.append(out["error"])
            print(f"   {time_step:8.4f} {out['energy']:12.6f} "
                  f"{out['error']:10.6f} {out['acceptance']:8.1%} "
                  f"{out['effective_time_step']/time_step:10.3f} "
                  f"{out['tau']:7.1f}")
        fit = extrapolate(time_steps, energies, errors)
        print(f"   extrapolated to dt = 0: {fit['intercept']:.6f} "
              f"+/- {fit['intercept_error']:.6f}   "
              f"(slope {fit['coefficients'][1]:+.5f})")
        print()
    print("   The bias is minute -- a few times 1e-4 even at dt = 0.5 -- and")
    print("   that is a statement about the trial function, not about the")
    print("   method.  The time-step error comes from the commutators that")
    print("   the short-time factorisation drops, and those vanish when Psi_T")
    print("   is an eigenfunction.  At alpha = 1 the local energy of this dot")
    print("   is very nearly constant, so there is almost nothing for the")
    print("   factorisation to get wrong.  Note also that the correlation")
    print("   time in *steps* grows as dt shrinks: what decorrelates a walker")
    print("   is imaginary time, so halving dt doubles the number of steps")
    print("   needed for the same statistics.")

    print()
    print("=" * 74)
    print("4. The same calculation from a deliberately poor guiding function")
    print("=" * 74)
    print("alpha = 0.80, beta = 0.10 is far off the variational minimum.  In")
    print("chapter 13 this would simply be a worse answer.  Here it is not:")
    print("the projection removes the error in the wave function, and only")
    print("the efficiency and the time-step sensitivity suffer.")
    print()
    poor = [0.1, 0.05, 0.02, 0.01]
    print(f"   {'dt':>8s} {'energy':>12s} {'error':>10s} {'tau':>7s}")
    energies, errors = [], []
    for time_step in poor:
        out = dmc(alpha=0.80, beta=0.10, n_walkers=2000, n_steps=3000,
                  time_step=time_step, rng=np.random.default_rng(2024))
        energies.append(out["energy"])
        errors.append(out["error"])
        print(f"   {time_step:8.4f} {out['energy']:12.6f} "
              f"{out['error']:10.6f} {out['tau']:7.1f}")
    fit = extrapolate(poor, energies, errors)
    print(f"   extrapolated to dt = 0: {fit['intercept']:.6f} "
          f"+/- {fit['intercept_error']:.6f}   "
          f"(slope {fit['coefficients'][1]:+.5f})")
    print()
    print("   Compare with the variational energy of the same trial function,")
    print("   which is 3.363877(5019).  Diffusion Monte Carlo turns a bad")
    print("   ansatz into the right answer; it charges for it in error bar --")
    print("   here some twenty-five times the error of the good guiding")
    print("   function, for four times the walkers -- and in a time-step")
    print("   dependence about five hundred times stronger, so that the")
    print("   largest usable dt is smaller by the same measure.")

    print()
    print("=" * 74)
    print("5. It is not a population-control bias")
    print("=" * 74)
    print("Worth checking, because a residual that does not move with the")
    print("time step often turns out to move with the number of walkers.")
    print()
    print(f"   {'walkers':>9s} {'good guide':>22s} {'poor guide':>22s}")
    for n_walkers in (250, 500, 1000, 2000):
        good = dmc(n_walkers=n_walkers, n_steps=3000, time_step=0.1,
                   rng=np.random.default_rng(2024))
        bad = dmc(alpha=0.80, beta=0.10, n_walkers=n_walkers, n_steps=3000,
                  time_step=0.1, rng=np.random.default_rng(2024))
        print(f"   {n_walkers:9d} "
              f"{good['energy']:14.6f}({good['error']*1e6:5.0f}) "
              f"{bad['energy']:14.6f}({bad['error']*1e6:5.0f})")
    print()
    print("   Neither column drifts: the error bars shrink as 1/sqrt(N) and")
    print("   the central values do not move.  The residual in the poor-guide")
    print("   column at dt = 0.1 is therefore time-step error, and section 3")
    print("   confirms it by extrapolating it away.")

    print()
    print("=" * 74)
    print("6. Where this leaves the quantum dot")
    print("=" * 74)
    print("   Hartree-Fock, 42 orbitals   3.161921     (table 11.3)")
    print("   MP2, 42 orbitals            3.027038     (table 11.4)")
    print("   CCSD, 42 orbitals           3.013626     (table 11.4)")
    print(f"   VMC, optimised              {E_VMC:.6f}     (chapter 14)")
    print("   DMC, extrapolated           see section 3")
    print(f"   exact (Taut)                {TAUT_ENERGY:.6f}")
    print()
    print("   The two-electron singlet ground state is nodeless, so there is")
    print("   no fixed-node error here and the projection is exact in")
    print("   principle.  What remains is the time step, the finite")
    print("   population, and the statistics -- all three controllable, none")
    print("   of them a property of the trial function.  That is the whole")
    print("   difference from chapters 13 and 14.")


if __name__ == "__main__":
    _demo()
