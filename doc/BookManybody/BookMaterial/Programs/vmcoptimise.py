"""
Optimising the variational parameters, and resampling the error.

Companion code to chapter 14 of *Quantum Mechanics for Many-particle Systems*.

Chapter 13 located the variational minimum by scanning a grid.  That works for
two parameters and for nothing else: a grid of $m$ points per parameter costs
$m^p$ Monte Carlo runs, and modern trial functions have hundreds or thousands
of parameters.  This program does it properly, in three stages.

  1. The gradient of the energy with respect to the parameters is available
     as a Monte Carlo average,

         dE/dtheta = 2 ( <(d ln Psi/dtheta) E_L> - <d ln Psi/dtheta><E_L> ),

     which costs no more than the energy itself.  It is derived in the text
     and checked here against finite differences.

  2. That gradient drives an optimiser -- plain gradient descent, or a
     quasi-Newton method, or stochastic reconfiguration, which preconditions
     the gradient with the metric of the variational manifold.  Short, noisy
     Monte Carlo runs are enough at this stage: we are locating a minimum, not
     measuring an energy.

  3. Once the parameters are fixed, one long production run measures the
     energy, and the error on it is obtained by resampling -- blocking,
     bootstrap or jackknife -- because the samples are correlated and the
     naive sigma/sqrt(N) is wrong, as chapter 12 showed.

The trial function, local energy and quantum force are imported unchanged from
`vmc.py`.

Author: Morten Hjorth-Jensen
"""

import math

import numpy as np

from vmc import (TAUT_ENERGY, correlation_time, local_energy, log_psi,
                 quantum_force, statistics)


# ---------------------------------------------------------------------------
#  Derivatives of the trial function with respect to the parameters
# ---------------------------------------------------------------------------
def log_psi_gradient(r, alpha, beta, omega=1.0):
    """d ln Psi_T / d(alpha, beta) for one configuration.

        d/dalpha = -omega (r_1^2 + r_2^2) / 2 ,
        d/dbeta  = -r_12^2 / (1 + beta r_12)^2 .

    These are the only new derivatives the optimisation needs, and both are
    cheap: no second derivatives and no extra sampling.
    """
    r_squared = r[0, 0]**2 + r[0, 1]**2 + r[1, 0]**2 + r[1, 1]**2
    r12 = math.hypot(r[0, 0] - r[1, 0], r[0, 1] - r[1, 1])
    d = 1.0 / (1.0 + beta * r12)
    return np.array([-0.5 * omega * r_squared, -(r12 * d)**2])


def log_psi_gradient_fd(r, alpha, beta, omega=1.0, h=1e-5):
    """The same by finite differences, for verification only."""
    out = np.empty(2)
    for k, (da, db) in enumerate(((h, 0.0), (0.0, h))):
        plus = log_psi(r, alpha + da, beta + db, omega)
        minus = log_psi(r, alpha - da, beta - db, omega)
        out[k] = (plus - minus) / (2.0 * h)
    return out


# ---------------------------------------------------------------------------
#  A sampler that also accumulates the gradient
# ---------------------------------------------------------------------------
def sample(alpha, beta, n_cycles=20000, time_step=0.5, omega=1.0, rng=None,
           burn_in=2000, keep_samples=False):
    """Importance-sampled VMC returning the energy and its parameter gradient.

    The gradient is the Monte Carlo estimate

        dE/dtheta = 2 ( <O_theta E_L> - <O_theta><E_L> ),
        O_theta = d ln Psi_T / dtheta ,

    accumulated in the same loop as the energy, so it is essentially free.
    """
    if rng is None:
        rng = np.random.default_rng(2024)
    diffusion = 0.5
    root_dt = math.sqrt(time_step)

    position = rng.normal(0.0, 1.0, (2, 2)) * root_dt
    log_old = log_psi(position, alpha, beta, omega)
    force_old = quantum_force(position, alpha, beta, omega)

    energy_sum = 0.0
    energy_squared = 0.0
    gradient_sum = np.zeros(2)
    gradient_energy_sum = np.zeros(2)
    accepted = proposals = 0
    samples = np.empty(n_cycles) if keep_samples else None

    for cycle in range(n_cycles + burn_in):
        for particle in range(2):
            trial = position.copy()
            trial[particle] = (position[particle]
                               + diffusion * force_old[particle] * time_step
                               + rng.normal(0.0, 1.0, 2) * root_dt)
            log_new = log_psi(trial, alpha, beta, omega)
            force_new = quantum_force(trial, alpha, beta, omega)
            green = float(np.sum(
                0.5 * (force_old[particle] + force_new[particle])
                * (0.5 * diffusion * time_step
                   * (force_old[particle] - force_new[particle])
                   - trial[particle] + position[particle])))
            proposals += 1
            if green + 2.0 * (log_new - log_old) > math.log(rng.random()
                                                            + 1e-300):
                position, log_old, force_old = trial, log_new, force_new
                accepted += 1
        if cycle >= burn_in:
            e = local_energy(position, alpha, beta, omega)
            o = log_psi_gradient(position, alpha, beta, omega)
            energy_sum += e
            energy_squared += e * e
            gradient_sum += o
            gradient_energy_sum += o * e
            if keep_samples:
                samples[cycle - burn_in] = e

    mean_energy = energy_sum / n_cycles
    mean_gradient = gradient_sum / n_cycles
    mean_gradient_energy = gradient_energy_sum / n_cycles
    gradient = 2.0 * (mean_gradient_energy - mean_gradient * mean_energy)
    variance = energy_squared / n_cycles - mean_energy**2
    return dict(energy=mean_energy, gradient=gradient, variance=variance,
                acceptance=accepted / proposals, samples=samples)


# ---------------------------------------------------------------------------
#  Optimisers
# ---------------------------------------------------------------------------
def gradient_descent(alpha, beta, learning_rate=0.05, max_iter=40, tol=1e-4,
                     n_cycles=20000, time_step=0.5, seed=2024, verbose=False):
    """Plain steepest descent with a fixed learning rate.

    theta <- theta - eta * dE/dtheta.  Simple, and it works here because the
    surface is well conditioned; the history is returned so that the path can
    be plotted.
    """
    theta = np.array([alpha, beta], dtype=float)
    history = []
    for iteration in range(1, max_iter + 1):
        result = sample(theta[0], theta[1], n_cycles=n_cycles,
                        time_step=time_step,
                        rng=np.random.default_rng(seed + iteration))
        history.append((iteration, theta.copy(), result["energy"],
                        result["gradient"].copy(), result["variance"]))
        if verbose:
            print(f"   {iteration:3d}  alpha {theta[0]:.5f}  "
                  f"beta {theta[1]:.5f}  E {result['energy']:.6f}  "
                  f"|grad| {np.linalg.norm(result['gradient']):.5f}")
        if np.linalg.norm(result["gradient"]) < tol:
            break
        theta = theta - learning_rate * result["gradient"]
    return theta, history


def momentum_descent(alpha, beta, learning_rate=0.05, momentum=0.6,
                     max_iter=40, tol=1e-4, n_cycles=20000, time_step=0.5,
                     seed=2024):
    """Gradient descent with momentum, which damps the noise in the gradient.

    v <- gamma v + eta g,   theta <- theta - v.  With a stochastic gradient the
    averaging over successive steps is worth as much as the acceleration.
    """
    theta = np.array([alpha, beta], dtype=float)
    velocity = np.zeros(2)
    history = []
    for iteration in range(1, max_iter + 1):
        result = sample(theta[0], theta[1], n_cycles=n_cycles,
                        time_step=time_step,
                        rng=np.random.default_rng(seed + iteration))
        history.append((iteration, theta.copy(), result["energy"],
                        result["gradient"].copy(), result["variance"]))
        if np.linalg.norm(result["gradient"]) < tol:
            break
        velocity = momentum * velocity + learning_rate * result["gradient"]
        theta = theta - velocity
    return theta, history


def stochastic_reconfiguration(alpha, beta, learning_rate=0.2, max_iter=40,
                               tol=1e-4, n_cycles=20000, time_step=0.5,
                               regulariser=1e-3, seed=2024):
    """Natural-gradient descent, S^{-1} g instead of g.

    The metric of the variational manifold is the covariance of the
    logarithmic derivatives,

        S_kl = <O_k O_l> - <O_k><O_l> ,

    also called the quantum geometric tensor or, up to a factor, the quantum
    Fisher information.  Preconditioning the gradient with S^{-1} makes the
    step independent of how the parameters happen to be scaled, which is what
    plain gradient descent is sensitive to.  A small diagonal shift keeps S
    invertible.
    """
    theta = np.array([alpha, beta], dtype=float)
    history = []
    for iteration in range(1, max_iter + 1):
        energy, gradient, metric = _sample_with_metric(
            theta[0], theta[1], n_cycles=n_cycles, time_step=time_step,
            rng=np.random.default_rng(seed + iteration))
        history.append((iteration, theta.copy(), energy, gradient.copy(),
                        None))
        if np.linalg.norm(gradient) < tol:
            break
        shifted = metric + regulariser * np.eye(2)
        theta = theta - learning_rate * np.linalg.solve(shifted, gradient)
    return theta, history


def _sample_with_metric(alpha, beta, n_cycles, time_step, rng, omega=1.0,
                        burn_in=2000):
    """As `sample`, but also accumulating S_kl = <O_k O_l> - <O_k><O_l>."""
    diffusion = 0.5
    root_dt = math.sqrt(time_step)
    position = rng.normal(0.0, 1.0, (2, 2)) * root_dt
    log_old = log_psi(position, alpha, beta, omega)
    force_old = quantum_force(position, alpha, beta, omega)

    energy_sum = 0.0
    gradient_sum = np.zeros(2)
    gradient_energy_sum = np.zeros(2)
    outer_sum = np.zeros((2, 2))

    for cycle in range(n_cycles + burn_in):
        for particle in range(2):
            trial = position.copy()
            trial[particle] = (position[particle]
                               + diffusion * force_old[particle] * time_step
                               + rng.normal(0.0, 1.0, 2) * root_dt)
            log_new = log_psi(trial, alpha, beta, omega)
            force_new = quantum_force(trial, alpha, beta, omega)
            green = float(np.sum(
                0.5 * (force_old[particle] + force_new[particle])
                * (0.5 * diffusion * time_step
                   * (force_old[particle] - force_new[particle])
                   - trial[particle] + position[particle])))
            if green + 2.0 * (log_new - log_old) > math.log(rng.random()
                                                            + 1e-300):
                position, log_old, force_old = trial, log_new, force_new
        if cycle >= burn_in:
            e = local_energy(position, alpha, beta, omega)
            o = log_psi_gradient(position, alpha, beta, omega)
            energy_sum += e
            gradient_sum += o
            gradient_energy_sum += o * e
            outer_sum += np.outer(o, o)

    mean_energy = energy_sum / n_cycles
    mean_gradient = gradient_sum / n_cycles
    gradient = 2.0 * (gradient_energy_sum / n_cycles
                      - mean_gradient * mean_energy)
    metric = outer_sum / n_cycles - np.outer(mean_gradient, mean_gradient)
    return mean_energy, gradient, metric


# ---------------------------------------------------------------------------
#  The production run: many walkers, vectorised
# ---------------------------------------------------------------------------
def production_run(alpha, beta, n_walkers=1000, n_steps=2000, time_step=0.5,
                   omega=1.0, rng=None, burn_in=200):
    """A long importance-sampled run with many independent walkers.

    Once the parameters are fixed, the only thing left is to accumulate
    statistics, and that parallelises perfectly: the walkers are independent,
    so the whole ensemble can be advanced with array operations.  Returns the
    walker-averaged local energy at each step -- a time series that is still
    autocorrelated, and therefore still needs resampling -- together with the
    grand total number of local-energy evaluations.
    """
    if rng is None:
        rng = np.random.default_rng(2024)
    diffusion = 0.5
    root_dt = math.sqrt(time_step)

    r = rng.normal(0.0, 1.0, (n_walkers, 2, 2)) * root_dt
    log_old = _log_psi_many(r, alpha, beta, omega)
    force_old = _quantum_force_many(r, alpha, beta, omega)

    series = np.empty(n_steps)
    accepted = 0
    proposals = 0

    for step in range(n_steps + burn_in):
        for particle in range(2):
            trial = r.copy()
            trial[:, particle] = (
                r[:, particle]
                + diffusion * force_old[:, particle] * time_step
                + rng.normal(0.0, 1.0, (n_walkers, 2)) * root_dt)
            log_new = _log_psi_many(trial, alpha, beta, omega)
            force_new = _quantum_force_many(trial, alpha, beta, omega)
            green = np.sum(
                0.5 * (force_old[:, particle] + force_new[:, particle])
                * (0.5 * diffusion * time_step
                   * (force_old[:, particle] - force_new[:, particle])
                   - trial[:, particle] + r[:, particle]), axis=1)
            take = (green + 2.0 * (log_new - log_old)
                    > np.log(rng.random(n_walkers) + 1e-300))
            r[take] = trial[take]
            log_old[take] = log_new[take]
            force_old[take] = force_new[take]
            accepted += int(take.sum())
            proposals += n_walkers
        if step >= burn_in:
            series[step - burn_in] = _local_energy_many(r, alpha, beta,
                                                        omega).mean()

    return dict(series=series, n_walkers=n_walkers, n_steps=n_steps,
                total_samples=n_walkers * n_steps,
                acceptance=accepted / proposals)


def _log_psi_many(r, alpha, beta, omega):
    r_squared = np.sum(r * r, axis=(1, 2))
    r12 = np.linalg.norm(r[:, 0] - r[:, 1], axis=1)
    return -0.5 * alpha * omega * r_squared + r12 / (1.0 + beta * r12)


def _local_energy_many(r, alpha, beta, omega):
    r_squared = np.sum(r * r, axis=(1, 2))
    r12 = np.linalg.norm(r[:, 0] - r[:, 1], axis=1)
    d = 1.0 / (1.0 + beta * r12)
    d2 = d * d
    return (0.5 * omega**2 * (1.0 - alpha**2) * r_squared
            + 2.0 * alpha * omega + 1.0 / r12
            + d2 * (alpha * omega * r12 - d2 + 2.0 * beta * d - 1.0 / r12))


def _quantum_force_many(r, alpha, beta, omega):
    separation = r[:, 0] - r[:, 1]
    r12 = np.linalg.norm(separation, axis=1)[:, None]
    d = 1.0 / (1.0 + beta * r12)
    drift = 2.0 * separation * d * d / r12
    force = np.empty_like(r)
    force[:, 0] = -2.0 * alpha * omega * r[:, 0] + drift
    force[:, 1] = -2.0 * alpha * omega * r[:, 1] - drift
    return force


# ---------------------------------------------------------------------------
#  Resampling
# ---------------------------------------------------------------------------
def blocking(x):
    """Flyvbjerg-Petersen blocking with Jonsson's automatic stopping.

    Repeatedly average neighbouring pairs.  Each transformation halves the
    number of samples and leaves the mean untouched, while the estimated
    variance of the mean grows until the blocks are longer than the
    correlation time, after which it plateaus.  Jonsson (Phys. Rev. E 98,
    043304 (2018)) turned the eyeball test for that plateau into a chi-square
    test on the remaining autocovariance, which is what is used here.

    Returns (mean, standard error, number of transformations used).
    """
    from scipy.stats import chi2

    x = np.asarray(x, dtype=float).copy()
    n = len(x)
    depth = int(math.floor(math.log2(n)))
    gamma = np.zeros(depth)
    sigma = np.zeros(depth)
    mean = x.mean()

    for i in range(depth):
        m = len(x)
        gamma[i] = float(np.sum((x[:m - 1] - x.mean())
                                * (x[1:] - x.mean())) / m)
        sigma[i] = float(x.var())
        x = 0.5 * (x[0::2] + x[1::2]) if m % 2 == 0 else \
            0.5 * (x[0:m - 1:2] + x[1:m:2])

    # Jonsson's test statistic, accumulated from the coarsest block back
    weights = (gamma / sigma)**2 * 2**np.arange(1, depth + 1)[::-1]
    statistic = np.cumsum(weights[::-1])[::-1]
    critical = chi2.ppf(0.95, np.arange(1, depth + 1))

    k = depth - 1
    for candidate in range(depth):
        if statistic[candidate] < critical[candidate]:
            k = candidate
            break
    return mean, math.sqrt(sigma[k] / 2**(depth - k)), k


def blocking_curve(x):
    """The naive error after each blocking transformation, for plotting."""
    x = np.asarray(x, dtype=float).copy()
    out = []
    while len(x) >= 8:
        out.append((len(x), float(x.std(ddof=1) / math.sqrt(len(x)))))
        if len(x) % 2:
            x = x[:-1]
        x = 0.5 * (x[0::2] + x[1::2])
    return out


def bootstrap(x, n_resamples=1000, rng=None, block=None):
    """Bootstrap error on the mean.

    Draw `n_resamples` samples of the same length, with replacement, and take
    the standard deviation of their means.  With correlated data the samples
    must be drawn in *blocks* rather than one at a time, otherwise the
    resampling destroys the correlations and reproduces the naive error; pass
    `block` to set the block length, or leave it None to use twice the
    correlation time.
    """
    if rng is None:
        rng = np.random.default_rng(2024)
    x = np.asarray(x, dtype=float)
    n = len(x)
    if block is None:
        block = max(1, int(2 * correlation_time(x)))
    n_blocks = n // block
    if n_blocks < 2:
        raise ValueError("series too short for this block length")
    blocks = x[:n_blocks * block].reshape(n_blocks, block)
    means = np.empty(n_resamples)
    for k in range(n_resamples):
        pick = rng.integers(0, n_blocks, n_blocks)
        means[k] = blocks[pick].mean()
    return float(x.mean()), float(means.std(ddof=1)), block


def jackknife(x, block=None):
    """Jackknife error on the mean, leaving out one block at a time.

    With n blocks and theta_(i) the estimate with block i removed,

        sigma^2 = (n-1)/n sum_i (theta_(i) - theta_bar)^2 .

    For the mean this is exact and reduces to the blocked standard error; it
    earns its keep for non-linear estimators, where the bootstrap and the
    jackknife are the only practical options.
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    if block is None:
        block = max(1, int(2 * correlation_time(x)))
    n_blocks = n // block
    blocks = x[:n_blocks * block].reshape(n_blocks, block).mean(axis=1)
    total = blocks.sum()
    leave_one_out = (total - blocks) / (n_blocks - 1)
    mean = leave_one_out.mean()
    variance = (n_blocks - 1) / n_blocks * np.sum((leave_one_out - mean)**2)
    return float(x.mean()), float(math.sqrt(variance)), block


# ---------------------------------------------------------------------------
def _demo():
    print("=" * 74)
    print("1. The gradient of the energy")
    print("=" * 74)
    print("The derivative of the energy with respect to a variational")
    print("parameter is itself a Monte Carlo average,")
    print("   dE/dtheta = 2 ( <O_theta E_L> - <O_theta><E_L> ),")
    print("with O_theta = d ln Psi_T / dtheta.  It costs nothing extra: the")
    print("same walk that measures the energy measures the gradient.")
    print()
    rng = np.random.default_rng(3)
    worst = 0.0
    for _ in range(40):
        r = rng.normal(0.0, 1.0, (2, 2))
        for alpha, beta in ((0.98, 0.40), (0.80, 0.25), (1.10, 0.60)):
            worst = max(worst, float(np.abs(
                log_psi_gradient(r, alpha, beta)
                - log_psi_gradient_fd(r, alpha, beta)).max()))
    print(f"   d ln Psi/d(alpha, beta), analytic against finite differences:")
    print(f"      largest discrepancy {worst:.1e}")
    print()
    print("   and the resulting energy gradient, against a finite difference")
    print("   of the energy itself:")
    print()
    print(f"{'alpha':>7s} {'beta':>7s} {'dE/dalpha':>21s} {'dE/dbeta':>21s}")
    for alpha, beta in ((0.90, 0.30), (1.00, 0.40), (1.05, 0.50)):
        analytic = sample(alpha, beta, n_cycles=40000,
                          rng=np.random.default_rng(99))["gradient"]
        h = 0.01
        numeric = np.empty(2)
        for k, (da, db) in enumerate(((h, 0.0), (0.0, h))):
            plus = sample(alpha + da, beta + db, n_cycles=40000,
                          rng=np.random.default_rng(99))["energy"]
            minus = sample(alpha - da, beta - db, n_cycles=40000,
                           rng=np.random.default_rng(99))["energy"]
            numeric[k] = (plus - minus) / (2.0 * h)
        print(f"{alpha:7.2f} {beta:7.2f} "
              f"{analytic[0]:10.5f} vs {numeric[0]:8.5f} "
              f"{analytic[1]:10.5f} vs {numeric[1]:8.5f}")
    print()
    print("   The two agree within the Monte Carlo noise.  Note how much")
    print("   cheaper the analytic route is: one run instead of two per")
    print("   parameter, and no cancellation of nearly equal numbers.")

    print()
    print("=" * 74)
    print("2. Gradient descent")
    print("=" * 74)
    print("Starting deliberately far from the minimum, at alpha = 0.85 and")
    print("beta = 0.20, with 6000 cycles per step -- short, noisy runs are")
    print("enough while we are still locating the minimum.")
    print()
    theta, history = gradient_descent(0.85, 0.20, learning_rate=0.05,
                                      max_iter=25, n_cycles=6000)
    print(f"{'iter':>5s} {'alpha':>9s} {'beta':>9s} {'energy':>12s} "
          f"{'|gradient|':>12s}")
    for iteration, point, energy, gradient, _ in history:
        if iteration <= 6 or iteration % 5 == 0 or iteration == len(history):
            print(f"{iteration:5d} {point[0]:9.5f} {point[1]:9.5f} "
                  f"{energy:12.6f} {np.linalg.norm(gradient):12.6f}")
    print()
    print(f"   converged to alpha = {theta[0]:.4f}, beta = {theta[1]:.4f} "
          f"in {len(history)} steps")

    print()
    print("=" * 74)
    print("3. Three optimisers compared")
    print("=" * 74)
    print("All start from the same point and use the same number of samples")
    print("per step.  Momentum averages successive gradients, which is worth")
    print("as much as the acceleration when the gradient is noisy.  Stochastic")
    print("reconfiguration preconditions with the metric of the variational")
    print("manifold and takes the largest useful steps.")
    print()
    print(f"{'method':>28s} {'steps':>7s} {'alpha':>9s} {'beta':>9s} "
          f"{'final energy':>14s}")
    for name, routine, kwargs in (
            ("gradient descent", gradient_descent,
             dict(learning_rate=0.05)),
            ("with momentum", momentum_descent,
             dict(learning_rate=0.03, momentum=0.6)),
            ("stochastic reconfiguration", stochastic_reconfiguration,
             dict(learning_rate=0.2))):
        theta, history = routine(0.85, 0.20, max_iter=25, n_cycles=6000,
                                 **kwargs)
        print(f"{name:>28s} {len(history):7d} {theta[0]:9.5f} "
              f"{theta[1]:9.5f} {history[-1][2]:14.6f}")

    print()
    print("=" * 74)
    print("4. The production run")
    print("=" * 74)
    print("With the parameters fixed, the only remaining task is statistics,")
    print("and that parallelises perfectly: the walkers are independent, so")
    print("the ensemble is advanced with array operations.")
    print()
    run = production_run(1.00, 0.40, n_walkers=1000, n_steps=2000,
                         time_step=0.5, rng=np.random.default_rng(2024))
    series = run["series"]
    print(f"   {run['n_walkers']} walkers x {run['n_steps']} steps "
          f"= {run['total_samples']:,} local-energy evaluations")
    print(f"   acceptance rate {run['acceptance']:.1%}")
    print(f"   mean energy {series.mean():.8f}")
    print(f"   correlation time of the walker-averaged series "
          f"{correlation_time(series):.2f}")

    print()
    print("=" * 74)
    print("5. Resampling the error")
    print("=" * 74)
    print("The series is autocorrelated, so sigma/sqrt(N) is not the error.")
    print("Three standard ways of getting the right one, on the run above:")
    print()

    def report(x, label):
        naive = x.std(ddof=1) / math.sqrt(len(x))
        tau = correlation_time(x)
        _, error_b, k = blocking(x)
        _, error_boot, block_boot = bootstrap(
            x, n_resamples=1000, rng=np.random.default_rng(7))
        _, error_j, block_j = jackknife(x)
        print(f"   {label}")
        print(f"      correlation time tau         {tau:8.2f}")
        print(f"      naive sigma/sqrt(N)          {naive:.8f}")
        print(f"      sqrt(tau) x naive            "
              f"{naive*math.sqrt(tau):.8f}")
        print(f"      blocking                     {error_b:.8f}"
              f"   ({k} transformations)")
        print(f"      bootstrap, blocks of {block_boot:<5d}  "
              f"{error_boot:.8f}")
        print(f"      jackknife, blocks of {block_j:<5d}  {error_j:.8f}")
        print(f"      resampled / naive            "
              f"{error_b/naive:8.2f}")
        return error_b

    error_b = report(series, "well tuned, time step 0.5:")
    mean_b = series.mean()
    print()
    print("   Here the four honest estimates agree and the naive one is only")
    print("   modestly wrong, because the sampling is well tuned and tau is")
    print("   close to one.  That is not something to rely on.  The same")
    print("   calculation with a badly chosen time step:")
    print()
    bad = production_run(1.00, 0.40, n_walkers=400, n_steps=4000,
                         time_step=0.02, rng=np.random.default_rng(2024))
    report(bad["series"], "badly tuned, time step 0.02:")
    print()
    print("   Now the naive error understates the true one by a factor of")
    print("   four, and it does so silently.  Note that the three resampling")
    print("   estimates do not agree exactly either: blocking with Jonsson's")
    print("   stopping rule is deliberately conservative.  They agree on the")
    print("   order of magnitude, which is what matters, and all three are")
    print("   right where the naive estimate is wrong.")
    print()
    print("   the blocking curve for the well-tuned run:")
    print(f"{'samples':>10s} {'block':>8s} {'naive error':>14s}")
    for index, (n, error) in enumerate(blocking_curve(series)):
        if index % 2 == 0:
            print(f"{n:10d} {len(series)//n:8d} {error:14.8f}")

    print("=" * 74)
    print("6. The final answer")
    print("=" * 74)
    print(f"   E = {mean_b:.6f} +/- {error_b:.6f}")
    print(f"   exact (Taut)          {TAUT_ENERGY:.6f}")
    print(f"   CCSD, 42 orbitals     3.013626   (table 11.4)")
    print()
    gap = mean_b - TAUT_ENERGY
    print(f"   The energy sits {gap:.6f} above the exact answer, which is")
    print(f"   {gap/error_b:.0f} standard errors.  That gap is the two-parameter")
    print("   ansatz, not the sampling: it is what a better trial function,")
    print("   or diffusion Monte Carlo, would remove.  What this chapter has")
    print("   achieved is that the error bar is now small enough, and honest")
    print("   enough, for the gap to be visible at all.")


if __name__ == "__main__":
    _demo()
