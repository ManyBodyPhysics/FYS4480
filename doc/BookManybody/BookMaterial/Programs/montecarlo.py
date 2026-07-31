"""
Elements from statistics: random numbers, Markov chains and the Metropolis
algorithm.

Companion code to chapter 12 of *Quantum mechanics for Many-particle Systems*.

The chapters up to now expanded the wave function in a basis and the error was
a truncation.  From here on the error is statistical instead, and controlling
it needs a different set of tools: probability distributions and their
moments, the central limit theorem, the correlation time that tells us how
many of our samples were really independent, and the Metropolis algorithm that
generates samples from a distribution we can only evaluate up to a constant.

This program works through all of them, and ends with a first variational
Monte Carlo calculation on the two-electron quantum dot of chapter 11.  That
last section is the point of the chapter: the trial function with alpha = 1 is
the lowest oscillator orbital, so its energy must come out at the
minimal-basis Hartree-Fock value 3.25331414 of table 11.1, and it does, to
within the statistical error.

Author: Morten Hjorth-Jensen
"""

import math

import numpy as np


# ---------------------------------------------------------------------------
#  Random number generators
# ---------------------------------------------------------------------------
class LinearCongruential:
    """The textbook generator x_{n+1} = (a x_n + c) mod M.

    Included to be looked at rather than used.  The parameters below are those
    of the notorious RANDU, which passes every one-dimensional test and fails
    catastrophically in three dimensions: its triples lie on fifteen planes.
    That failure mode -- perfectly uniform marginals, badly correlated
    tuples -- is exactly the one a Monte Carlo integration is sensitive to.
    """

    def __init__(self, seed=1, a=65539, c=0, modulus=2**31):
        self.state = seed
        self.a, self.c, self.modulus = a, c, modulus

    def __call__(self):
        self.state = (self.a * self.state + self.c) % self.modulus
        return self.state / self.modulus

    def sample(self, n):
        return np.array([self() for _ in range(n)])


def uniformity_test(x, bins=10):
    """Chi-square-like measure of how evenly x fills [0, 1)."""
    counts, _ = np.histogram(x, bins=bins, range=(0.0, 1.0))
    expected = len(x) / bins
    return float(np.sum((counts - expected) ** 2 / expected))


def serial_correlation(x, lag=1):
    """Correlation between x_i and x_{i+lag}; should be ~0 for a good RNG."""
    a, b = x[:-lag], x[lag:]
    return float(np.corrcoef(a, b)[0, 1])


def planes_test(generator, n=3000):
    """Fraction of consecutive triples lying on the RANDU planes.

    For RANDU, x_{n+2} = 6 x_{n+1} - 9 x_n (mod 1) exactly, so the quantity
    below is an integer for every triple.  A good generator gives no such
    relation.
    """
    x = generator.sample(n)
    residual = (x[2:] - 6.0 * x[1:-1] + 9.0 * x[:-2]) % 1.0
    return float(np.mean(np.minimum(residual, 1.0 - residual) < 1e-12))


# ---------------------------------------------------------------------------
#  Monte Carlo integration
# ---------------------------------------------------------------------------
def brute_force(f, n, rng, a=0.0, b=1.0):
    """Integral of f over [a, b] by uniform sampling, with its error.

    Returns (estimate, error).  The error is sigma/sqrt(n) with sigma the
    sample standard deviation -- valid only because the samples are
    independent, which is exactly what stops being true once we use a Markov
    chain.
    """
    x = a + (b - a) * rng.random(n)
    values = (b - a) * f(x)
    return float(values.mean()), float(values.std(ddof=1) / math.sqrt(n))


def importance_sampled(f, p, sampler, n, rng):
    """Integral of f by sampling from p and averaging f/p."""
    x = sampler(n, rng)
    values = f(x) / p(x)
    return float(values.mean()), float(values.std(ddof=1) / math.sqrt(n))


def exponential_sampler(rate=1.0):
    """Inverse-transform sampling of rate*exp(-rate x) on [0, infinity)."""
    def sample(n, rng):
        return -np.log(1.0 - rng.random(n)) / rate
    return sample


def box_muller(n, rng):
    """The normal distribution from two uniforms, by the polar transformation.

    x = sqrt(-2 ln u1) cos(2 pi u2) is the classic pair of transformations
    that turns the uniform generator into a Gaussian one -- the
    two-dimensional version of the inverse-CDF trick, used because the
    Gaussian CDF cannot be inverted in closed form.
    """
    u1, u2 = rng.random(n), rng.random(n)
    return np.sqrt(-2.0 * np.log(1.0 - u1)) * np.cos(2.0 * math.pi * u2)


def acceptance_rejection(f, n, rng, a=0.0, b=3.0, ceiling=None):
    """Von Neumann's method: throw darts at the box and count the hits."""
    if ceiling is None:
        ceiling = float(np.max(f(np.linspace(a, b, 1000))))
    x = a + (b - a) * rng.random(n)
    y = ceiling * rng.random(n)
    hits = np.count_nonzero(y < f(x))
    area = (b - a) * ceiling
    fraction = hits / n
    error = area * math.sqrt(fraction * (1.0 - fraction) / n)
    return area * fraction, error


def six_dimensional_brute_force(n, rng, cutoff=5.0):
    """int exp(-x^2-y^2)(x-y)^2 d^3x d^3y by uniform sampling in a box."""
    points = cutoff * (2.0 * rng.random((n, 6)) - 1.0)
    x, y = points[:, :3], points[:, 3:]
    weight = np.exp(-np.sum(x * x, axis=1) - np.sum(y * y, axis=1))
    values = weight * np.sum((x - y) ** 2, axis=1)
    volume = (2.0 * cutoff) ** 6
    return (float(volume * values.mean()),
            float(volume * values.std(ddof=1) / math.sqrt(n)))


def six_dimensional_importance(n, rng):
    """The same integral, sampling the Gaussian and averaging (x-y)^2.

    Writing exp(-x^2) = sqrt(pi) times the normal density of variance 1/2,
    the integral becomes pi^3 <(x-y)^2> over that density, and the integrand
    that is left is a low-order polynomial.
    """
    scale = 1.0 / math.sqrt(2.0)
    x = rng.normal(0.0, scale, (n, 3))
    y = rng.normal(0.0, scale, (n, 3))
    values = math.pi ** 3 * np.sum((x - y) ** 2, axis=1)
    return float(values.mean()), float(values.std(ddof=1) / math.sqrt(n))


def six_dimensional_exact():
    """3 pi^3, from <x^2> = <y^2> = 3/2 and <x.y> = 0."""
    return 3.0 * math.pi ** 3


# ---------------------------------------------------------------------------
#  Errors: the central limit theorem, correlation time, blocking
# ---------------------------------------------------------------------------
def autocorrelation(x, max_lag=None):
    """kappa_d = f_d / sigma^2, the normalised autocorrelation function."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    if max_lag is None:
        max_lag = min(n // 4, 2000)
    deviation = x - x.mean()
    variance = float(np.dot(deviation, deviation) / n)
    return np.array([float(np.dot(deviation[:n - d], deviation[d:]) / n
                           / variance) for d in range(max_lag)])


def correlation_time(x, max_lag=None):
    """tau = 1 + 2 sum kappa_d, truncated at the first negative kappa.

    The sum over all lags is dominated by noise at large d, so it is standard
    to cut it where the autocorrelation first goes negative.  The effective
    number of independent samples is then n / tau.
    """
    kappa = autocorrelation(x, max_lag)
    cut = np.argmax(kappa[1:] < 0.0) + 1 if np.any(kappa[1:] < 0.0) else len(kappa)
    return 1.0 + 2.0 * float(np.sum(kappa[1:cut]))


def naive_error(x):
    """sigma / sqrt(n) -- correct only for independent samples."""
    x = np.asarray(x, dtype=float)
    return float(x.std(ddof=1) / math.sqrt(len(x)))


def corrected_error(x):
    """sqrt(tau) times the naive error."""
    return naive_error(x) * math.sqrt(correlation_time(x))


def blocking(x):
    """Flyvbjerg-Petersen blocking: the error as data are blocked together.

    Repeatedly average neighbouring pairs.  The naive error of the blocked
    data grows and then plateaus; the plateau is the true error, and no
    autocorrelation function has to be computed to find it.
    """
    x = np.asarray(x, dtype=float)
    out = []
    while len(x) >= 4:
        out.append((len(x), naive_error(x)))
        if len(x) % 2:
            x = x[:-1]
        x = 0.5 * (x[0::2] + x[1::2])
    return out


# ---------------------------------------------------------------------------
#  Random walks, the diffusion equation and Markov chains
# ---------------------------------------------------------------------------
def random_walk(n_walkers, n_steps, rng, step=1.0, limit=None):
    """Positions of n_walkers after each of n_steps unbiased jumps.

    With `limit` set the lattice is finite, from -limit to +limit, with
    periodic boundaries: a walker leaving one end reappears at the other.  On
    an infinite line the distribution spreads for ever; on a finite ring it
    reaches a steady state, which is what makes the entropy plateau.
    """
    moves = step * np.where(rng.random((n_steps, n_walkers)) < 0.5, -1.0, 1.0)
    positions = np.cumsum(moves, axis=0)
    if limit is not None:
        span = 2 * limit + 1
        positions = ((positions + limit) % span) - limit
    return positions


def walk_entropy(positions, n_bins=101, limit=50.5):
    """S = -sum w ln w for the walker distribution at each time step."""
    edges = np.linspace(-limit, limit, n_bins + 1)
    out = []
    for row in positions:
        counts, _ = np.histogram(row, bins=edges)
        w = counts / counts.sum()
        w = w[w > 0]
        out.append(float(-np.sum(w * np.log(w))))
    return np.array(out)


def transition_matrix(n_sites):
    """The one-dimensional walk as a Markov matrix with reflecting ends."""
    W = np.zeros((n_sites, n_sites))
    for j in range(n_sites):
        for i in (j - 1, j + 1):
            if 0 <= i < n_sites:
                W[i, j] += 0.5
        if j == 0:
            W[0, 0] += 0.5
        if j == n_sites - 1:
            W[-1, -1] += 0.5
    return W


def steady_state(W, tol=1e-12, max_iter=100000):
    """Iterate w <- W w until it stops changing; also the eigenvector route."""
    w = np.zeros(W.shape[0])
    w[0] = 1.0
    for iteration in range(1, max_iter + 1):
        new = W @ w
        if np.abs(new - w).max() < tol:
            break
        w = new
    values, vectors = np.linalg.eig(W)
    k = int(np.argmin(np.abs(values - 1.0)))
    eigenvector = np.real(vectors[:, k])
    eigenvector = eigenvector / eigenvector.sum()
    return new, eigenvector, iteration


EXAMPLE_W = np.array([[1/4, 1/9, 3/8, 1/3],
                      [2/4, 2/9, 0.0, 1/3],
                      [0.0, 1/9, 3/8, 0.0],
                      [1/4, 5/9, 2/8, 1/3]])


# ---------------------------------------------------------------------------
#  The Metropolis algorithm
# ---------------------------------------------------------------------------
def metropolis(log_target, x0, n_samples, rng, step=1.0, burn_in=1000):
    """Sample from a distribution known only up to normalisation.

    Symmetric proposals, so the acceptance ratio is just the ratio of the
    target densities and the normalisation never appears.  Returns the chain
    and the acceptance rate.
    """
    x = np.atleast_1d(np.asarray(x0, dtype=float)).copy()
    current = log_target(x)
    chain = np.empty((n_samples,) + x.shape)
    accepted = 0
    for k in range(n_samples + burn_in):
        trial = x + step * (2.0 * rng.random(x.shape) - 1.0)
        proposed = log_target(trial)
        if proposed - current > math.log(rng.random() + 1e-300):
            x, current = trial, proposed
            accepted += 1
        if k >= burn_in:
            chain[k - burn_in] = x
    return chain, accepted / (n_samples + burn_in)


def detailed_balance_residual(W, w):
    """max |W_ij w_j - W_ji w_i|, zero when detailed balance holds."""
    flux = W * w[None, :]
    return float(np.abs(flux - flux.T).max())


# ---------------------------------------------------------------------------
#  Variational Monte Carlo on the two-electron quantum dot of chapter 11
# ---------------------------------------------------------------------------
#  Trial function
#
#      Psi_T = exp(-alpha omega (r1^2 + r2^2) / 2)
#              * exp( a r12 / (1 + beta r12) )
#
#  The first factor is the lowest oscillator orbital for each electron, scaled
#  by a variational parameter; the second is a Jastrow factor with the cusp
#  parameter a = 1, the value the Coulomb cusp condition fixes for two
#  electrons of opposite spin in two dimensions.  With `cusp = 0` the Jastrow
#  factor is switched off entirely, and then alpha = 1 makes the trial state
#  exactly the minimal-basis Hartree-Fock determinant of chapter 11.  Note
#  that beta = 0 does *not* switch the Jastrow factor off: it leaves
#  exp(r12), which drives the electrons apart.
# ---------------------------------------------------------------------------
def log_psi(positions, alpha, beta, omega=1.0, cusp=0.0):
    """ln|Psi_T| for one configuration, positions of shape (2, 2)."""
    r2 = np.sum(positions * positions)
    value = -0.5 * alpha * omega * r2
    if cusp:
        r12 = math.dist(positions[0], positions[1])
        value += cusp * r12 / (1.0 + beta * r12)
    return value


def local_energy_simple(positions, alpha, omega=1.0):
    """E_L for the pure Gaussian trial function, analytically.

        E_L = 2 alpha omega + (1/2) omega^2 (1 - alpha^2) (r1^2 + r2^2)
              + 1 / r12

    At alpha = 1 the first two terms collapse to 2 omega and the local energy
    is 2 omega + 1/r12: the trial state is an eigenstate of the
    non-interacting problem and all the variance comes from the interaction.
    """
    r2 = float(np.sum(positions * positions))
    r12 = math.dist(positions[0], positions[1])
    return (2.0 * alpha * omega
            + 0.5 * omega ** 2 * (1.0 - alpha ** 2) * r2
            + 1.0 / r12)


def local_energy_numerical(positions, alpha, beta, omega=1.0, h=1e-4,
                           cusp=0.0):
    """E_L by finite differences, valid for any trial function.

        E_L = -(1/2) sum_i grad_i^2 Psi / Psi + V

    Computed from ln Psi so that the wave function itself never overflows;
    the identity is  grad^2 Psi / Psi = grad^2 ln Psi + (grad ln Psi)^2.
    """
    laplacian = 0.0
    gradient_squared = 0.0
    centre = log_psi(positions, alpha, beta, omega, cusp)
    for i in range(positions.shape[0]):
        for d in range(positions.shape[1]):
            shifted = positions.copy()
            shifted[i, d] += h
            plus = log_psi(shifted, alpha, beta, omega, cusp)
            shifted[i, d] -= 2.0 * h
            minus = log_psi(shifted, alpha, beta, omega, cusp)
            laplacian += (plus + minus - 2.0 * centre) / h ** 2
            gradient_squared += ((plus - minus) / (2.0 * h)) ** 2
    kinetic = -0.5 * (laplacian + gradient_squared)
    r2 = float(np.sum(positions * positions))
    r12 = math.dist(positions[0], positions[1])
    return kinetic + 0.5 * omega ** 2 * r2 + 1.0 / r12


def vmc(alpha, beta=0.0, jastrow=False, n_samples=200000, rng=None, omega=1.0,
        step=1.2, burn_in=5000):
    """Brute-force variational Monte Carlo for the two-electron dot.

    With ``jastrow=False`` the trial function is the pure Gaussian and the
    local energy is the analytic one; with ``jastrow=True`` the cusp factor is
    switched on and the local energy is evaluated by finite differences, which
    costs more but cannot be got algebraically wrong.

    Returns a dictionary with the energy, its naive and correlation-corrected
    errors, the variance of the local energy, the correlation time and the
    acceptance rate.  The variance is as informative as the energy: it
    vanishes when the trial function is exact.
    """
    if rng is None:
        rng = np.random.default_rng(2024)
    cusp = 1.0 if jastrow else 0.0

    def energy(p):
        return (local_energy_numerical(p, alpha, beta, omega, cusp=cusp)
                if jastrow else local_energy_simple(p, alpha, omega))

    positions = rng.normal(0.0, 1.0, (2, 2))
    current = 2.0 * log_psi(positions, alpha, beta, omega, cusp)
    samples = np.empty(n_samples)
    accepted = 0
    for k in range(n_samples + burn_in):
        trial = positions + step * (2.0 * rng.random((2, 2)) - 1.0)
        proposed = 2.0 * log_psi(trial, alpha, beta, omega, cusp)
        if proposed - current > math.log(rng.random() + 1e-300):
            positions, current = trial, proposed
            accepted += 1
        if k >= burn_in:
            samples[k - burn_in] = energy(positions)

    tau = correlation_time(samples, max_lag=500)
    return dict(energy=float(samples.mean()),
                naive_error=naive_error(samples),
                error=naive_error(samples) * math.sqrt(tau),
                variance=float(samples.var(ddof=1)),
                tau=tau,
                acceptance=accepted / (n_samples + burn_in),
                samples=samples)


def hartree_fock_minimal_basis(omega=1.0):
    """2 hbar omega + sqrt(pi/2) sqrt(hbar omega), the first row of table 11.1."""
    return 2.0 * omega + math.sqrt(math.pi / 2.0) * math.sqrt(omega)


def gaussian_energy_exact(alpha, omega=1.0):
    """<E> for the pure Gaussian trial function, in closed form.

    Under |Psi_T|^2 = exp(-alpha omega (r1^2+r2^2)) the two averages needed
    are <r1^2+r2^2> = 2/(alpha omega) and, by scaling lengths,
    <1/r12> = sqrt(pi/2) sqrt(alpha omega).  Hence

        E(alpha) = omega alpha + omega/alpha + sqrt(pi/2) sqrt(alpha omega),

    which reduces to the minimal-basis Hartree-Fock energy at alpha = 1 and
    has its minimum at alpha ~ 0.7628.  Every VMC number in the alpha scan can
    therefore be checked against an exact formula.
    """
    return (omega * alpha + omega / alpha
            + math.sqrt(math.pi / 2.0) * math.sqrt(alpha * omega))


def gaussian_optimum(omega=1.0):
    """The alpha minimising gaussian_energy_exact, by bisection on dE/dalpha."""
    def derivative(a):
        return (omega - omega / a ** 2
                + 0.5 * math.sqrt(math.pi / 2.0) * math.sqrt(omega / a))
    lo, hi = 0.3, 1.5
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if derivative(mid) > 0.0:
            hi = mid
        else:
            lo = mid
    alpha = 0.5 * (lo + hi)
    return alpha, gaussian_energy_exact(alpha, omega)


# ---------------------------------------------------------------------------
def _demo():
    rng = np.random.default_rng(2024)

    print("=" * 74)
    print("1. Random numbers")
    print("=" * 74)
    print("A generator can be perfectly uniform and still be useless.  RANDU")
    print("passes every one-dimensional test; its triples lie on planes.")
    print()
    randu = LinearCongruential(seed=1, a=65539, c=0, modulus=2**31)
    x_randu = randu.sample(20000)
    x_good = rng.random(20000)
    print(f"{'generator':>12s} {'chi2 (10 bins)':>16s} {'lag-1 corr':>13s} "
          f"{'on the planes':>15s}")
    for name, x, generator in (("RANDU", x_randu, LinearCongruential(1)),
                               ("PCG64", x_good, None)):
        planes = (planes_test(generator) if generator is not None
                  else float(np.mean(np.abs(
                      (x[2:] - 6*x[1:-1] + 9*x[:-2]) % 1.0) < 1e-12)))
        print(f"{name:>12s} {uniformity_test(x):16.2f} "
              f"{serial_correlation(x):13.5f} {planes:14.1%}")
    print()
    print("Both look uniform and uncorrelated at lag one.  Every RANDU triple")
    print("satisfies x_{n+2} = 6 x_{n+1} - 9 x_n exactly, so the points fill")
    print("only fifteen planes of the unit cube; a three-dimensional integral")
    print("sampled with it is systematically wrong.  We use numpy's PCG64")
    print("throughout, seeded explicitly so that every run is reproducible.")

    print()
    print("=" * 74)
    print("2. Monte Carlo integration and the square-root law")
    print("=" * 74)
    print("The error of a Monte Carlo estimate falls as 1/sqrt(N) regardless")
    print("of the number of dimensions.  That single fact is why the method")
    print("wins in high dimensions: a grid needs N^d points for the same")
    print("accuracy in d dimensions, Monte Carlo needs N.")
    print()
    exact = 4.0 * math.atan(1.0)

    def f(x):
        return 4.0 / (1.0 + x * x)

    print(f"{'N':>9s} {'estimate':>13s} {'error':>11s} {'|est - pi|':>12s} "
          f"{'error*sqrt(N)':>15s}")
    for n in (10**3, 10**4, 10**5, 10**6):
        value, error = brute_force(f, n, rng)
        print(f"{n:9d} {value:13.6f} {error:11.6f} {abs(value-exact):12.6f} "
              f"{error*math.sqrt(n):15.4f}")
    print()
    print("The last column is constant: the error really is sigma/sqrt(N),")
    print(f"with sigma the standard deviation of the integrand ({0.6:.1f} here).")

    print()
    print("=" * 74)
    print("3. Changing variables: the inverse cumulative distribution")
    print("=" * 74)
    print("A generator gives uniform numbers on [0,1); the inverse of the")
    print("cumulative distribution turns them into anything integrable.")
    print()
    sample = exponential_sampler(1.0)(200000, rng)
    print(f"   exponential, rate 1: mean {sample.mean():.4f} (exact 1), "
          f"variance {sample.var():.4f} (exact 1)")
    gauss = box_muller(200000, rng)
    print(f"   Box-Muller normal:   mean {gauss.mean():+.4f} (exact 0), "
          f"variance {gauss.var():.4f} (exact 1)")
    value, error = acceptance_rejection(np.exp, 200000, rng, 0.0, 3.0)
    print(f"   acceptance-rejection for int_0^3 exp(x) dx: "
          f"{value:.4f} +/- {error:.4f}   (exact {math.exp(3)-1:.4f})")

    print()
    print("=" * 74)
    print("4. Importance sampling")
    print("=" * 74)
    print("If the integrand is peaked, uniform sampling wastes most of its")
    print("effort.  Sampling from a density that resembles the integrand and")
    print("averaging f/p instead flattens what is actually averaged, and the")
    print("variance drops.  The integral is")
    print("   int exp(-x^2-y^2) (x-y)^2 d^3x d^3y = 3 pi^3 = "
          f"{six_dimensional_exact():.6f}")
    print()
    print(f"{'N':>9s} {'brute force':>26s} {'importance sampling':>26s}")
    for n in (10**4, 10**5, 10**6):
        b_value, b_error = six_dimensional_brute_force(n, rng)
        i_value, i_error = six_dimensional_importance(n, rng)
        print(f"{n:9d} {b_value:14.4f} +/- {b_error:7.4f} "
              f"{i_value:14.4f} +/- {i_error:7.4f}")
    print()
    b_value, b_error = six_dimensional_brute_force(10**6, rng)
    i_value, i_error = six_dimensional_importance(10**6, rng)
    print(f"At N = 10^6 importance sampling is {b_error/i_error:.0f} times")
    print("more accurate, which is a factor of about "
          f"{(b_error/i_error)**2:.0f} in computer time.")

    print()
    print("=" * 74)
    print("5. The central limit theorem")
    print("=" * 74)
    print("Averages of m samples from any distribution with a finite variance")
    print("are normally distributed with variance sigma^2/m.  Here the parent")
    print("distribution is uniform on [0,1), which is not remotely Gaussian.")
    print()
    print(f"{'m':>6s} {'mean of averages':>18s} {'variance':>12s} "
          f"{'sigma^2/m':>12s} {'ratio':>8s}")
    parent_variance = 1.0 / 12.0
    for m in (1, 2, 5, 10, 50):
        averages = rng.random((20000, m)).mean(axis=1)
        predicted = parent_variance / m
        print(f"{m:6d} {averages.mean():18.6f} {averages.var():12.6f} "
              f"{predicted:12.6f} {averages.var()/predicted:8.4f}")

    print()
    print("=" * 74)
    print("6. Correlated data and the true error")
    print("=" * 74)
    print("The 1/sqrt(N) law assumes independent samples.  A Markov chain")
    print("does not produce them, and using sigma/sqrt(N) then understates the")
    print("error -- the one mistake that makes a Monte Carlo result look")
    print("better than it is.")
    print()
    chain, acceptance = metropolis(
        lambda x: -0.5 * float(np.dot(x, x)), np.zeros(1), 100000, rng,
        step=0.5)
    values = chain[:, 0]
    tau = correlation_time(values)
    print(f"   Metropolis sampling of a unit Gaussian, step 0.5")
    print(f"   acceptance rate     : {acceptance:.3f}")
    print(f"   mean                : {values.mean():+.5f}  (exact 0)")
    print(f"   variance            : {values.var():.5f}   (exact 1)")
    print(f"   correlation time    : {tau:.2f}")
    print(f"   naive error         : {naive_error(values):.6f}")
    print(f"   corrected error     : {corrected_error(values):.6f}")
    print(f"   effective samples   : {len(values)/tau:.0f} of {len(values)}")
    print()
    print("   blocking: the error grows and then plateaus at the true value")
    print(f"{'block size':>12s} {'samples':>10s} {'error':>12s}")
    for size, (n, error) in enumerate(blocking(values)):
        if size % 3 == 0 and n >= 16:
            print(f"{len(values)//n:12d} {n:10d} {error:12.6f}")

    print()
    print("=" * 74)
    print("7. Random walks and the diffusion equation")
    print("=" * 74)
    print("An unbiased walk has <x> = 0 and <x^2> = n l^2, so the walkers")
    print("spread as sqrt(t) and the distribution tends to a Gaussian.  This")
    print("is the diffusion equation with D = l^2 / 2 epsilon.")
    print()
    positions = random_walk(4000, 400, rng)
    print(f"{'steps':>7s} {'<x>':>10s} {'<x^2>':>10s} {'n l^2':>10s} "
          f"{'ratio':>8s}")
    for n in (25, 50, 100, 200, 400):
        row = positions[n - 1]
        print(f"{n:7d} {row.mean():+10.4f} {(row**2).mean():10.3f} "
              f"{float(n):10.3f} {(row**2).mean()/n:8.4f}")
    print()
    print("   On a finite ring the spreading has to stop.  With 4000 walkers")
    print("   on a lattice from -50 to +50 and periodic boundaries, the")
    print("   entropy S = -sum w ln w rises and then flattens:")
    print()
    entropy = walk_entropy(random_walk(4000, 3000, rng, limit=50))
    for n in (1, 10, 100, 500, 1000, 2000, 3000):
        print(f"      after {n:4d} steps: S = {entropy[n-1]:.4f}")
    print(f"      uniform distribution: S = {math.log(101):.4f}")
    print()
    print("   The walkers have spread over all available sites and the")
    print("   system has reached a steady state.  Equilibrium is not imposed")
    print("   anywhere -- it happens because the walk explores every state,")
    print("   which is the ergodic hypothesis, and it is the reason a Markov")
    print("   chain can be used to sample a distribution at all.")

    print()
    print("=" * 74)
    print("8. Markov chains and the steady state")
    print("=" * 74)
    print("A Markov chain is a transition matrix acting on a probability")
    print("vector.  Its steady state is the eigenvector with eigenvalue one,")
    print("and iterating from any starting vector converges to it.")
    print()
    final, eigenvector, iterations = steady_state(EXAMPLE_W)
    print("   the 4x4 example matrix, starting from w = (1,0,0,0):")
    print(f"      iterated to convergence in {iterations} steps:")
    print("        " + "  ".join(f"{v:.6f}" for v in final))
    print("      eigenvector of W with eigenvalue 1:")
    print("        " + "  ".join(f"{v:.6f}" for v in eigenvector))
    print(f"      agree: {np.allclose(final, eigenvector, atol=1e-8)}")
    print()
    W = transition_matrix(11)
    final, eigenvector, iterations = steady_state(W)
    print("   the eleven-site random walk with reflecting ends:")
    print(f"      steady state is uniform: "
          f"{np.allclose(final, np.ones(11)/11, atol=1e-8)}")
    print(f"      detailed balance residual: "
          f"{detailed_balance_residual(W, final):.1e}")

    print()
    print("=" * 74)
    print("9. The Metropolis algorithm")
    print("=" * 74)
    print("We can rarely sample a distribution directly, but we can almost")
    print("always evaluate it up to a constant.  Metropolis needs only the")
    print("ratio w_j/w_i, so the normalisation never has to be computed --")
    print("which for a many-body system is the whole difficulty.")
    print()
    print("   sampling exp(-x^2/2 - x^4/10), normalisation unknown:")

    def log_target(x):
        value = float(x[0])
        return -0.5 * value ** 2 - 0.1 * value ** 4

    for step in (0.1, 0.5, 2.0, 8.0, 32.0):
        chain, acceptance = metropolis(log_target, np.zeros(1), 200000, rng,
                                       step=step)
        values = chain[:, 0]
        tau = correlation_time(values)
        print(f"      step {step:5.1f}: acceptance {acceptance:5.1%}, "
              f"tau {tau:7.2f}, <x^2> = {(values**2).mean():.4f} "
              f"+/- {corrected_error(values**2):.4f}")
    print()
    print("   A small step is accepted almost always but moves nowhere; a")
    print("   large one is nearly always rejected and the walker stays put.")
    print("   Both give long correlation times, and the best step is in")
    print("   between.  The usual rule of thumb is to tune the step for an")
    print("   acceptance rate somewhere near one half.")

    print()
    print("=" * 74)
    print("10. Variational Monte Carlo on the quantum dot of chapter 11")
    print("=" * 74)
    print("Everything above now pays off.  The trial function is")
    print("   Psi_T = exp(-alpha omega (r1^2+r2^2)/2) "
          "exp(r12 / (1 + beta r12)),")
    print("the energy is the average of the local energy E_L = H Psi_T / Psi_T")
    print("over |Psi_T|^2, and the average is taken by Metropolis.")
    print()
    print("First a check.  At alpha = 1 and beta = 0 the trial function IS")
    print("the minimal-basis Hartree-Fock determinant of chapter 11, so the")
    print("energy must be 2 hbar omega + sqrt(pi/2) sqrt(hbar omega).")
    print()
    reference = hartree_fock_minimal_basis()
    result = vmc(alpha=1.0, jastrow=False, n_samples=200000, rng=rng,
                 step=1.0)
    print(f"   VMC      : {result['energy']:.6f} +/- {result['error']:.6f}")
    print(f"   exact    : {reference:.6f}   (table 11.1, first row)")
    print(f"   deviation: {abs(result['energy']-reference)/result['error']:.2f}"
          " standard errors")
    print(f"   variance of E_L {result['variance']:.4f}, "
          f"acceptance {result['acceptance']:.2f}, tau {result['tau']:.1f}")
    print()
    print("Now vary alpha, still without the Jastrow factor.  Here the average")
    print("can be done analytically as well,")
    print("   E(alpha) = omega alpha + omega/alpha "
          "+ sqrt(pi/2) sqrt(alpha omega),")
    print("so every Monte Carlo number can be checked against a formula.")
    print()
    optimum, optimum_energy = gaussian_optimum()
    print(f"{'alpha':>8s} {'VMC':>12s} {'error':>10s} {'exact':>12s} "
          f"{'deviation':>10s} {'variance':>10s}")
    for alpha in (0.70, optimum, 0.85, 1.00):
        r = vmc(alpha=alpha, jastrow=False, n_samples=200000, rng=rng,
                step=1.0)
        exact_value = gaussian_energy_exact(alpha)
        print(f"{alpha:8.4f} {r['energy']:12.6f} {r['error']:10.6f} "
              f"{exact_value:12.6f} "
              f"{abs(r['energy']-exact_value)/r['error']:9.2f}s "
              f"{r['variance']:10.4f}")
    print()
    print(f"The best pure Gaussian has alpha = {optimum:.4f} and "
          f"E = {optimum_energy:.6f},")
    print("below the Hartree-Fock value because a scaled orbital beats the")
    print("bare oscillator one.  The variance stays around 3, and it is")
    print("dominated by the rare configurations in which the two electrons")
    print("come close and 1/r12 blows up.  That divergence is what the")
    print("Jastrow factor is for.")
    print()
    print("Switching it on, with the cusp parameter fixed at a = 1 by the")
    print("Coulomb cusp condition and beta variational:")
    print()
    print(f"{'alpha':>7s} {'beta':>7s} {'energy':>12s} {'error':>10s} "
          f"{'variance':>11s}")
    for alpha, beta in ((0.98, 0.30), (0.98, 0.40), (0.98, 0.48),
                        (0.96, 0.44), (1.00, 0.44)):
        r = vmc(alpha=alpha, beta=beta, jastrow=True, n_samples=60000,
                rng=np.random.default_rng(2024), step=1.0)
        print(f"{alpha:7.2f} {beta:7.2f} {r['energy']:12.6f} "
              f"{r['error']:10.6f} {r['variance']:11.5f}")
    print()
    print("   exact (Taut)          : 3.000000")
    print("   CCSD, 42 orbitals     : 3.013626   (table 11.4)")
    print("   Hartree-Fock, 42 orb. : 3.161921   (table 11.3)")
    print()
    print("Two things to take away.  The variance has fallen by three orders")
    print("of magnitude, from about 3 to about 0.002: it would vanish for the")
    print("exact wave function, so it measures the quality of the trial")
    print("function directly, and it is a more sensitive diagnostic than the")
    print("energy.  And the energy itself is now below the CCSD result of")
    print("chapter 11 -- not because the many-body treatment is better, but")
    print("because there is no basis set here at all.  What we have instead")
    print("is a statistical error, and sections 5, 6 and 9 are what keep it")
    print("honest.")


if __name__ == "__main__":
    _demo()
