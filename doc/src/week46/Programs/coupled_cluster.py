from functools import partial

import jax
import jax.numpy as np
import numpy as onp
from jax import Array, jit

jax.config.update("jax_enable_x64", True)


def init_pairing_v(g: float, pnum: int, hnum: int) -> tuple[Array, Array, Array]:
    """Returns potential matrices of the pairing model in three relevant channels.

    Args:
        g (float): Strength of the pairing interaction, as in Eq. (8.42).
        pnum (int): Number of particle states.
        hnum (int): Number of hole states.

    Returns:
        tuple: A tuple containing three numpy arrays:
        The interaction as a 4-indexed tensor in the:
        - v_pppp (np.array): particle-particle channel (pnum, pnum, pnum, pnum).
        - v_pphh (np.array): particle-hole channel (pnum, pnum, hnum, hnum).
        - v_hhhh (np.array): hole-hole channel (hnum, hnum, hnum, hnum).
    """

    v_pppp = onp.zeros((pnum, pnum, pnum, pnum))
    v_pphh = onp.zeros((pnum, pnum, hnum, hnum))
    v_hhhh = onp.zeros((hnum, hnum, hnum, hnum))

    gval = -0.5 * g

    # Even indices for particles and holes
    p_even = onp.arange(0, pnum, 2)
    h_even = onp.arange(0, hnum, 2)

    # For v_pppp
    for a in p_even:
        for b in p_even:
            v_pppp[a, a + 1, b, b + 1] = gval
            v_pppp[a + 1, a, b, b + 1] = -gval
            v_pppp[a, a + 1, b + 1, b] = -gval
            v_pppp[a + 1, a, b + 1, b] = gval

    # For v_pphh
    for a in p_even:
        for i in h_even:
            v_pphh[a, a + 1, i, i + 1] = gval
            v_pphh[a + 1, a, i, i + 1] = -gval
            v_pphh[a, a + 1, i + 1, i] = -gval
            v_pphh[a + 1, a, i + 1, i] = gval

    # For v_hhhh
    for j in h_even:
        for i in h_even:
            v_hhhh[j, j + 1, i, i + 1] = gval
            v_hhhh[j + 1, j, i, i + 1] = -gval
            v_hhhh[j, j + 1, i + 1, i] = -gval
            v_hhhh[j + 1, j, i + 1, i] = gval
    v_pppp, v_pphh, v_hhhh = map(np.array, (v_pppp, v_pphh, v_hhhh))

    return v_pppp, v_pphh, v_hhhh


def init_pairing_fock(
    delta: float, g: float, pnum: int, hnum: int
) -> tuple[Array, Array]:
    """Initialize the Fock matrix of the pairing model.

    Args:
        delta (float): Single-particle spacing, as in Eq. (8.41).
        g (float): Pairing strength, as in Eq. (8.42).
        pnum (int): Number of particle states.
        hnum (int): Number of hole states.

    Returns:
        tuple: A tuple containing two numpy arrays:
        The Fock matrix in two channels as numpy arrays
        np.array(pnum, pnum), np.array(hnum, hnum).
    """
    deltaval = 0.5 * delta
    gval = -0.5 * g
    f_pp = onp.zeros((pnum, pnum))
    f_hh = onp.zeros((hnum, hnum))

    # Set f_hh diagonals
    i_even = onp.arange(0, hnum, 2)
    f_hh[i_even, i_even] = deltaval * i_even + gval
    f_hh[i_even + 1, i_even + 1] = deltaval * i_even + gval

    # Set f_pp diagonals
    a_even = onp.arange(0, pnum, 2)
    f_pp[a_even, a_even] = deltaval * (hnum + a_even)
    f_pp[a_even + 1, a_even + 1] = deltaval * (hnum + a_even)

    f_pp, f_hh = np.array(f_pp), np.array(f_hh)
    return f_pp, f_hh


def init_t2(v_pphh: Array, f_pp: Array, f_hh: Array) -> Array:
    """Initialize T2 amplitudes as in MBPT2.

    Args:
        v_pphh (Array): Particle-hole interaction tensor.
        f_pp (Array): Fock matrix for particle states.
        f_hh (Array): Fock matrix for hole states.

    Returns:
        Array: Initial T2 amplitudes.
    """
    Fhh_diag = np.diag(f_hh)
    Fpp_diag = np.diag(f_pp)

    # denominator shape: (pnum, pnum, hnum, hnum)
    denom = (
        Fhh_diag[None, None, :, None]
        + Fhh_diag[None, None, None, :]
        - Fpp_diag[:, None, None, None]
        - Fpp_diag[None, :, None, None]
    )

    t2_new = v_pphh / denom
    return t2_new


@jit
def ccd_iter(
    v_pppp: Array, v_pphh: Array, v_hhhh: Array, f_pp: Array, f_hh: Array, t2: Array
) -> Array:
    """Perform one iteration of the CCD equations (8.34).

    Args:
        v_pppp (Array): Pairing tensor in pppp channel.
        v_pphh (Array): Pairing tensor in pphh channel.
        v_hhhh (Array): Pairing tensor in hhhh channel.
        f_pp (Array): Fock matrix in pp channel.
        f_hh (Array): Fock matrix in hh channel.
        t2 (Array): Initial T2 amplitude.

    Returns:
        Array: New T2 amplitude.
    """
    Hbar_pphh = (
        v_pphh
        + np.einsum("bc,acij->abij", f_pp, t2)
        - np.einsum("ac,bcij->abij", f_pp, t2)
        - np.einsum("abik,kj->abij", t2, f_hh)
        + np.einsum("abjk,ki->abij", t2, f_hh)
        + 0.5 * np.einsum("abcd,cdij->abij", v_pppp, t2)
        + 0.5 * np.einsum("abkl,klij->abij", t2, v_hhhh)
    )

    # hh intermediate, see (8.47)
    chi_hh = 0.5 * np.einsum("cdkl,cdjl->kj", v_pphh, t2)

    Hbar_pphh = Hbar_pphh - (
        np.einsum("abik,kj->abij", t2, chi_hh) - np.einsum("abik,kj->abji", t2, chi_hh)
    )

    # pp intermediate, see (8.46)
    chi_pp = -0.5 * np.einsum("cdkl,bdkl->cb", v_pphh, t2)

    Hbar_pphh = Hbar_pphh + (
        np.einsum("acij,cb->abij", t2, chi_pp) - np.einsum("acij,cb->baij", t2, chi_pp)
    )

    # hhhh intermediate, see (8.48)
    chi_hhhh = 0.5 * np.einsum("cdkl,cdij->klij", v_pphh, t2)

    Hbar_pphh = Hbar_pphh + 0.5 * np.einsum("abkl,klij->abij", t2, chi_hhhh)

    # phph intermediate, see (8.49)
    chi_phph = +0.5 * np.einsum("cdkl,dblj->bkcj", v_pphh, t2)

    Hbar_pphh = Hbar_pphh + (
        np.einsum("bkcj,acik->abij", chi_phph, t2)
        - np.einsum("bkcj,acik->baij", chi_phph, t2)
        - np.einsum("bkcj,acik->abji", chi_phph, t2)
        + np.einsum("bkcj,acik->baji", chi_phph, t2)
    )

    Fhh_diag = np.diag(f_hh)
    Fpp_diag = np.diag(f_pp)
    denom = (
        Fhh_diag[None, None, :, None]
        + Fhh_diag[None, None, None, :]
        - Fpp_diag[:, None, None, None]
        - Fpp_diag[None, :, None, None]
    )

    t2_new = t2 + Hbar_pphh / denom
    return t2_new


@jit
def ccd_energy(v_pphh: Array, t2: Array) -> Array:
    """Compute CCD energy.

    Args:
        v_pphh (Array): Pairing tensor in pphh channel.
        t2 (Array): T2 amplitude.

    Returns:
        Array: CCD energy.
    """
    return 0.25 * np.sum(v_pphh * t2)


def calculate_ccd_energy(
    pnum: int, hnum: int, delta: float, g: float, verbose: bool = False
) -> float:
    """Calculate the CCD energy for the pairing model.

    Args:
        pnum (int): Number of particle states.
        hnum (int): Number of hole states.
        delta (float): Single-particle spacing.
        g (float): Pairing strength.
        verbose (bool): Whether to print the results.

    Returns:
        float: The CCD energy.
    """
    if verbose:
        print("Parameters:")
        print(f"{delta = }, {g = }, {pnum = }, {hnum = }")

    # Initialize pairing matrix elements and Fock matrix
    v_pppp, v_pphh, v_hhhh = init_pairing_v(g, pnum, hnum)
    f_pp, f_hh = init_pairing_fock(delta, g, pnum, hnum)

    # Initialize T2 amplitudes from MBPT2
    t2 = init_t2(v_pphh, f_pp, f_hh)
    mpbt_erg = ccd_energy(v_pphh, t2)

    # Exact MBPT2 for comparison
    exact_mbpt2 = -0.25 * g**2 * (1.0 / (2.0 + g) + 2.0 / (4.0 + g) + 1.0 / (6.0 + g))

    # iterate CCD equations niter times
    niter = 200
    mix = 0.3
    erg_old = 0.0
    eps = np.float64(1.0e-10)
    for iter_count in range(niter):
        t2_new = ccd_iter(v_pppp, v_pphh, v_hhhh, f_pp, f_hh, t2)
        erg = ccd_energy(v_pphh, t2_new)
        myeps = np.abs(erg - erg_old) / np.abs(erg)
        if myeps < eps:
            t2 = t2_new
            num_iter = iter_count
            break
        erg_old = erg
        t2 = mix * t2_new + (1.0 - mix) * t2

    if verbose:
        print(f"MPBT2 energy = {mpbt_erg}, compared to exact: {exact_mbpt2}")
        print(f"Converged in {str(num_iter):3s} iterations, Energy = {erg}")
    return erg


def CCD_energies(g_values: float) -> Array:
    pnum = 4  # number of particle states
    hnum = 4  # number of hole states
    delta = 1.0

    par_main = partial(calculate_ccd_energy, pnum, hnum, delta)

    energies = onp.zeros_like(g_values)
    for i, gval in enumerate(g_values):
        delta_erg = par_main(g=gval)
        energies[i] = delta_erg

    return 2 - g_values + energies


if __name__ == "__main__":
    # set parameters as for model
    pnum = 4  # number of particle states
    hnum = 4  # number of hole states
    delta = 1.0

    par_main = partial(calculate_ccd_energy, pnum, hnum, delta, verbose=True)

    g = np.linspace(-1, 1, 100)
    energies = onp.zeros_like(g)
    for i, gval in enumerate(g):
        delta_erg = par_main(g=gval)
        energies[i] = delta_erg

    import matplotlib.pyplot as plt

    plt.plot(g, 2 - g + energies)
    plt.xlabel("g")
    plt.ylabel("Energy")
    plt.show()
