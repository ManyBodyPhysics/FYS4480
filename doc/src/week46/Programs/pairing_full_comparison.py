#!/usr/bin/env python3
"""
Coupled Cluster Doubles (CCD) for the Pairing Model

Comparison of:
1. Exact Diagonalization
2. MBPT2 
3. CCD

CORRECTED VERSION with proper single-particle energies from Fock matrix
"""

from functools import partial
import jax
import jax.numpy as jnp
import numpy as np
from jax import Array, jit
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.linalg import eigh

jax.config.update("jax_enable_x64", True)

# ==============================================================================
# 1. EXACT DIAGONALIZATION (CORRECTED)
# ==============================================================================

def exact_diagonalization_pairing(pnum, hnum, delta, g):
    """
    Exact diagonalization with CORRECTED single-particle energies.
    
    Key fix: Single-particle energies must match Fock matrix (includes mean-field pairing)
    """
    n_levels = 4  # 4 doubly degenerate levels
    n_pairs = 2   # 4 fermions = 2 pairs
    
    # Basis: all ways to place 2 pairs in 4 levels
    basis_states = list(combinations(range(n_levels), n_pairs))
    n_basis = len(basis_states)  # C(4,2) = 6
    
    # CORRECTED: Single-particle energies from Fock matrix
    deltaval = 0.5 * delta
    gval = -0.5 * g  # Mean-field pairing contribution
    
    epsilon = np.zeros(n_levels)
    # Hole levels
    epsilon[0] = deltaval * 0 + gval
    epsilon[1] = deltaval * 2 + gval
    # Particle levels
    epsilon[2] = deltaval * (hnum + 0)
    epsilon[3] = deltaval * (hnum + 2)
    
    # Construct Hamiltonian
    H = np.zeros((n_basis, n_basis))
    
    for i, state_i in enumerate(basis_states):
        for j, state_j in enumerate(basis_states):
            # Diagonal: sum of pair energies
            if i == j:
                H[i, j] = sum(2 * epsilon[k] for k in state_i)
            
            # Off-diagonal: pairing interaction
            for k in range(n_levels):
                for l in range(n_levels):
                    if l in state_j and k not in state_j:
                        state_temp = list(state_j)
                        state_temp.remove(l)
                        state_temp.append(k)
                        state_temp = tuple(sorted(state_temp))
                        if state_temp == state_i:
                            H[i, j] -= g / 2.0
    
    eigenvalues, eigenvectors = eigh(H)
    return eigenvalues, eigenvectors, basis_states, H


# ==============================================================================
# 2. CCD IMPLEMENTATION (Original Code)
# ==============================================================================

def init_pairing_v(g, pnum, hnum):
    """Initialize pairing interaction tensors."""
    v_pppp = np.zeros((pnum, pnum, pnum, pnum))
    v_pphh = np.zeros((pnum, pnum, hnum, hnum))
    v_hhhh = np.zeros((hnum, hnum, hnum, hnum))
    gval = -0.5 * g
    p_even = np.arange(0, pnum, 2)
    h_even = np.arange(0, hnum, 2)
    
    for a in p_even:
        for b in p_even:
            v_pppp[a, a + 1, b, b + 1] = gval
            v_pppp[a + 1, a, b, b + 1] = -gval
            v_pppp[a, a + 1, b + 1, b] = -gval
            v_pppp[a + 1, a, b + 1, b] = gval
    
    for a in p_even:
        for i in h_even:
            v_pphh[a, a + 1, i, i + 1] = gval
            v_pphh[a + 1, a, i, i + 1] = -gval
            v_pphh[a, a + 1, i + 1, i] = -gval
            v_pphh[a + 1, a, i + 1, i] = gval
    
    for j in h_even:
        for i in h_even:
            v_hhhh[j, j + 1, i, i + 1] = gval
            v_hhhh[j + 1, j, i, i + 1] = -gval
            v_hhhh[j, j + 1, i + 1, i] = -gval
            v_hhhh[j + 1, j, i + 1, i] = gval
    
    v_pppp, v_pphh, v_hhhh = map(jnp.array, (v_pppp, v_pphh, v_hhhh))
    return v_pppp, v_pphh, v_hhhh


def init_pairing_fock(delta, g, pnum, hnum):
    """Initialize Fock matrix (includes mean-field pairing)."""
    deltaval = 0.5 * delta
    gval = -0.5 * g
    f_pp = np.zeros((pnum, pnum))
    f_hh = np.zeros((hnum, hnum))
    
    i_even = np.arange(0, hnum, 2)
    f_hh[i_even, i_even] = deltaval * i_even + gval
    f_hh[i_even + 1, i_even + 1] = deltaval * i_even + gval
    
    a_even = np.arange(0, pnum, 2)
    f_pp[a_even, a_even] = deltaval * (hnum + a_even)
    f_pp[a_even + 1, a_even + 1] = deltaval * (hnum + a_even)
    
    f_pp, f_hh = jnp.array(f_pp), jnp.array(f_hh)
    return f_pp, f_hh


def init_t2(v_pphh, f_pp, f_hh):
    """Initialize T2 amplitudes from MBPT2."""
    Fhh_diag = jnp.diag(f_hh)
    Fpp_diag = jnp.diag(f_pp)
    denom = (Fhh_diag[None, None, :, None] + Fhh_diag[None, None, None, :]
             - Fpp_diag[:, None, None, None] - Fpp_diag[None, :, None, None])
    return v_pphh / denom


@jit
def ccd_energy(v_pphh, t2):
    """Compute CCD correlation energy."""
    return 0.25 * jnp.sum(v_pphh * t2)


@jit
def ccd_iter(v_pppp, v_pphh, v_hhhh, f_pp, f_hh, t2):
    """One CCD iteration."""
    Hbar_pphh = (v_pphh + jnp.einsum("bc,acij->abij", f_pp, t2)
                 - jnp.einsum("ac,bcij->abij", f_pp, t2)
                 - jnp.einsum("abik,kj->abij", t2, f_hh)
                 + jnp.einsum("abjk,ki->abij", t2, f_hh)
                 + 0.5 * jnp.einsum("abcd,cdij->abij", v_pppp, t2)
                 + 0.5 * jnp.einsum("abkl,klij->abij", t2, v_hhhh))
    
    chi_hh = 0.5 * jnp.einsum("cdkl,cdjl->kj", v_pphh, t2)
    Hbar_pphh -= jnp.einsum("abik,kj->abij", t2, chi_hh) - jnp.einsum("abik,kj->abji", t2, chi_hh)
    
    chi_pp = -0.5 * jnp.einsum("cdkl,bdkl->cb", v_pphh, t2)
    Hbar_pphh += jnp.einsum("acij,cb->abij", t2, chi_pp) - jnp.einsum("acij,cb->baij", t2, chi_pp)
    
    chi_hhhh = 0.5 * jnp.einsum("cdkl,cdij->klij", v_pphh, t2)
    Hbar_pphh += 0.5 * jnp.einsum("abkl,klij->abij", t2, chi_hhhh)
    
    chi_phph = 0.5 * jnp.einsum("cdkl,dblj->bkcj", v_pphh, t2)
    Hbar_pphh += (jnp.einsum("bkcj,acik->abij", chi_phph, t2)
                  - jnp.einsum("bkcj,acik->baij", chi_phph, t2)
                  - jnp.einsum("bkcj,acik->abji", chi_phph, t2)
                  + jnp.einsum("bkcj,acik->baji", chi_phph, t2))
    
    Fhh_diag = jnp.diag(f_hh)
    Fpp_diag = jnp.diag(f_pp)
    denom = (Fhh_diag[None, None, :, None] + Fhh_diag[None, None, None, :]
             - Fpp_diag[:, None, None, None] - Fpp_diag[None, :, None, None])
    
    return t2 + Hbar_pphh / denom


def calculate_ccd_energy(pnum, hnum, delta, g, verbose=False):
    """Calculate CCD and MBPT2 energies."""
    v_pppp, v_pphh, v_hhhh = init_pairing_v(g, pnum, hnum)
    f_pp, f_hh = init_pairing_fock(delta, g, pnum, hnum)
    
    t2 = init_t2(v_pphh, f_pp, f_hh)
    mbpt2_erg = ccd_energy(v_pphh, t2)
    
    niter, mix, erg_old, eps = 200, 0.3, 0.0, jnp.float64(1.0e-10)
    
    for iter_count in range(niter):
        t2_new = ccd_iter(v_pppp, v_pphh, v_hhhh, f_pp, f_hh, t2)
        erg = ccd_energy(v_pphh, t2_new)
        if jnp.abs(erg - erg_old) / jnp.abs(erg) < eps:
            if verbose:
                print(f"CCD converged in {iter_count} iterations")
            break
        erg_old = erg
        t2 = mix * t2_new + (1.0 - mix) * t2
    
    if verbose:
        print(f"MBPT2 correlation: {mbpt2_erg:.10f}")
        print(f"CCD correlation:   {erg:.10f}")
    
    return float(erg), float(mbpt2_erg)


# ==============================================================================
# 3. UNITARY COUPLED CLUSTER DOUBLES (UCCD) WITH TROTTERIZATION
# ==============================================================================

def calculate_uccd_energy(pnum, hnum, delta, g, n_trotter=4, verbose=False):
    """
    Calculate UCCD energy with Trotterization.
    
    UCCD uses the unitary ansatz: |ψ⟩ = exp(T - T†)|ψ_ref⟩
    where T - T† is anti-Hermitian (ensuring unitarity).
    
    Trotterization: exp(T - T†) ≈ [exp((T - T†)/n)]^n
    
    Args:
        pnum: Number of particle states
        hnum: Number of hole states
        delta: Single-particle spacing
        g: Pairing strength
        n_trotter: Number of Trotter steps
        verbose: Print convergence info
    
    Returns:
        E_uccd_corr: UCCD correlation energy
    """
    from scipy.linalg import expm
    from scipy.optimize import minimize
    
    # Get exact Hamiltonian in full basis
    eigenvalues, eigenvectors, basis_states, H_full = exact_diagonalization_pairing(
        pnum, hnum, delta, g)
    
    n_basis = len(basis_states)
    
    # Reference state: pairs in levels (0,1)
    psi_ref = np.zeros(n_basis)
    psi_ref[0] = 1.0
    
    # Generate double excitation operators
    double_excitations = []
    
    for i in range(4):  # 4 levels
        for j in range(4):
            if i == j:
                continue
            
            # T_ij: move pair from level j to level i
            T_op = np.zeros((n_basis, n_basis), dtype=complex)
            
            for idx_bra, state_bra in enumerate(basis_states):
                for idx_ket, state_ket in enumerate(basis_states):
                    if j in state_ket and i not in state_ket:
                        state_new = list(state_ket)
                        state_new.remove(j)
                        state_new.append(i)
                        state_new = tuple(sorted(state_new))
                        
                        if state_new == state_bra:
                            T_op[idx_bra, idx_ket] = 1.0
            
            if np.any(np.abs(T_op) > 1e-10):
                double_excitations.append((i, j, T_op))
    
    n_amplitudes = len(double_excitations)
    
    # UCCD energy function
    def uccd_energy_trotter(amplitudes):
        # Anti-Hermitian operator: T - T†
        T_anti = np.zeros((n_basis, n_basis), dtype=complex)
        
        for idx, (i, j, T_op) in enumerate(double_excitations):
            T_anti += amplitudes[idx] * (T_op - T_op.conj().T)
        
        # Trotterization
        step_op = expm(T_anti / n_trotter)
        
        psi_uccd = psi_ref.copy()
        for _ in range(n_trotter):
            psi_uccd = step_op @ psi_uccd
        
        psi_uccd = psi_uccd / np.linalg.norm(psi_uccd)
        
        E = np.real(np.vdot(psi_uccd, H_full @ psi_uccd))
        return E
    
    # Optimize
    initial_amps = np.random.uniform(-0.1, 0.1, n_amplitudes)
    
    result = minimize(
        uccd_energy_trotter,
        initial_amps,
        method='L-BFGS-B',
        bounds=[(-5, 5)] * n_amplitudes,
        options={'maxiter': 500, 'disp': False}
    )
    
    E_uccd_total = result.fun
    E_ref = H_full[0, 0]
    E_uccd_corr = E_uccd_total - E_ref
    
    if verbose:
        print(f"UCCD (Trotter={n_trotter}) converged in {result.nit} iterations")
        print(f"UCCD correlation: {E_uccd_corr:.10f}")
    
    return E_uccd_corr


# ==============================================================================
# 4. MAIN EXECUTION
# ==============================================================================

if __name__ == "__main__":
    # Parameters
    pnum, hnum, delta = 4, 4, 1.0
    
    print("="*80)
    print("PAIRING MODEL: EXACT vs MBPT2 vs CCD vs UCCD")
    print("="*80)
    print(f"\nSystem: {hnum} fermions, 4 doubly degenerate levels")
    print(f"Parameters: pnum={pnum}, hnum={hnum}, delta={delta}")
    
    # Single point test at g=0.5
    g_test = 0.5
    
    print(f"\n" + "="*80)
    print(f"SINGLE POINT TEST: g = {g_test}")
    print("="*80)
    
    eigenvalues, eigenvectors, basis_states, H = exact_diagonalization_pairing(
        pnum, hnum, delta, g_test)
    
    print(f"\nBasis dimension: {len(basis_states)}")
    print(f"Basis states: {basis_states}")
    print(f"\nHamiltonian ({H.shape[0]}×{H.shape[1]}):")
    print(H)
    
    print(f"\nEigenvalues:")
    for i, E in enumerate(eigenvalues):
        print(f"  E_{i}: {E:+.10f}")
    
    E_exact = eigenvalues[0]
    E_ref = H[0, 0]
    
    # Calculate all methods
    E_ccd_corr, E_mbpt2_corr = calculate_ccd_energy(pnum, hnum, delta, g_test, verbose=True)
    
    # UCCD with different Trotter steps
    trotter_steps = [1, 2, 4, 8]
    uccd_results = {}
    
    print(f"\n{'='*80}")
    print("UCCD WITH DIFFERENT TROTTER STEPS")
    print("="*80)
    
    for n_trotter in trotter_steps:
        E_uccd_corr = calculate_uccd_energy(pnum, hnum, delta, g_test, 
                                           n_trotter=n_trotter, verbose=True)
        uccd_results[n_trotter] = E_uccd_corr
    
    print(f"\n{'='*80}")
    print("COMPARISON")
    print("="*80)
    print(f"\n{'Method':20s} | {'Total E':>14s} | {'Corr E':>14s} | {'Error':>12s}")
    print("-"*75)
    print(f"{'Exact':20s} | {E_exact:+14.10f} | {E_exact-E_ref:+14.10f} | {'--':>12s}")
    print(f"{'MBPT2':20s} | {E_ref+E_mbpt2_corr:+14.10f} | {E_mbpt2_corr:+14.10f} | {abs(E_ref+E_mbpt2_corr-E_exact):12.2e}")
    print(f"{'CCD':20s} | {E_ref+E_ccd_corr:+14.10f} | {E_ccd_corr:+14.10f} | {abs(E_ref+E_ccd_corr-E_exact):12.2e}")
    
    for n_trotter in trotter_steps:
        E_uccd = E_ref + uccd_results[n_trotter]
        label = f"UCCD (Trotter={n_trotter})"
        print(f"{label:20s} | {E_uccd:+14.10f} | {uccd_results[n_trotter]:+14.10f} | {abs(E_uccd-E_exact):12.2e}")
    
    # Energy scan
    print(f"\n{'='*80}")
    print("ENERGY SCAN")
    print("="*80)
    
    g_values = np.linspace(-1, 1, 50)
    energies_exact = np.zeros_like(g_values)
    energies_ccd = np.zeros_like(g_values)
    energies_mbpt2 = np.zeros_like(g_values)
    energies_uccd = {}
    for nt in [1, 2, 4]:
        energies_uccd[nt] = np.zeros_like(g_values)
    
    print(f"\nComputing energies for {len(g_values)} values of g...")
    for i, g in enumerate(g_values):
        evals, _, _, H_mat = exact_diagonalization_pairing(pnum, hnum, delta, g)
        energies_exact[i] = evals[0]
        E_ccd_c, E_mbpt2_c = calculate_ccd_energy(pnum, hnum, delta, g)
        E_ref = H_mat[0, 0]
        energies_ccd[i] = E_ref + E_ccd_c
        energies_mbpt2[i] = E_ref + E_mbpt2_c
        
        # UCCD with different Trotter steps
        for nt in [1, 2, 4]:
            E_uccd_c = calculate_uccd_energy(pnum, hnum, delta, g, n_trotter=nt)
            energies_uccd[nt][i] = E_ref + E_uccd_c
        
        if (i+1) % 10 == 0:
            print(f"  Progress: {i+1}/{len(g_values)}")
    
    # Statistics
    error_ccd = np.abs(energies_ccd - energies_exact)
    error_mbpt2 = np.abs(energies_mbpt2 - energies_exact)
    error_uccd = {}
    for nt in [1, 2, 4]:
        error_uccd[nt] = np.abs(energies_uccd[nt] - energies_exact)
    
    print(f"\n{'='*80}")
    print("STATISTICS")
    print("="*80)
    print(f"\n{'Method':25s} | {'Max Error':>15s} | {'Mean Error':>15s} | {'RMS Error':>15s}")
    print("-"*80)
    print(f"{'MBPT2':25s} | {np.max(error_mbpt2):15.6e} | {np.mean(error_mbpt2):15.6e} | {np.sqrt(np.mean(error_mbpt2**2)):15.6e}")
    print(f"{'CCD':25s} | {np.max(error_ccd):15.6e} | {np.mean(error_ccd):15.6e} | {np.sqrt(np.mean(error_ccd**2)):15.6e}")
    
    for nt in [1, 2, 4]:
        label = f"UCCD (Trotter={nt})"
        print(f"{label:25s} | {np.max(error_uccd[nt]):15.6e} | {np.mean(error_uccd[nt]):15.6e} | {np.sqrt(np.mean(error_uccd[nt]**2)):15.6e}")
    
    improvement_ccd = np.mean(error_mbpt2) / np.mean(error_ccd)
    improvement_uccd = np.mean(error_mbpt2) / np.mean(error_uccd[4])
    print(f"\nCCD is {improvement_ccd:.1f}× more accurate than MBPT2")
    print(f"UCCD (Trotter=4) is {improvement_uccd:.1f}× more accurate than MBPT2")
    
    # Visualization
    print(f"\n{'='*80}")
    print("GENERATING PLOTS")
    print("="*80)
    
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # Plot 1: All energies
    ax1 = fig.add_subplot(gs[0, :2])
    ax1.plot(g_values, energies_exact, 'k-', lw=2.5, label='Exact', zorder=5)
    ax1.plot(g_values, energies_ccd, 'b--', lw=2, label='CCD', zorder=4)
    ax1.plot(g_values, energies_mbpt2, 'r:', lw=2, label='MBPT2', zorder=3)
    ax1.plot(g_values, energies_uccd[1], 'g-.', lw=1.5, alpha=0.7, label='UCCD (n=1)', zorder=2)
    ax1.plot(g_values, energies_uccd[2], 'm-.', lw=1.5, alpha=0.7, label='UCCD (n=2)', zorder=2)
    ax1.plot(g_values, energies_uccd[4], 'c-.', lw=1.5, alpha=0.7, label='UCCD (n=4)', zorder=2)
    ax1.set_xlabel('Pairing strength g', fontsize=12)
    ax1.set_ylabel('Energy', fontsize=12)
    ax1.set_title('Ground State Energy vs Pairing Strength', fontweight='bold', fontsize=13)
    ax1.legend(fontsize=10, ncol=2)
    ax1.grid(alpha=0.3)
    
    # Plot 2: Errors (log scale)
    ax2 = fig.add_subplot(gs[0, 2])
    ax2.semilogy(g_values, error_mbpt2, 'r-', lw=2, label='MBPT2')
    ax2.semilogy(g_values, error_ccd, 'b-', lw=2, label='CCD')
    ax2.semilogy(g_values, error_uccd[1], 'g--', lw=1.5, alpha=0.7, label='UCCD (n=1)')
    ax2.semilogy(g_values, error_uccd[4], 'c--', lw=1.5, alpha=0.7, label='UCCD (n=4)')
    ax2.set_xlabel('g', fontsize=11)
    ax2.set_ylabel('|E - E_exact|', fontsize=11)
    ax2.set_title('Absolute Error', fontweight='bold', fontsize=12)
    ax2.legend(fontsize=9)
    ax2.grid(alpha=0.3, which='both')
    
    # Plot 3: CCD vs UCCD comparison
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.plot(g_values, energies_exact, 'k-', lw=2, label='Exact')
    ax3.plot(g_values, energies_ccd, 'b--', lw=2, label='CCD')
    ax3.plot(g_values, energies_uccd[4], 'c-.', lw=2, label='UCCD (n=4)')
    ax3.set_xlabel('g', fontsize=11)
    ax3.set_ylabel('Energy', fontsize=11)
    ax3.set_title('CCD vs UCCD', fontweight='bold', fontsize=12)
    ax3.legend(fontsize=10)
    ax3.grid(alpha=0.3)
    
    # Plot 4: UCCD convergence with Trotter steps
    ax4 = fig.add_subplot(gs[1, 1])
    # Pick a specific g value
    g_idx = 30  # g ≈ 0.2
    trotter_range = [1, 2, 4, 8]
    uccd_energies_at_g = []
    for nt in trotter_range:
        if nt in energies_uccd:
            uccd_energies_at_g.append(energies_uccd[nt][g_idx])
        else:
            # Calculate for n=8
            evals, _, _, H_mat = exact_diagonalization_pairing(pnum, hnum, delta, g_values[g_idx])
            E_ref = H_mat[0, 0]
            E_uccd_c = calculate_uccd_energy(pnum, hnum, delta, g_values[g_idx], n_trotter=nt)
            uccd_energies_at_g.append(E_ref + E_uccd_c)
    
    ax4.plot(trotter_range, uccd_energies_at_g, 'o-', lw=2, markersize=8)
    ax4.axhline(energies_exact[g_idx], color='k', ls='--', lw=1.5, label='Exact')
    ax4.axhline(energies_ccd[g_idx], color='b', ls='--', lw=1.5, label='CCD')
    ax4.set_xlabel('Trotter steps', fontsize=11)
    ax4.set_ylabel('Energy', fontsize=11)
    ax4.set_title(f'UCCD Trotter Convergence (g={g_values[g_idx]:.2f})', fontweight='bold', fontsize=12)
    ax4.legend(fontsize=9)
    ax4.grid(alpha=0.3)
    
    # Plot 5: Energy differences
    ax5 = fig.add_subplot(gs[1, 2])
    ax5.plot(g_values, energies_exact - energies_mbpt2, 'r-', lw=2, label='Exact - MBPT2')
    ax5.plot(g_values, energies_exact - energies_ccd, 'b-', lw=2, label='Exact - CCD')
    ax5.plot(g_values, energies_exact - energies_uccd[4], 'c-', lw=2, label='Exact - UCCD(n=4)')
    ax5.axhline(0, color='k', ls='--', alpha=0.5)
    ax5.set_xlabel('g', fontsize=11)
    ax5.set_ylabel('Energy difference', fontsize=11)
    ax5.set_title('Energy Differences', fontweight='bold', fontsize=12)
    ax5.legend(fontsize=9)
    ax5.grid(alpha=0.3)
    
    # Plot 6: Improvement factors
    ax6 = fig.add_subplot(gs[2, 0])
    improvement_ccd_vs_g = error_mbpt2 / (error_ccd + 1e-12)
    improvement_uccd_vs_g = error_mbpt2 / (error_uccd[4] + 1e-12)
    ax6.plot(g_values, improvement_ccd_vs_g, 'b-', lw=2, label='CCD vs MBPT2')
    ax6.plot(g_values, improvement_uccd_vs_g, 'c-', lw=2, label='UCCD(n=4) vs MBPT2')
    ax6.axhline(1, color='k', ls='--', alpha=0.5)
    ax6.set_xlabel('g', fontsize=11)
    ax6.set_ylabel('Improvement factor', fontsize=11)
    ax6.set_title('Error Improvement over MBPT2', fontweight='bold', fontsize=12)
    ax6.legend(fontsize=9)
    ax6.grid(alpha=0.3)
    
    # Plot 7: CCD vs UCCD error ratio
    ax7 = fig.add_subplot(gs[2, 1])
    ratio_uccd_ccd = error_uccd[4] / (error_ccd + 1e-12)
    ax7.plot(g_values, ratio_uccd_ccd, 'purple', lw=2)
    ax7.axhline(1, color='k', ls='--', alpha=0.5, label='Equal accuracy')
    ax7.set_xlabel('g', fontsize=11)
    ax7.set_ylabel('Error(UCCD) / Error(CCD)', fontsize=11)
    ax7.set_title('UCCD vs CCD Accuracy', fontweight='bold', fontsize=12)
    ax7.legend(fontsize=9)
    ax7.grid(alpha=0.3)
    
    # Plot 8: All methods error comparison
    ax8 = fig.add_subplot(gs[2, 2])
    methods = ['MBPT2', 'CCD', 'UCCD\n(n=1)', 'UCCD\n(n=2)', 'UCCD\n(n=4)']
    mean_errors = [
        np.mean(error_mbpt2),
        np.mean(error_ccd),
        np.mean(error_uccd[1]),
        np.mean(error_uccd[2]),
        np.mean(error_uccd[4])
    ]
    colors = ['red', 'blue', 'lightgreen', 'mediumseagreen', 'cyan']
    ax8.bar(methods, mean_errors, color=colors, alpha=0.7, edgecolor='black')
    ax8.set_ylabel('Mean absolute error', fontsize=11)
    ax8.set_title('Average Error Comparison', fontweight='bold', fontsize=12)
    ax8.set_yscale('log')
    ax8.grid(alpha=0.3, axis='y')
    
    plt.suptitle('Pairing Model: Exact vs MBPT2 vs CCD vs UCCD', 
                 fontsize=14, fontweight='bold', y=0.995)
    plt.savefig('pairing_full_comparison.png', dpi=150, bbox_inches='tight')
    print("\n✓ Plot saved: pairing_full_comparison.png")
    plt.show()
    
    print(f"\n{'='*80}")
    print("✓ ANALYSIS COMPLETE")
    print("="*80)
