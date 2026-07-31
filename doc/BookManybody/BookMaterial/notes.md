## Notes about possible content for many-body lecture notes (covers more than FYS4480)

###  Intro chapter with basic definitions and simple examples and mathematics of many-body functions
     - Definitions of SDs etc, permutation operators,linear algebra reminder including reminder about determinants,
       vector and mtx algebra, tensor products, representations, unitary transformations, link to quantities like
     - one-body and two-body densities, rms radii etc. Discuss ansatze for wave functions and more.
     - Include in this part math on basis sets.
     - Discuss ansazte for wave functions, discuss Dicke states
###  2nd quantization for bosons and fermions and more
     - Commutation rules and definition of creation and annihilation operators
     - Proof of wick's theorem -- done in chapter3.tex, section "Proof of Wick's theorem"
     - Discuss Wick's generalized theorem -- done in chapter3.tex, statement plus a full
       proof (internal contractions vanish, two groups by induction on the appending
       lemma, then any number of groups), with the Fermi-vacuum form and numerical
       verification in wick.py
     - particle-hole picture
     - interaction, Schroedinger and Heisenberg pictures, pros and cons
     - time dependent wick's theorem
     - Gell-Man and Low's theorem
     - Adiabatic switching
     - Derivation of expressions for different parts of Hamiltonians, 1b, 2b, 3b etc
     - Add material on how to use sympy alternatively develop exercises where one writes codes for using Wick's theorem
     - Wigner-Jordan transformation and 2nd quantization -- done, chapter3.tex
       sec:jordanwigner, from week8.tex: parity strings, n_p = (1-Z_p)/2, the
       hopping and pair identities (which differ only in the sign of the YY
       channel), the eight weight-four strings of a general two-body term, and
       the O(M^4) string count.  Note: week8.tex's four-term formula for the
       pair-transfer operator is incomplete -- the correct decomposition has
       eight strings, and week8's "p != q" in the compact pairing Hamiltonian
       should read "p < q".  Both corrected and verified in jordanwigner.py
     - Chapter 3 now has an exercise section (six exercises with answers,
       covering normal ordering, Wick, the generalised theorem, the Fermi
       vacuum and Jordan-Wigner)
     - Baker-Campbell-Hausdorf (BCH) -- now developed in chapter6.tex
     - Suzuki-Trotter as an approximation to BCH -- now in chapter6.tex
### Physical systems and models  (chapter4.tex, done)
     - Lipkin (quasispin SU(2)), pairing, Fermi-Hubbard, Calogero-Sutherland
     - Pairing plus a particle-hole term, i.e. the chapter-7 model introduced
       here rather than arriving unannounced in chapter 7.  Includes the TikZ
       level scheme fig:4-levels with the reference determinant, the
       seniority-conserving 2p-2h pair transfer, and the 1p-1h and
       pair-breaking 2p-2h excitations that only V_ph reaches.  Key point:
       the q = r terms of V_ph are just the pairing term again, so g -> g + 2f
       in every seniority-conserving matrix element (verified in models.py)
     - Quantum Heisenberg model, as the strong-coupling limit of Hubbard and
       in its own right: 4-site ring E_0 = -2J checked against the Hubbard
       table from both sides, Bethe ansatz J(1/4 - ln 2) with finite-size
       scaling to N = 16, Jordan-Wigner back to spinless fermions (spectra
       verified identical, not just ground states)
     - The models on a quantum computer (sec:qubitmodels): the general
       two-body Hamiltonian as Pauli strings, the pairing Hamiltonian in the
       full 2L-qubit encoding and in the compact one-qubit-per-level
       (quasispin) encoding, and the fact that the particle-hole term makes
       the compact encoding illegal -- leakage out of the seniority-zero
       space is exactly f/2.  String counts in tab:4-paulicounts
     - Still to add if wanted: two-dimensional Heisenberg / Hubbard, the
       XXZ anisotropy, spin gap and the Haldane conjecture; Bravyi-Kitaev and
       parity encodings worked out for the pairing model

### FCI and diagrams and particle-hole representations
     - Basics of FCI
     - Rewriting in terms of a particle-hole picture
     - Discuss slater determinants and similarity transformations and algorithms for solving eigenvalue problems
     - May discuss eigenvector continuation
     - Introduce a diagrammatic representation
### Mean-field theories  (chapter6.tex, done)
     - Hartree-Fock in coordinate space and 2nd quantization
     - Variational calculus and Lagrange multipliers
     - The density matrix and the self-consistent field algorithm
     - Thouless theorem
     - Baker-Campbell-Hausdorff and Suzuki-Trotter (moved here from the
       2nd-quantization chapter, since Thouless' theorem needs BCH and the
       same algebra returns in coupled cluster and in quantum simulation)
     - Slater dets in HF theory
     - Stability of the HF solution; the stability matrix is the RPA matrix.
       Worked instability: the Lipkin model at chi = 1.  Worked non-example:
       the pairing model, where HF gains nothing.
     - DFT links
     - The electron gas as example: closed-form HF, band width, the vanishing
       effective mass, E_0/N = 2.21/r_s^2 - 0.916/r_s and the bound minimum
     - FCI and HF, diagrammatic representations and critical discussions

### Mean-field applications: TDA and RPA  (chapter7.tex, done)
     - The pairing + particle-hole model (Hjorth-Jensen, Dean, Hagen, Kvaal,
       J. Phys. G 37, 064035 (2010), Sec. 2.1); reduces to the chapter-4
       pairing model at f = 0
     - Exact diagonalization in the S_z = 0 sector
     - The equations of motion; TDA as CIS; RPA and the quasiboson approximation
     - RPA matrix = chapter-6 stability matrix in a different arrangement:
       real frequencies <=> stable mean field.  Soft mode at the Lipkin
       transition.
     - The RPA correlation energy and the RPA ground state as exp(Z_2)|HF>
     - BCS, the pairing gap and the (spurious) finite-system transition
     - QRPA: the Goldstone mode from broken U(1), identified by its overlap
       with the number operator, and the Delta N = 0 vs pairing-vibration
       mode content
     - Still to add if wanted: number projection, renormalised QRPA, HFB

### Density functional theory  (chapter8.tex, done)
     - Born-Oppenheimer and the electronic Hamiltonian
     - One- and two-body reduced density matrices in coordinate space; the
       Slater-determinant factorisation of rho^(2) and the N-representability
       problem
     - Hohenberg-Kohn theorems I and II with the reductio ad absurdum proof;
       v-representability and the Levy-Lieb constrained search
     - Thomas-Fermi from the electron gas of chapter 6, and why it fails
     - The Kohn-Sham construction, the total-energy formula, KS vs HF
     - LDA built explicitly from the uniform gas (analytic in 3D, numerical
       for the 1D model interaction); the exchange hole and its sum rule
     - The self-interaction error and the locality error, both measured on
       the trap system of chapters 2 and 6
     - Still to add if wanted: GGA construction from exact constraints,
       time-dependent DFT as the density-functional analogue of RPA,
       DFT for nuclei (Skyrme functionals)
###  Many-body perturbation theory  (chapter9.tex, done -- algebraic only)
     - The exact result Delta E = <Phi_0|H_I|Psi_0>, P/Q projections, the
       resolvent expansion with a free parameter omega
     - Brillouin-Wigner (omega = E), its exact resummation, and its failure
       of size extensivity
     - Rayleigh-Schroedinger (omega = W_0) to first, second and third order,
       with the renormalisation term; the recursion to arbitrary order
     - The wave operator and the link to FCI: C_ij^ab to first order
     - All four models to third order (pairing, Lipkin, Hubbard, pairing+ph)
       with comparison against HF, CID, RPA and FCI
     - Convergence of the series, and the Moller-Plesset vs Epstein-Nesbet
       partition
     - Still to add: the diagrammatic representation and the linked-diagram
       theorem from time-dependent theory (deliberately left out for now)
###  Coupled cluster theories, standard and unitary  (chapter10.tex, done)
     - The exponential ansatz, C vs T amplitudes, and why the exponential
       cures the size-consistency failure of CID
     - The similarity-transformed Hamiltonian; non-hermiticity as the price
       of a terminating BCH expansion (exactly at the fourfold commutator)
     - Full algebraic derivation of the CCD equations (three ladder/ring
       linear terms, four quadratic topologies) and of the CCSD singles and
       doubles equations with the Stanton-Gauss-Watts-Bartlett intermediates
     - CC as a resummation of MBPT: the order-by-order expansion reproduces
       MP2, MP3 and (for the pairing model) MP4 exactly
     - Pairing model worked example with the CCD code from week48.ipynb;
       CCSD on the pairing+particle-hole model where singles are non-zero
     - Size extensivity verified to machine precision for 2 and 3 subsystems
     - Comparison against HF, MBPT2/3, CID, RPA and FCI
     - Unitary coupled cluster (UCCD and UCCSD): the variational property,
       why the BCH series no longer terminates, the pair generators for the
       pairing model, Trotterisation linked to chapter 6, and the two key
       numerical facts -- re-optimisation absorbs the splitting error, and
       the energy error falls as 1/n^2 while the state error falls as 1/n
     - The UnitaryCC class is the VQE ansatz, reused in chapter11.tex
     - Still to add: equation-of-motion CC, the diagrammatic derivation,
       CCSD(T)
###  Applications: quantum dots  (chapter11.tex, done)
     - Two-dimensional parabolic trap with the bare Coulomb interaction, from
       qdhfmbptccsd_uccsd.ipynb sections 1-3 (Hamiltonian, polar basis,
       Coulomb matrix element) plus the HF/MP2/CC results
     - Coulomb integrals in closed form (Anisimova-Matulis); checks:
       <00;00|v|00;00> = sqrt(pi/2) sqrt(hw) to 12 digits, the sqrt(hw)
       scaling, m_p+m_q = m_r+m_s, and the particle-swap/hermiticity
       symmetries
     - HF, MP2, CCD, CCSD reproduce the notebook's HF and MP2 numbers exactly
       (E_HF = 3.16192140 and 20.72025707; E_MP2 = -0.13488329 and
       -0.41769581 at 42 orbitals)
     - IMPORTANT, two errors in the source notebook, both corrected here:
       (i) it claims Brillouin's theorem makes T1 = 0 so that CCSD = CCD.
       False: the CCSD singles equation is driven by the doubles through the
       <ma||ef> and <nm||ei> terms, which survive at t1 = 0.  Here
       max|t1| = 5.7e-3 (N=2) and 1.0e-2 (N=6).
       (ii) its CC energies are wrong (E_CCD = -0.09725, E_CCSD = -0.05858 for
       N=2, explained away as a "DIIS-path artefact").  The correct values are
       -0.14799908 and -0.14829527, verified because for two electrons CCSD
       must equal FCI, and it does to 3e-13 in every basis
     - Benchmarks: two-electron FCI in every basis (dim 861 at 42 orbitals),
       Taut's analytic E = 3 at hw = 1, and a full six-electron FCI in the
       12-orbital basis (dim 924) where MP2/CCD/CCSD/UCCSD recover
       92.55/98.13/98.47/99.20 % of the correlation energy
     - UCCSD beats CCSD even when evaluated at the borrowed CCSD amplitudes,
       and stays above FCI as the variational principle requires.  Note the
       sign convention: UnitaryCC's doubles generator is a+_a a+_b a_i a_j,
       so t2 must be imported with a minus (verified against a full
       variational optimisation, cosine 1.0000)
     - Still to add if wanted: perturbative triples, open shells (N=4),
       twelve electrons, and the Monte Carlo comparison promised in the text

###  Green's function theory and parquet theory
     - Notes ready but may not teach
###  SRG and IMSRG
     - Notes ready but may not teach
###  Monte Carlo methods
     - Notes for FYS4411
###  Quantum computing
     - VQE and unitary CC, notes ready
###  time-dependent many-body theory
     - Notes to be developed
###  Applications to different systems like the electron gass, Lipkin model, Pairing model, infinite nuclear matter, and more
     - Notes ready, some exercises and projects to be developed.