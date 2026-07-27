## Notes about possible content for many-body lecture notes (covers more than FYS4480)

###  Intro chapter with basic definitions and simple examples and mathematics of many-body functions
     - Definitions of SDs etc, permutation operators,linear algebra reminder including reminder about determinants,
       vector and mtx algebra, tensor products, representations, unitary transformations, link to quantities like
     - one-body and two-body densities, rms radii etc. Discuss ansatze for wave functions and more.
     - Include in this part math on basis sets.
     - Discuss ansazte for wave functions, discuss Dicke states
###  2nd quantization for bosons and fermions and more
     - Commutation rules and definition of creation and annihilation operators
     - Proof of wick's theorem
     - Discuss Wick's generalized theorem
     - particle-hole picture
     - interaction, Schroedinger and Heisenberg pictures, pros and cons
     - time dependent wick's theorem
     - Gell-Man and Low's theorem
     - Adiabatic switching
     - Derivation of expressions for different parts of Hamiltonians, 1b, 2b, 3b etc
     - Add material on how to use sympy alternatively develop exercises where one writes codes for using Wick's theorem
     - Wigner-Jordan transformation and 2nd quantization
     - Baker-Campbell-Hausdorf (BCH) -- now developed in chapter6.tex
     - Suzuki-Trotter as an approximation to BCH -- now in chapter6.tex
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
###  Many-body perturbation theory
     - Time dependent and time-independent representation
     - Brillouin-Wigner and Rayleigh-Schrødinger pert theory
     - Diagrammatic representation
     - Linked-diagram theorem based on time-dependent theory
###  Coupled cluster theories, standard and unitary
     - Derivation of equations for singles and doubles, reminder on unitary transformations
     - non-hermiticity
     - Specialize to CCD case and compare with FCI and MBPT
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