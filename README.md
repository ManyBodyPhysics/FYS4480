# Quantum Mechanics for Many-Particle systems

This course gives an introduction to the quantum mechanics of
many-body systems and the methods relevant for many-body
problems in such diverse areas as atomic, molecular, solid-state and
nuclear physics, chemistry and materials science. A theoretical
understanding of the behavior of quantum-mechanical many-body systems, that is, systems containing many interacting particles - is a
considerable challenge in that, normally, no exact solution can be found.
Instead, reliable methods are needed for approximate but accurate
simulations of such systems.

The aim of this course is to present some of the most widely used
many-body methods, starting with the underlying formalism of second
quantization. The topics
covered include second quantization with creation and annivilation operators, Wick's theorem, Feynman diagram rules, microscopic mean-field theories
(Hartree-Fock and Kohn-Sham theories), many-body perturbation theory,
large-scale diagonalization methods, coupled cluster theory,
algorithms from quantum computing, and
Green's function approaches. Both fermionic and bosonic systems are
discussed, depending on the interests of the participants.
Selected physical systems from various fields such as quantum
chemistry, solid-state physics and nuclear physics are studied,
depending on the background and interests of the participants.


## Instructor information
* _Name_: Morten Hjorth-Jensen
* _Email_: morten.hjorth-jensen@fys.uio.no
* _Phone_: +47-48257387
* _Office_: Department of Physics, University of Oslo, Eastern wing, room FØ470 
* _Name_: Øyvind Sigmundson Schøyen
* _Email_: oyvinssc@student.matnat.uio.no
* _Office_: Department of Physics, University of Oslo, Eastern wing, room FØ452 

## Practicalities

1. Four lectures per week, Fall semester, 10 ECTS. The lectures will be recorded and linked to this site;
2. Two hours of exercise sessions for work on projects and exercises;
3. Two projects which are graded and count 30% each of the final grade and a final oral exam that counts 40% of the final grade;
4. Weekly assignments;
5. The course is offered as a so-called _cloned_ course,  FYS4480 at the Master of Science level and FYS9480 as a PhD course;
6. Weekly email with summary of activities will be mailed to all participants;

## Textbooks

_Recommended textbooks_:
In the folder https://github.com/ManyBodyPhysics/FYS4480/tree/master/doc/Literature you will find the textbooks we will be following.
Weekly reading assignments based on these texts will be sent before each week.
In particular we recommend the texts by
- Negele and Orland at https://www.taylorfrancis.com/books/mono/10.1201/9780429497926/quantum-many-particle-systems-john-negele-henri-orland
- Szabo and Ostlund at for example https://www.amazon.com/Modern-Quantum-Chemistry-Introduction-Electronic/dp/0486691861?asin=0486691861&revisionId=&format=4&depth=1
- Shavitt and Bartlett at https://www.cambridge.org/core/books/manybody-methods-in-chemistry-and-physics/D12027E4DAF75CE8214671D842C6B80C

## Topics (not all will be discussed)

- Intro chapter with basic definitions and simple examples and mathematics of many-body functions
     - Definitions of SDs etc, permutation operators,linear algebra reminder including reminder about determinants,
       vector and mtx algebra, tensor products, representations, unitary transformations, link to quantities like
     - one-body and two-body densities, rms radii etc. Discuss ansatze for wave functions and more.
     - Ansaztes for wave functions
- 2nd quantization for bosons and fermions and more
     - Commutation rules and definition of creation and annihilation operators
     - Proof of wick's theorem
     - Discuss Wick's generalized theorem
     - particle-hole picture
     - interaction, Schroedinger and Heisenberg pictures, pros and cons
     - time dependent wick's theorem
     - Gell-Man and Low's theorem
     - Adiabatic switching
     - Derivation of expressions for different parts of Hamiltonians, 1b, 2b, 3b etc
     - Wigner-Jordan transformation and 2nd quantization
     - Baker-Campbell-Hausdorf (BCH)
     - Suzuki-Trotter as an approximation to BCH
- FCI and diagrams and particle-hole representations
     - Basics of FCI
     - Rewriting in terms of a particle-hole picture
     - Discuss slater determinants and similarity transformations and algorithms for solving eigenvalue problems
     - Eigenvector continuation
     - Introduce a diagrammatic representation
- Mean-field theories
     - Hartree-Fock in coordinate space and 2nd quantization
     - Thouless theorem
     - Slater dets in HF theory
     - DFT links
     - The electron gas as example
     - FCI and HF, diagrammatic representations and critical discussions
- Many-body perturbation theory
     - Time dependent and time-independent representation
     - Brillouin-Wigner and Rayleigh-Schrødinger pert theory
     - Diagrammatic representation
     - Linked-diagram theorem based on time-dependent theory
- Coupled cluster theories, standard and unitary
     - Derivation of equations for singles and doubles, reminder on unitary transformations
     - non-hermiticity
     - Specialize to CCD case and compare with FCI and MBPT
- Green's function theory and parquet theory
- SRG and IMSRG
- Monte Carlo methods
     - Taught in FYS4411
- Quantum computing
     - VQE and unitary CC
- Time-dependent many-body theory
- Applications to different systems like the electron gass, Lipkin model, Pairing model, infinite nuclear matter, and more



## Teaching schedule with links to material

###  Week 34, August 22-26, 2022

- Topics to be covered
  - Thursday: Introduction to many-body physics, notations and definitions
  - Video of lecture at https://www.youtube.com/watch?v=ScHgSRCmq-0
  - Friday: Discussion of notations and important properties in many-body physics
  - Video of lecture at https://www.youtube.com/watch?v=M14xN2BAULg
- Lecture Material at https://manybodyphysics.github.io/FYS4480/doc/pub/secondquant/html/secondquant-bs.html
- Recommended reading: Szabo and Ostlund, chapters 1 and 2, see https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/Literature/szaboostlund.pdf
###  Week 35, August 29-September 2, 2022

- Topics to be covered
  - Thursday: Fermion state functions and computation of expectation values in first quantization
  - Video of lecture at https://youtu.be/-AkzRdrq3Qc
  - Friday: Introduction of second quantization
  - Video of lecture at https://youtu.be/Nbzp3wr0vyA
- Lecture Material at https://manybodyphysics.github.io/FYS4480/doc/pub/secondquant/html/secondquant-bs.html
- Recommended reading: Szabo and Ostlund, chapters 1 and 2, see https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/Literature/szaboostlund.pdf
- First exercise set at https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/Exercises/2022/ExercisesWeek35.pdf

### Week 36, September 5-9

- Topics to be covered
  - Thursday: Second quantization, operators in second quantization and diagrammatic representation
  - Video of lecture at https://youtu.be/ofow8wjvP6U
  - Friday: Second quantization and Wick's theorem
  - Video of lecture at https://youtu.be/qylybY6N5rQ
- Lecture Material: These slides, handwritten notes and	Szabo and Ostlund sections 2.3 and 2.4. Sections 3.1-3.3 of Shavitt and Bartlett covers most of the material discussed this week.
- Second exercise set at https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/Exercises/2022/ExercisesWeek36.pdf




### Week 37, September 12-16

- Topics to be covered
  - Thursday: Exercises for week 37
  - Friday: Wrappping up Wick's theorem, Wick's generalized theorem and diagrammatic representation, Particle-hole formalism and definition of new reference state, normalordering of operators
  - Video of lecture at https://youtu.be/LzLS_C7GMEM
- Lecture Material: These slides, handwritten notes and chapter 3 of Shavitt and Bartlett covers most of the material discussed this week.
- Third exercise set at https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/Exercises/2022/ExercisesWeek37.pdf


### Week 38, September 19-23

- Topics to be covered
  - Thursday: Discussion of particle-hole formalism with examples
  - Video of lecture at https://youtu.be/zfaCSdyJlis
  - Friday: Particle-hole formalism and definition of new reference state, normalordering of operators
  - Video of lecture at	https://youtu.be/NfkG8isV4wc
- Lecture Material: These slides, handwritten notes and chapter 3 and 4 of Shavitt and Bartlett covers most of the material discussed this week.
- Fourth exercise set at https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/Exercises/2022/ExercisesWeek38.pdf


###  Week 39, September 26-30
- Topics to be covered
  - Thursday:
    - Repetition  of particle-hole formalism
    - Diagrammatic representation
    - Introduction of full configuration interaction theory
    - Video of Lecture at https://youtu.be/_VtGqoXK6b8
  - Friday: 
    - Full configuration interaction (FCI) theory
    - Lipkin model as an example of applications of FCI theory
    - Video of lecture https://youtu.be/6goQRQv62gY
- Lecture Material: These slides, handwritten notes
- Fifth exercise set at https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/Exercises/2022/ExercisesWeek39.pdf
