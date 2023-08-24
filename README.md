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

## Practicalities

1. Four lectures per week, Fall semester, 10 ECTS. Thursday 1015am-12pm and Friday 1015-12pm. The lectures will be recorded and linked to this site;
2. Two hours of exercise sessions for work on projects and exercises, Friday 1215pm-2pm;
3. Two projects which are graded and count 30% each of the final grade and a final oral exam that counts 40% of the final grade;
4. Weekly assignments which are not graded;
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
- FCI, diagrams and particle-hole representations
     - Basics of FCI
     - Rewriting in terms of a particle-hole picture
     - Slater determinants and similarity transformations and algorithms for solving eigenvalue problems
     - Eigenvector continuation
     - Introduce a diagrammatic representation
- Mean-field theories
     - Hartree-Fock in coordinate space and 2nd quantization
     - Thouless theorem
     - Slater dets in HF theory
     - Density functional theory
     - The electron gas as example
     - FCI and HF, diagrammatic representations and critical discussions
- Many-body perturbation theory
     - Time dependent and time-independent representations
     - Brillouin-Wigner and Rayleigh-Schrødinger pert theory
     - Diagrammatic representation
     - Linked-diagram theorem based on time-dependent theory
- Coupled cluster theories, standard and unitary
     - Derivation of equations for singles and doubles, reminder on unitary transformations
     - Unitary coupled cluster theory
- Green's function theory and parquet theory
- SRG and IMSRG
- Monte Carlo methods
     - Taught in FYS4411
- Quantum computing
     - VQE and unitary CC
- Time-dependent many-body theory
- Applications to different systems like the electron gass, Lipkin model, Pairing model, infinite nuclear matter, and more




## Teaching schedule with links to material

###  Week 34, August 21-25

- Topics to be covered
  - Thursday: Introduction to many-body physics, notations and definitions
    - Video of lecture at https://youtu.be/uKku_uAnmRM
  - Friday: Discussion of notations and important properties in many-body physics
- Lecture Material at https://manybodyphysics.github.io/FYS4480/doc/pub/secondquant/html/secondquant-bs.html

###  Week 35, August 28-September 1

- Topics to be covered
  - Thursday: Fermion state functions and computation of expectation values in first quantization
  - Friday: Introduction of second quantization
- Lecture Material at https://manybodyphysics.github.io/FYS4480/doc/pub/secondquant/html/secondquant-bs.html
- First exercise set at https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/Exercises/2023/ExercisesWeek35.pdf

### Week 36, September 4-8

- Topics to be covered
  - Thursday: Second quantization, operators in second quantization and diagrammatic representation
  - Friday: Second quantization and Wick's theorem
- Lecture Material: These slides, handwritten notes and	Szabo and Ostlund sections 2.3 and 2.4. Sections 3.1-3.3 of Shavitt and Bartlett covers most of the material discussed this week.
- Second exercise set at https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/Exercises/2023/ExercisesWeek36.pdf



### Week 37, September 11-15 

- Topics to be covered
  - Thursday: Exercises for week 37
  - Friday: Wrappping up Wick's theorem, Wick's generalized theorem and diagrammatic representation, Particle-hole formalism and definition of new reference state, normalordering of operators
- Lecture Material: These slides, handwritten notes and chapter 3 of Shavitt and Bartlett covers most of the material discussed this week.
- Third exercise set at https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/Exercises/2023/ExercisesWeek37.pdf

### Week 38, September 18-22

- Topics to be covered
  - Thursday: Discussion of particle-hole formalism with examples
  - Friday: Particle-hole formalism and definition of new reference state, normalordering of operators
- Lecture Material: These slides, handwritten notes and chapter 3 and 4 of Shavitt and Bartlett covers most of the material discussed this week.
- Fourth exercise set at https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/Exercises/2023/ExercisesWeek38.pdf

###  Week 39, September 25-29
- Topics to be covered
  - Thursday:
    - Repetition  of particle-hole formalism
    - Diagrammatic representation
    - Introduction of full configuration interaction theory
  - Friday: 
    - Full configuration interaction (FCI) theory
    - Lipkin model as an example of applications of FCI theory
- Lecture Material: These slides, handwritten notes
- Fifth exercise set at https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/Exercises/2023/ExercisesWeek39.pdf

### Week 40, October 2-6
- Topics to be covered
  - Thursday:
    - Repetition  of Full Configuration Interaction theory
    - Start Hartree-Fock theory
  - Friday: 
    - Hartree-Fock theory and stability of equations
- Lecture Material: These slides, handwritten notes and Szabo and Ostlund, sections 3.1-3.4
- Sixth exercise set at https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/Exercises/2023/ExercisesWeek40.pdf

### Week 41, October 9-13, 2023
- Topics to be covered
  - Thursday:
    - Repetition  of derivation of Hartree-Fock theory
    - Stability of Hartree-Fock theory and Thouless theorem
  - Friday:
    - Work on first midterm
- Lecture Material: Lecture slides week 41, handwritten notes and Szabo and Ostlund, chapter 3
- First midterm set at https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/Exercises/2023/FirstMidterm2023.pdf


### Week 42, October 16-20, 2023
- Topics to be covered
  - Thursday:
    - Repetition  of derivation of Hartree-Fock theory
    - Discussion of implementation of Hartree-Fock code
    - Stability of Hartree-Fock theory and Thouless theorem
  o Friday:
    * Work on first midterm
- Lecture Material: These slides, handwritten notes and Szabo and Ostlund, chapter 3
- First midterm set at https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/Exercises/2023/FirstMidterm2023.pdf




###  Week 43, October 23-27, 2023
- Topics to be covered
  - Thursday:
    - Linking FCI and Hartree-Fock theory with many-body perturbation theory
    - Time-independent perturbation theory
  - Friday: 
    - Time-dependent perturbation theory and diagrammatic representation
- Lecture Material: These slides, handwritten notes and Shavitt and Bartlett chapter 4 and sections 5.1-5.3
- Seventh exercise set at https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/Exercises/2023/ExercisesWeek43.pdf



###  Week 44, October 30-November 2, 2023
- Topics to be covered
  - Thursday:
    - Time-independent perturbation theory and diagrammatic representation
    - Examples of contributions to perturbation theory
  - Friday: 
    - Time-dependent perturbation theory and diagrammatic representation
    - Discussion of the Heisenberg, Schroedinger and Interaction pictures
- Lecture Material: These slides, handwritten notes and Shavitt and Bartlett chapter 4-6
- Eight exercise set at https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/Exercises/2023/ExercisesWeek44.pdf

###  Week 45, November 6-10, 2023
- Topics to be covered
  - Friday:
    - Adiabatic hypothesis and linked diagram theorem
    - Examples of diagrams and discussion of Pauli violating diagrams 
    - Summary of perturbation theory and start discussion of Coupled Cluster theory
- Lecture Material: Slides at https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/src/week45/LatexSlides/slidesweek45.pdf, handwritten notes and Shavitt and Bartlett chapters 4-6
- Ninth exercise set at https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/Exercises/2023/ExercisesWeek45.pdf

###  Week 46, November 13-17, 2023
- Topics to be covered
  - Thursday:
    - Summary many-body perturbation theory and discussion of Coupled Cluster theory
  - Friday:
    - Coupled Cluster theory with doubles excitations
    - Diagrammatic representation of Coupled Cluster theory
- Lecture Material: Lecture slides at https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/src/week46/LatexSlides/slidesweek46.pdf, handwritten notes and Shavitt and Bartlett chapter 9
- Work on second midterm set at https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/Exercises/2023/SecondMidterm2023.pdf

###  Week 47, November 20-24, 2023
- Topics to be covered
  - Thursday:
    - Coupled Cluster theory, singles and doubles excitations
  - Friday:
    - Coupled Cluster theory and summary of course
- Lecture Material: Slides at https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/src/week47/LatexSlides/slidesweek47.pdf, handwritten notes and Shavitt and Bartlett chapters 9 and 10
- Work on second midterm set at https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/Exercises/2023/SecondMidterm2023.pdf

###  Week 48, November 27-December 1, 2023
- Topics to be covered
  - Thursday:
    - Coupled Cluster theory, singles and doubles excitations
  - Friday:
    - Coupled Cluster theory and summary of course
- Lecture Material: Slides at https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/src/week47/LatexSlides/slidesweek47.pdf, handwritten notes and Shavitt and Bartlett chapters 9 and 10
- Work on second midterm set at https://github.com/ManyBodyPhysics/FYS4480/blob/master/doc/Exercises/2023/SecondMidterm2023.pdf


###  Week 50, December 19 and 21, 2023
- Final oral exam, day to be determined

