# doc/BookManybody — the textbook

Source for *Quantum mechanics for Many-particle Systems: from standard methods
to quantum computing and machine learning* (Morten Hjorth-Jensen), a Springer
monograph. The same material is also maintained as a Jupyter-book in
`../LectureNotes/` and drives the weekly FYS4480/9480 course.

## Files

```
book.tex           root document — svmono class, all macros, \part/\input structure
chapter1.tex       Linear algebra: notation, then numerical linear algebra,
                   eigenvalue problems and the SVD
chapter2.tex       From linear algebra to many-body physics: Slater
                   determinants, the energy functional, Monte Carlo updates
chapter3.tex       Second quantization, the particle-hole formalism and
                   Wick's theorem, with proofs of both the ordinary and the
                   generalised theorem, the bit representation, the
                   Jordan-Wigner transformation, and an exercise section
chapter4.tex       Physical systems and models: Lipkin, pairing, pairing plus
                   a particle-hole term (the model of chapter 7, with a TikZ
                   level scheme showing 1p-1h and 2p-2h excitations),
                   Fermi-Hubbard, Heisenberg, Calogero-Sutherland, and the
                   qubit encodings of all of them (sec:qubitmodels)
chapter5.tex       Full configuration interaction: the CI expansion, the
                   eigenvalue problem, the exponential wall, truncations
chapter6.tex       Mean-field approaches: Hartree-Fock, Thouless' theorem,
                   BCH and Trotter, stability, the infinite electron gas
chapter7.tex       Mean-field applications: the pairing + particle-hole model,
                   TDA, RPA, BCS and QRPA
chapter8.tex       Density functional theory: density matrices, the
                   Hohenberg-Kohn theorems, Kohn-Sham, the LDA
chapter9.tex       Many-body perturbation theory: Brillouin-Wigner and
                   Rayleigh-Schroedinger to third order, four models
chapter10.tex      Coupled cluster theory: the exponential ansatz, full CCD
                   and CCSD derivations, the pairing model, and unitary
                   coupled cluster with Trotterisation (the VQE ansatz)
chapter11.tex      Applications: the two-dimensional quantum dot. Coulomb
                   integrals in closed form, then HF, MP2, CCD, CCSD and
                   UCCSD against exact diagonalisation. Source material:
                   qdhfmbptccsd_uccsd.ipynb
add.tex  linalg.tex  eig.tex  vmcfermion2015.do.txt  version2024.do.txt
calogero.tex  SecondMidterm2025.tex  Hubbard.ipynb  combinatorial.tex
week40.do.txt  pmm.tex  dft.do.txt  mbpt.ipynb  cctheory.tex  week48.ipynb
pairing1.ipynb  pairing_ph_exact_diag.ipynb  pairing_ph_rpa_qrpa.ipynb
FYS4480*.pdf
                   staging files and handwritten notes, merged into the
                   chapters above; not included by book.tex
variationalmc.tex  Variational Monte Carlo
Greensfunction.tex Green's function theory
preface.tex  acknow.tex  dedic.tex  acronym.tex
svmono.cls  spphys.bst  spbasic.bst  spmpsci.bst   (Springer class + bib styles)
Exercises/         weekly exercise sheets and answers (see below)
BookMaterial/      legacy DocOnce sources (*.do.txt), notes.md outline, Programs/
```

`BookMaterial/notes.md` is the **planned outline** for the full book — consult
it before proposing new chapter content, and update it when the plan changes.

## Building

`book.tex` compiles with `pdflatex` + `bibtex` + `makeindex` (style `svind.ist`).
It uses `[dvips]` options on `epsfig`/`graphicx` and loads `pstricks`, so
figure handling is DVI-oriented in places.

The build is clean: no errors, no undefined references or citations.

### Things that are deliberate, not bugs

1. The front matter (`dedic`, `preface`, `acknow`, `acronym`) is commented out
   in `book.tex`, which is why the PDF is shorter than the content suggests.
2. Most `\part`s in `book.tex` have their `\input` lines commented out —
   work in progress.
3. `book.tex` loads `pstricks` and `[dvips]{epsfig}`, though nothing in
   chapters 1–10 uses them. Dropping them would allow a plain `pdflatex` route.

## LaTeX conventions

- Chapter files start directly with `\chapter{...}`, no preamble (see
  `chapter2.tex`). All packages and macros belong in `book.tex`.
- **Never redefine a macro locally that `book.tex` already defines.** The book
  preamble defines, among others:
  - `\be \ee \bea \eea \beN \eeN \beaN \eeaN \bdm \edm` (equation shorthands)
  - `\op{} \Op{} \OP{} \vec{} \matr{} \bfv{} \uvec{} \det{}`
  - `\mean{}` and sized variants `\meanb \meanbb \meanbbb \meanbbbb`
  - `\brab \ketb \overlap \bracket \projection` + `b/bb/bbb/bbbb` sized variants
  - `\ud{} \udd{} \bigO \real{} \prob \cov \var \PsiT \Det{}`
  - `Python` and `C++` `lstnewenvironment`s for code listings
- Equations: numbered `equation` environments with labels in the pattern
  `\label{eq:2-1a}` (chapter number, running index, letter suffix). Reference
  as `Eq.~(\ref{eq:2-1a})`.
- All ten chapters use numbered `\section{...}`. Chapters 1--3 end with an
  `Exercises` section containing warm-up exercises followed by the weekly
  FYS4480 sessions and their answers; chapters 4--10 relegate the long
  derivations to their exercises, each with a worked solution written as
  `\subsection*{Exercise N: title}` followed by `\paragraph{Answer to a).}`.
- `\Det{...}` expands to `{|\bfv{#1}|}`, so it takes a *matrix symbol*, not an
  arbitrary expression — `\Det{\bm{A}}` is fine, `\Det{\sum_k c_k b_{jk}}`
  breaks. Use `\det(...)` for expressions. `\bra`/`\ket` are not defined;
  use `\bbra`/`\bket`, `\braket{}{}` and `\bracket{}{}{}`.
- svmono provides `theorem`, used in chapter 6 for Thouless' theorem.
- `book.tex` also defines `\wmark`, `\wline`, `\wlineb` and `\wstrut` for the
  Wick contraction lines of Chapter 3, `\levelpanel` plus the `pairarrow` and
  `pharrow` TikZ styles for the single-particle level schemes of Chapter 4,
  and a `notebox` environment for the many-body and quantum-computing
  connection boxes.
- `\levelpanel{tag}{x offset}{occupied slots}{caption}` draws one panel of
  four doubly degenerate levels and leaves the coordinates `tag-p-u` and
  `tag-p-d` behind, so excitation arrows are drawn *after* all the panels, in
  the same `tikzpicture`. See `fig:4-levels`.
- Add `\index{...}` entries for new terminology — the book has `\makeindex`.
- Citations use `\cite{}` with `\bibliographystyle{unsrt}`; Springer styles
  (`spphys`, `spbasic`, `spmpsci`) are available if the publisher requires them.
- Notation follows `chapter2.tex`: `a^\dagger_\alpha`, `|0\rangle`,
  `|\alpha_1\dots\alpha_n\rangle_{\mathrm{AS}}`, `\{A,B\}` for anticommutators.

## Exercises/

Weekly sheets for the course. Two document styles coexist:

- **Exercise sheets** (`ExercisesWeekNN.tex`, `FirstMidterm2025.tex`,
  `SecondMidterm2025.tex`): `\documentclass[prc]{revtex4}`, title
  `Exercises FYS4480/9480, week NN, <Month D-D, YYYY>`, body organised as
  `\subsection*{Exercise 1}`, `\subsection*{Exercise 2}`, …
- **Answer sheets** (`AnswersWeekNN.tex`, `AnswerLipkin.tex`):
  `\documentclass[a4paper, 11pt, notitlepage, english]{article}` with a fuller
  preamble including `listings` set up for Python.

Both define their own local macros (`\bra`, `\ket`, `\element`, `\normord`,
`\Heff`, `\Veff`, `\Span`, `\tr`, `\diag`, …) — copy the preamble from the
nearest existing week rather than inventing one. Each file is standalone.

Recurring model systems across the weeks: the **Lipkin model**, the **pairing
model**, the **electron gas** (`electrongas.tex`), and simple harmonic-oscillator
systems. `week40.do.txt` (Hartree-Fock, Thouless, stability, the electron gas)
was merged into `chapter6.tex`, together with the BCH/Trotter material from
`pmm.tex`. The two notebooks `pairing_ph_exact_diag.ipynb` and
`pairing_ph_rpa_qrpa.ipynb` were merged into `chapter7.tex`, and
`dft.do.txt` into `chapter8.tex`, and `mbpt.ipynb` into `chapter9.tex`.
`cctheory.tex`, `week48.ipynb` and `pairing1.ipynb` were merged into
`chapter10.tex`. Chapters
9 and 10 are deliberately algebraic: the diagrammatic representation has not
been written yet. Note that the explicit CCSD equations in `cctheory.tex`
contain several index errors (misplaced amplitude labels, a spurious
`t_e^e`); chapter 10 uses the standard correct forms, validated numerically
against FCI. The `UnitaryCC` class of `coupledcluster.py` is written to be
reused by `chapter11.tex` and the quantum-computing chapters: it builds the
many-body Hamiltonian and the anti-Hermitian generators in the determinant
basis, and offers both the exact exponential and the Trotterised product
form.

## Jupyter-book target

`../LectureNotes/` holds the executable version: `_config.yml`, `_toc.yml`,
markdown front matter (`intro.md`, `teachers.md`, `textbooks.md`) and chapter
notebooks. `linearalgebra.ipynb` is the companion to `chapter1.tex` and
`manybodybasics.ipynb` to `chapter2.tex`, `wicktheorem.ipynb` to the Wick
sections of `chapter3.tex`, `jordanwigner.ipynb` to the last section of
`chapter3.tex` and `sec:qubitmodels` of `chapter4.tex`,
and `models.ipynb` plus `Hubbard.ipynb` to
`chapter4.tex`, `fullci.ipynb` to `chapter5.tex`, `hartreefock.ipynb` to
`chapter6.tex`, `tdarpa.ipynb` to `chapter7.tex`, `dft.ipynb` to
`chapter8.tex`, `mbptheory.ipynb` to `chapter9.tex` and
`coupledcluster.ipynb` to `chapter10.tex` and `quantumdots.ipynb` to
`chapter11.tex`; the older set is `notation.ipynb`,
`secondquant.ipynb`, `fci.ipynb`, `hfock.ipynb`, `mbpt.ipynb`, `cc.ipynb`,
`vmcdmc.ipynb`, … New book material that includes runnable code should have a
notebook counterpart, and new notebooks must be registered in `_toc.yml`.

Note: `_toc.yml` lists `greensfunctions.ipynb`, which does not exist yet.

## BookMaterial/

Legacy **DocOnce** sources (`*.do.txt`) that predate the LaTeX chapters —
`secondquant.do.txt`, `fci.do.txt`, `hfock.do.txt`, `cc.do.txt`, `mbpt.do.txt`,
`notation.do.txt`, `quantumcomputing.do.txt`. The `.dlog` files are DocOnce
translation logs; `*.p.tex` and `*~` are generated/backup artefacts — never edit
those. `chapter2.tex` was derived from `secondquant.do.txt`, so that pair is the
reference example for converting DocOnce to book LaTeX.

`BookMaterial/Programs/` holds the programs referenced by the chapters. Those
written for chapters 1 and 2 are:

| File | Chapter | Contents |
|---|---|---|
| `direct_solvers.py` | 1 | Gaussian elimination, LU, Cholesky, tridiagonal |
| `iterative_solvers.py` | 1 | Jacobi, Gauss-Seidel, SOR, conjugate gradient |
| `householder.py` | 1 | Householder tridiagonalisation, `tqli`, power method |
| `lanczos.py` | 1 | Lanczos with reorthogonalisation, pairing model |
| `schrodinger_diagonalization.py` | 1 | Schrödinger by diagonalisation |
| `svd.py` | 1 | SVD, pseudoinverse, Schmidt decomposition, two-body factorisation |
| `slaterdeterminant.py` | 2 | Slater determinants, energy functional, minimal SCF |
| `slater_update.py` | 2 | Ratio $R$, Sherman-Morrison updates, nodal-surface stability |
| `wick.py` | 3 | Vacuum expectation values by anticommutation, by Wick contractions and by the generalised theorem; the generalised theorem also checked as an operator identity in a small Fock space |
| `jordanwigner.py` | 3, 4 | Jordan-Wigner operators as matrices, exact Pauli decomposition by traces, the anticommutator check, the hopping/pair/two-body identities, and the pairing and pairing+p-h Hamiltonians in both the 2L-qubit and the compact L-qubit encodings |
| `models.py` | 4 | Lipkin, pairing, pairing+particle-hole, Hubbard, Heisenberg and Calogero: matrices, spectra, exact limits, Jordan-Wigner check |
| `fci.py` | 5 | Determinants as bit strings, pairing and Hubbard FCI, truncations, size consistency, Hilbert-space growth |
| `hartreefock.py` | 6 | SCF with the density matrix, Thouless rotations, the stability matrix, BCH and Trotter errors, the electron gas |
| `rpa.py` | 7 | Pairing + particle-hole model: exact diagonalization, TDA, RPA, BCS, QRPA, all matrices built as double commutators |
| `dft.py` | 8 | LDA built from the uniform gas, the exchange hole, a Kohn-Sham solver in the chapter-6 basis, self-interaction and locality errors |
| `mbpt.py` | 9 | Rayleigh-Schroedinger to any order by recursion, Brillouin-Wigner self-consistently, four models, size extensivity, convergence |
| `coupledcluster.py` | 10 | General spin-orbital CCSD with intermediates, CCD, a validating FCI in the same basis, the order-by-order expansion, and `UnitaryCC` (UCCD/UCCSD, exact and Trotterised) |
| `quantumdot.py` | 11 | 2D quantum dot: closed-form Coulomb integrals, the antisymmetrised two-body matrix, a fast two-electron FCI, and the whole chapter-10 hierarchy applied to it. Takes a few minutes to run |

All of these run on `numpy` alone except `hartreefock.py`, which also needs
`scipy.linalg` for `expm`/`logm`; `rpa.py`, which needs `scipy.linalg`,
`scipy.sparse` and `scipy.optimize`; and `dft.py`, which needs `scipy.linalg`
and `scipy.interpolate`. `mbpt.py` imports `fci.py`, `hartreefock.py` and
`rpa.py` and reuses their model builders, so those three must stay importable;
`coupledcluster.py` optionally imports `fci.py`, `rpa.py` and `mbpt.py` for its
comparison table and degrades gracefully without them.
The notebooks additionally use `matplotlib`. `rpa.py` takes about 20 s to run
in full, everything else a few seconds.

The one-dimensional trap (harmonic well, softened Coulomb repulsion, eight
oscillator orbitals, 401-point grid) recurs in chapters 1, 2, 6 and 8, and
the same energies must come out of `slaterdeterminant.py`, `hartreefock.py`
and `dft.py`. If you change the grid or the softening, re-run all three.

Every number quoted in chapters 1–10 comes from one of these; each file runs
as a script and prints the tables that appear in the text. **If you change a
program, re-run it and update the corresponding table.** The older programs
`CCD_PairingModel.py`, `NeutronMatterCCD_Ladders.py`, `ucc.py`, `ho1dim.py` and
`ho2dim.py` are also here; `sd.py` sits one level up in `BookMaterial/`.

## House rules

- Match the surrounding prose voice: first person plural ("We introduce…",
  "It follows that…"), derivations written out step by step.
- Prefer editing existing chapter files over creating new ones; when a new
  chapter is genuinely needed, add it to `book.tex` under the right `\part`
  and to `BookMaterial/notes.md`.
- Leave `*~` backup files, `.dlog` logs and `.DS_Store` untouched.
