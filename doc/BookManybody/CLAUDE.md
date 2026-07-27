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
                   Wick's theorem
chapter4.tex       Physical systems and models: Lipkin, pairing,
                   Fermi-Hubbard, Calogero-Sutherland
chapter5.tex       Full configuration interaction: the CI expansion, the
                   eigenvalue problem, the exponential wall, truncations
chapter6.tex       Mean-field approaches: Hartree-Fock, Thouless' theorem,
                   BCH and Trotter, stability, the infinite electron gas
chapter7.tex       Mean-field applications: the pairing + particle-hole model,
                   TDA, RPA, BCS and QRPA
add.tex  linalg.tex  eig.tex  vmcfermion2015.do.txt  version2024.do.txt
calogero.tex  SecondMidterm2025.tex  Hubbard.ipynb  combinatorial.tex
week40.do.txt  pmm.tex  pairing_ph_exact_diag.ipynb
pairing_ph_rpa_qrpa.ipynb  FYS4480*.pdf
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
   chapters 1–7 uses them. Dropping them would allow a plain `pdflatex` route.

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
- All seven chapters use numbered `\section{...}`. Chapters 1--3 end with an
  `Exercises` section containing warm-up exercises followed by the weekly
  FYS4480 sessions and their answers; chapters 4--7 relegate the long
  derivations to their exercises, each with a worked solution written as
  `\subsection*{Exercise N: title}` followed by `\paragraph{Answer to a).}`.
- `\Det{...}` expands to `{|\bfv{#1}|}`, so it takes a *matrix symbol*, not an
  arbitrary expression — `\Det{\bm{A}}` is fine, `\Det{\sum_k c_k b_{jk}}`
  breaks. Use `\det(...)` for expressions. `\bra`/`\ket` are not defined;
  use `\bbra`/`\bket`, `\braket{}{}` and `\bracket{}{}{}`.
- svmono provides `theorem`, used in chapter 6 for Thouless' theorem.
- `book.tex` also defines `\wmark`, `\wline`, `\wlineb` and `\wstrut` for the
  Wick contraction lines of Chapter 3, and a `notebox` environment for the
  many-body and quantum-computing connection boxes.
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
`pairing_ph_rpa_qrpa.ipynb` were merged into `chapter7.tex`.

## Jupyter-book target

`../LectureNotes/` holds the executable version: `_config.yml`, `_toc.yml`,
markdown front matter (`intro.md`, `teachers.md`, `textbooks.md`) and chapter
notebooks. `linearalgebra.ipynb` is the companion to `chapter1.tex` and
`manybodybasics.ipynb` to `chapter2.tex`, `wicktheorem.ipynb` to the Wick
sections of `chapter3.tex`, and `models.ipynb` plus `Hubbard.ipynb` to
`chapter4.tex`, `fullci.ipynb` to `chapter5.tex`, `hartreefock.ipynb` to
`chapter6.tex` and `tdarpa.ipynb` to `chapter7.tex`; the older set is
`notation.ipynb`,
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
| `wick.py` | 3 | Vacuum expectation values by anticommutation and by Wick contractions |
| `models.py` | 4 | Lipkin, pairing, Hubbard and Calogero: matrices, spectra, exact limits |
| `fci.py` | 5 | Determinants as bit strings, pairing and Hubbard FCI, truncations, size consistency, Hilbert-space growth |
| `hartreefock.py` | 6 | SCF with the density matrix, Thouless rotations, the stability matrix, BCH and Trotter errors, the electron gas |
| `rpa.py` | 7 | Pairing + particle-hole model: exact diagonalization, TDA, RPA, BCS, QRPA, all matrices built as double commutators |

All of these run on `numpy` alone except `hartreefock.py`, which also needs
`scipy.linalg` for `expm`/`logm` and `rpa.py`, which needs `scipy.linalg`,
`scipy.sparse` and `scipy.optimize`; the notebooks additionally use
`matplotlib`. `rpa.py` takes about 20 s to run in full.

Every number quoted in chapters 1–7 comes from one of these; each file runs
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
