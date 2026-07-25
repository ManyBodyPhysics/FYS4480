---
name: manybody-code
description: Write or refactor Python, C++ and Fortran programs and Jupyter notebooks for the many-body book and the FYS4480/9480 course — CCD, pairing model, Hartree-Fock, FCI, VMC, unitary coupled cluster, harmonic-oscillator basis code. Use when adding a program to doc/BookManybody/BookMaterial/Programs or doc/Programs, writing a notebook chapter for the Jupyter-book in doc/LectureNotes, or embedding code listings into LaTeX.
---

# Code and notebooks for the many-body material

## Where code lives

| Path | Role |
|---|---|
| `doc/BookManybody/BookMaterial/Programs/` | Programs tied to the book: `CCD_PairingModel.py`, `NeutronMatterCCD_Ladders.py`, `ucc.py`, `ho1dim.py`, `ho2dim.py` |
| `doc/BookManybody/BookMaterial/sd.py` | Slater-determinant helper, sits directly in `BookMaterial/` |
| `doc/Programs/Python/`, `doc/Programs/Cpp/` | General course programs |
| `doc/Programs/PairingPH/`, `TDHF/`, `PySCF/`, `QDcoulomb/`, `PairingResolution/` | Topic-specific programs and notebooks |
| `doc/LectureNotes/*.ipynb` | Jupyter-book chapters |

## Language conventions

Python, C++ and Fortran, all in an **object-oriented** style. Prefer classes
that mirror the physics — a `SlaterDeterminant`, a `Hamiltonian`, a
`CoupledCluster` solver — over long procedural scripts, so the same structure
can be shown in all three languages.

- **Python**: NumPy/SciPy for linear algebra, matplotlib for figures. Keep the
  numerics readable; the code is teaching material first and fast code second.
  Comment the physics, not the syntax.
- **C++**: modern C++ with classes and RAII; Armadillo or Eigen for linear
  algebra where the existing code already uses them (check `doc/Programs/Cpp/`
  before choosing).
- **Fortran**: free-form Fortran 90+ with modules and derived types.

Where a program illustrates a method covered in the book, name the variables
after the symbols in the corresponding chapter (`t_ijab` for CCD amplitudes,
`h_pq`/`v_pqrs` for one- and two-body matrix elements, etc.).

## Notebooks and the Jupyter-book

The executable book lives in `doc/LectureNotes/` and is built with
`jupyter-book`:

- `_config.yml` — title, repo URL, `execute_notebooks: auto`, 30 s cell timeout,
  `allow_errors: false`. Notebooks must therefore **run clean and fast**.
- `_toc.yml` — `format: jb-book`, `root: intro`, chapters grouped in parts
  ("About the course", "Basics of Many-Body Physics", "Post Hartree-Fock
  Methods", …). Any new notebook must be added here or it will not appear.
- Front matter is markdown: `intro.md` (the root), `teachers.md`,
  `textbooks.md`. Chapters are notebooks: `notation.ipynb`,
  `secondquant.ipynb`, `fci.ipynb`, `hfock.ipynb`, `mbpt.ipynb`, `cc.ipynb`,
  `vmcdmc.ipynb`, `gradientmethods.ipynb`, `resamplingmethods.ipynb`,
  `vectorization.ipynb`, `parallelization.ipynb`.
- Known gap: `_toc.yml` lists `greensfunctions.ipynb`, which does not exist —
  the build will fail on it until the notebook is written or the entry removed.

Build:

```bash
cd doc/LectureNotes && jupyter-book build .
```

Notebook shape: markdown cell with the theory (LaTeX math is fine), then a code
cell that implements it, then a cell that runs a small case and prints or plots
the result. Keep example sizes small enough to finish within the timeout.

Never edit `_build/` or `.ipynb_checkpoints/`.

## Embedding code in the LaTeX book

`book.tex` defines `Python` and `C++` listing environments:

```latex
\begin{Python}{}
def hartree_fock(...):
    ...
\end{Python}
```

Prefer including the real program file over retyping it, and keep listings
short — point to the full program in `BookMaterial/Programs/` for the rest.

## Verifying

- Run the program or execute the notebook end to end before declaring done.
- For a notebook: `jupyter nbconvert --to notebook --execute <file>.ipynb`.
- Sanity-check the physics: compare against a known limit (e.g. CCD vs FCI on
  the pairing model, HF energy vs exact diagonalization for small systems).
  Existing programs in the repo often contain such reference numbers.
