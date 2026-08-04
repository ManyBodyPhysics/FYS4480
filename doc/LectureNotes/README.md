# Lecture notes: the executable companion to the book

This directory is a [Jupyter Book](https://jupyterbook.org).  Each notebook is
the runnable counterpart of one chapter of *Quantum Mechanics for
Many-Particle Systems*, and each imports the corresponding program from
`../BookManybody/BookMaterial/Programs`, so the numbers on these pages and the
numbers in the book come from the same code.

## Building

```bash
pip install -r requirements.txt
jupyter-book build .
```

Open `_build/html/index.html`.  A full build from cold takes about four
minutes, almost all of it in `quantumdots.ipynb`.

The configuration uses `execute_notebooks: cache`, so the outputs are stored
in `_build/.jupyter_cache` and later builds only re-run the notebooks that
have changed.  To force everything to run again:

```bash
jupyter-book build --all .
```

To start completely clean:

```bash
rm -rf _build && jupyter-book build .
```

## Contents

| Notebook | Chapter | Companion program |
|---|---|---|
| `linearalgebra.ipynb` | 1, linear algebra and eigenvalue problems | `direct_solvers.py`, `iterative_solvers.py`, `householder.py`, `lanczos.py`, `svd.py` |
| `manybodybasics.ipynb` | 2, Slater determinants and the energy functional | `slaterdeterminant.py`, `slater_update.py` |
| `wicktheorem.ipynb` | 3, Wick's theorem and its generalisation | `wick.py` |
| `jordanwigner.ipynb` | 3, from fermions to qubits | `jordanwigner.py` |
| `models.ipynb` | 4, Lipkin, pairing, Hubbard, Heisenberg, Calogero | `models.py` |
| `fullci.ipynb` | 5, full configuration interaction | `fci.py` |
| `hartreefock.ipynb` | 6, Hartree-Fock and Thouless' theorem | `hartreefock.py` |
| `tdarpa.ipynb` | 7, TDA, RPA, BCS and QRPA | `rpa.py` |
| `dft.ipynb` | 8, density functional theory | `dft.py` |
| `mbptheory.ipynb` | 9, many-body perturbation theory | `mbpt.py` |
| `coupledcluster.ipynb` | 10, coupled cluster and unitary CC | `coupledcluster.py` |
| `quantumdots.ipynb` | 11, the two-dimensional quantum dot | `quantumdot.py` |
| `statistics.ipynb` | 12, statistics, random walks and Metropolis | `montecarlo.py` |

## Notes for anyone editing these

* **The notebooks are executed from this directory.**  They add
  `../BookManybody/BookMaterial/Programs` to `sys.path`, which only resolves
  if the working directory is `doc/LectureNotes`.  `run_in_temp` is therefore
  set to `false` in `_config.yml`; do not change it.

* **`allow_errors` is `false`.**  A cell that raises stops the build.  That is
  deliberate: these notebooks are the regression test for the programs, and a
  silent failure would leave a wrong number on a page.

* **Watch out for notebooks written without trailing newlines.**  A cell whose
  `source` list has lines that do not end in `\n` renders as one glued-together
  line and fails with a `SyntaxError`.  To check the whole directory:

  ```python
  import json, glob
  for f in sorted(glob.glob("*.ipynb")):
      nb = json.load(open(f))
      bad = sum(1 for c in nb["cells"] for l in c["source"][:-1]
                if not l.endswith("\n"))
      if bad:
          print(f, bad, "glued lines")
  ```

* **`only_build_toc_files` is `true`**, so notebooks that are not listed in
  `_toc.yml` are ignored.  `Hubbard.ipynb`, `notation.ipynb` and
  `gradientmethods.ipynb` predate the book and are currently left out; add
  them to `_toc.yml` if you want them on the site.

* **The two expensive notebooks** are `quantumdots.ipynb` (about ninety
  seconds, dominated by the six-electron coupled-cluster runs at forty-two
  orbitals) and `statistics.ipynb` (about thirty-five seconds, dominated by
  the variational Monte Carlo).  If a build times out, raise
  `execute.timeout` in `_config.yml` rather than shrinking the calculations.
