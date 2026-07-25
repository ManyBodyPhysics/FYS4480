# ManyBodyPhysics / FYS4480

Course and textbook repository for **FYS4480 (MSc) / FYS9480 (PhD)**, *Quantum
mechanics for many-particle systems*, University of Oslo. Author: Morten
Hjorth-Jensen (mhjensen@uio.no). Remote: https://github.com/ManyBodyPhysics/FYS4480

## Layout

| Path | Contents |
|---|---|
| `doc/BookManybody/` | **The textbook.** Springer `svmono` LaTeX source, exercises, book material. See `doc/BookManybody/CLAUDE.md`. |
| `doc/LectureNotes/` | Jupyter-book: `_config.yml`, `_toc.yml`, `*.ipynb` chapters, `_build/` |
| `doc/Programs/` | Standalone programs: `Python/`, `Cpp/`, `PairingPH/`, `TDHF/`, `PySCF/`, `QDcoulomb/`, `PairingResolution/` |
| `doc/Exercises/` | Course exercises by year (`2022/`–`2025/`, `Templates/`) — distinct from the book's `doc/BookManybody/Exercises/` |
| `doc/pub/` | Published weekly material (`week36/`–`week39/`, `fci/`, `intro/`) |
| `doc/slides/`, `doc/Notes/`, `doc/HandwrittenNotes/`, `doc/Literature/`, `doc/Images/`, `doc/Projects/` | Supporting material |

## Working conventions

- **Programming languages**: Python, C++ and Fortran, written in an
  object-oriented style. Prefer Jupyter notebooks and Jupyter-book for
  presenting and discussing code.
- **Most work happens in `doc/BookManybody/`.** Read
  `doc/BookManybody/CLAUDE.md` before editing anything there — it documents the
  LaTeX conventions, the macro set, and known build quirks.
- Do not commit or push unless asked.
- `_build/`, `*~` backup files and `.DS_Store` are noise; leave them alone.

## Skills

Project skills live in `.claude/skills/`:

- `book-chapter` — writing/editing textbook chapters in `doc/BookManybody/`
- `course-exercises` — weekly exercise and answer sets for FYS4480/9480
- `manybody-code` — Python/C++/Fortran programs and notebooks for the book
- `doconce-convert` — handling the legacy DocOnce (`*.do.txt`) material
