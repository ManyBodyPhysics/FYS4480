---
name: doconce-convert
description: Work with the legacy DocOnce sources in doc/BookManybody/BookMaterial (*.do.txt — secondquant, fci, hfock, cc, mbpt, notation, quantumcomputing). Use when converting DocOnce material into an svmono book chapter or a Jupyter notebook, salvaging text from the old drafts, or interpreting the .dlog logs and *.p.tex generated files.
---

# DocOnce material in BookMaterial/

`doc/BookManybody/BookMaterial/` holds the original **DocOnce** sources that
predate the LaTeX chapters. They are the raw material for both the book
chapters and the `doc/LectureNotes/` notebooks.

## What is what

| Pattern | Meaning |
|---|---|
| `*.do.txt` | DocOnce source — the thing to read and convert |
| `*.do.txt~` | Editor backup — ignore, never edit |
| `*.dlog` | DocOnce translation log (records the target format and exercise count) |
| `*.p.tex` | Generated LaTeX from DocOnce — never edit by hand |
| `notes.md` | The planned outline of the whole book |

Sources present: `secondquant.do.txt`, `fci.do.txt`, `hfock.do.txt`,
`cc.do.txt`, `mbpt.do.txt`, `notation.do.txt`, `quantumcomputing.do.txt`.

## DocOnce syntax you will meet

| DocOnce | LaTeX equivalent |
|---|---|
| `======= Title =======` | `\chapter{Title}` (7 `=`) |
| `===== Title =====` | `\section{Title}` (5 `=`) |
| `=== Title ===` | `\subsection{Title}` (3 `=`) |
| `!bt` … `!et` | delimiters around raw LaTeX math — drop them, keep the inside |
| `!bc` … `!ec` | code block — becomes a `Python`/`C++` listing environment |
| `label{eq:2-1a}` | `\label{eq:2-1a}` (add the backslash) |
| `ref{eq:2-1a}` | `\ref{eq:2-1a}` (add the backslash) |
| `!bsubex` / `!esubex`, `!bans` / `!eans` | exercise sub-parts and answers |
| `_bold_`, `*emphasis*` | `\textbf{}`, `\emph{}` |
| `URL:"text":"link"` | `\href{link}{text}` |

Inline math (`$...$`) and full `equation` environments pass through unchanged.

## Converting to a book chapter

`chapter2.tex` was produced from `secondquant.do.txt` — **diff those two to see
exactly what the conversion looks like** before converting anything else.

Procedure:

1. Replace the top-level `=======` heading with `\chapter{...}`, demote the
   rest of the heading levels accordingly.
2. Strip every `!bt`/`!et` pair, keeping the LaTeX between them.
3. Add the missing backslashes to `label{}` and `ref{}`; rewrite references as
   `Eq.~(\ref{eq:...})`.
4. Convert `!bc`/`!ec` blocks to the `Python` or `C++` listing environments
   defined in `book.tex`.
5. Replace ad-hoc notation with the book's macros (see the `book-chapter`
   skill) — but do not add a preamble; chapter files have none.
6. Move exercises out of the chapter into `doc/BookManybody/Exercises/` using
   the `course-exercises` skill, unless they are meant as in-text problems.
7. Add `\input{<name>}` to `book.tex` under the right `\part`.

## Converting to a notebook

DocOnce originally generated the `doc/LectureNotes/*.ipynb` chapters (the
`.dlog` files record `Translating doconce text in X.do.txt to ipynb`). If a
notebook for the topic already exists, **edit the notebook directly** rather
than regenerating — the notebooks have diverged from the DocOnce sources.

If `doconce` is installed, the historical command was of the form
`doconce format ipynb <file>.do.txt`; verify the tool is available before
relying on it, and treat hand-conversion as the default.

## Rules

- Never edit `*.p.tex`, `*.dlog` or `*~` files.
- The DocOnce sources are the older draft; where they disagree with an existing
  chapter or notebook, the chapter/notebook wins. Flag substantive conflicts
  instead of silently picking one.
