---
name: book-chapter
description: Write, extend or edit chapters of the many-body physics textbook in doc/BookManybody (book.tex, chapter1.tex, chapter2.tex, variationalmc.tex, Greensfunction.tex, preface.tex). Use when adding new sections or derivations, restructuring the book, adding a new chapter to a part, checking LaTeX macro/label consistency, or preparing material for the Springer svmono build.
---

# Writing textbook chapters

The book is *Quantum mechanics for Many-particle Systems: from standard methods
to quantum computing and machine learning*, a Springer monograph built from
`doc/BookManybody/book.tex` with the `svmono` class.

## Before writing

1. Read `doc/BookManybody/CLAUDE.md` — it lists the macro set and known build
   quirks.
2. Read `doc/BookManybody/BookMaterial/notes.md` — the planned outline for the
   whole book. New material should fit a bullet there; if it does not, say so
   and propose an outline update.
3. Open the neighbouring chapter (usually `chapter2.tex`) and match its voice,
   equation-label scheme and level of detail.
4. Check `doc/BookManybody/BookMaterial/*.do.txt` — a DocOnce draft of the topic
   may already exist. If so, use the `doconce-convert` skill.

## Chapter file anatomy

A chapter that is `\input` by `book.tex` has **no preamble**. It begins:

```latex
\chapter{Second quantization}

We introduce the time-independent operators ...
```

Everything else — packages, `\newcommand`s, `lstset` configuration — lives in
`book.tex`. Do not add `\documentclass`, `\begin{document}` or duplicate macro
definitions in a chapter file.

> `chapter1.tex` currently violates this (it is a standalone `article`). Treat
> `chapter2.tex` as the model.

## Macros available from book.tex

Use these rather than raw markup:

| Purpose | Macro |
|---|---|
| Equation shorthands | `\be … \ee`, `\bea … \eea`, starred `\beN/\eeN`, `\beaN/\eeaN`, `\bdm/\edm` |
| Operators | `\op{H}` (widehat), `\Op{H}` / `\OP{H}` (bold widehat) |
| Vectors / matrices | `\vec{r}`, `\bfv{r}`, `\matr{A}`, `\uvec{n}` |
| Expectation values | `\mean{H}`, sized: `\meanb \meanbb \meanbbb \meanbbbb` |
| Dirac notation | `\overlap{a}{b}`, `\bracket{a}{H}{b}`, `\projection{a}{b}`, `\brab{}`, `\ketb{}` plus `b/bb/bbb/bbbb` sizes |
| Misc | `\ud{x}`, `\udd{x}`, `\bigO`, `\real{z}`, `\prob`, `\cov`, `\var`, `\Det{A}`, `\PsiT` |
| Code | `\begin{Python}{}…\end{Python}`, `\begin{C++}{}…\end{C++}` |

If a needed macro is missing, add it to the macro block in `book.tex` — not to
the chapter.

## Equations and cross-references

- Numbered `equation` environments, labelled `\label{eq:<chap>-<n><letter>}`,
  e.g. `\label{eq:2-1a}`, `\label{eq:2-3b}`.
- Reference as `Eq.~(\ref{eq:2-1a})` — always with the tie and parentheses.
- Keep derivations explicit: state the assumption, show the algebra, state the
  result. The existing text walks through anticommutation relations line by
  line; match that granularity.

## Style

- First person plural: "We introduce…", "We define…", "It follows that…".
- Notation: `a^{\dagger}_\alpha`, `|0\rangle`,
  `|\alpha_1\dots\alpha_n\rangle_{\mathrm{AS}}`, `\{A,B\}` for anticommutators,
  `[A,B]` for commutators.
- Add `\index{term}` for each new concept (the book runs `\makeindex`).
- Cite with `\cite{}`; the bibliography style is `unsrt`, with Springer's
  `spphys`/`spbasic`/`spmpsci` available in the folder.

## Adding a new chapter

1. Create `doc/BookManybody/<name>.tex` starting with `\chapter{...}`.
2. Add `\input{<name>}` to `book.tex` under the correct `\part`:
   *Mathematical basis and essential computational elements* → *Wave function
   based methods* → *Stochastic methods* → *Machine Learning based approaches*
   → *Quantum Computing*.
3. Add a matching bullet to `BookMaterial/notes.md`.
4. If the chapter contains runnable code, plan a Jupyter-book counterpart in
   `doc/LectureNotes/` and register it in `_toc.yml` (see the `manybody-code`
   skill).

## Verifying

Try a compile before declaring done:

```bash
cd doc/BookManybody && pdflatex -interaction=nonstopmode book.tex
```

Expect pre-existing failures from `chapter1.tex` and the missing
`foreword.tex`. Compare against the state before the edit rather than demanding
a clean run, and report anything *new*.

If `pdflatex` is unavailable, at minimum check: balanced environments, every
`\ref` has a `\label`, no duplicate labels, no locally redefined macros.
