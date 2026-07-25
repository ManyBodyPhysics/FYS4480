---
name: course-exercises
description: Create or edit weekly exercise sheets, answer sets and midterms for FYS4480/FYS9480 in doc/BookManybody/Exercises (ExercisesWeekNN.tex, AnswersWeekNN.tex, FirstMidterm2025.tex). Use when writing new problems, drafting solutions, adapting a week's set to a new semester, or working with the Lipkin, pairing or electron-gas model exercises.
---

# FYS4480/9480 exercises and answers

Files live in `doc/BookManybody/Exercises/`. Every file is **standalone** — it
carries its own preamble and compiles on its own, unlike the book chapters.

## Two document classes

**Exercise sheets** — `ExercisesWeekNN.tex`, `FirstMidterm2025.tex`,
`SecondMidterm2025.tex`:

```latex
\documentclass[prc]{revtex4}
\usepackage[dvips]{graphicx}
\usepackage{mathrsfs}
\usepackage{amsfonts}
\usepackage{lscape}
\usepackage{epic,eepic}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage[dvips]{epsfig}
\usepackage[T1]{fontenc}
\usepackage{hyperref}
\usepackage{bezier}
\usepackage{pstricks}
\usepackage{dcolumn}
\usepackage{bm}

\newcommand{\One}{\hat{\mathbf{1}}}
\newcommand{\eff}{\text{eff}}
\newcommand{\Heff}{\hat{H}_\text{eff}}
\newcommand{\Veff}{\hat{V}_\text{eff}}
\newcommand{\braket}[1]{\langle#1\rangle}
\newcommand{\Span}{\operatorname{sp}}
\newcommand{\tr}{\operatorname{trace}}
\newcommand{\diag}{\operatorname{diag}}
\newcommand{\bra}[1]{\left\langle #1 \right|}
\newcommand{\ket}[1]{\left| #1 \right\rangle}
\newcommand{\element}[3]{\bra{#1}#2\ket{#3}}
\newcommand{\normord}[1]{\left\{#1\right\}}

\begin{document}
\title{Exercises FYS4480/9480, week 38, September 15-19, 2025}
\maketitle
```

**Answer sheets** — `AnswersWeekNN.tex`, `AnswerLipkin.tex`:
`\documentclass[a4paper, 11pt, notitlepage, english]{article}` with a heavier
preamble (`babel`, `geometry`, `framed`, `multicol`, and `listings` configured
for Python: `\lstset{language=python}`, `frame=single`, red bold keywords,
italic blue comments).

**Do not invent a preamble.** Copy it verbatim from the nearest existing week
and adjust only the title.

## Structure

- Title line: `Exercises FYS4480/9480, week NN, <Month D-D, YYYY>` — week number
  and date range must match the UiO autumn-semester schedule. Existing sets run
  weeks 34–45, with midterms as separate files.
- Problems as `\subsection*{Exercise 1}`, `\subsection*{Exercise 2}`, …
- Weeks build on each other; a sheet often opens with "This exercise is a
  continuation of the exercises from last week on …". Read the previous week's
  file before writing a new one so notation and the model setup carry over.
- Answers mirror the exercise numbering one-to-one.

## Recurring model systems

- **Lipkin model** — spin operators `J_\pm`, `J_z`, states `\ket{J,J_z}`, four
  fermions in two levels. `AnswerLipkin.tex` is the reference treatment.
- **Pairing model** — `pairingfig.tex` holds the level diagrams.
- **Electron gas** — `electrongas.tex`.
- Simple harmonic-oscillator systems (`BookMaterial/ho1dim.py`, `ho2dim.py`).

Reuse the established notation, e.g.
`\ket{\Phi_{J_z=-2}} = a_{1-}^{\dagger}a_{2-}^{\dagger}a_{3-}^{\dagger}a_{4-}^{\dagger}\ket{0}`,
with the raising/lowering relations spelled out explicitly for the students.

## Writing new problems

- State the physical setup, then ask for a derivation, then (often) ask for a
  numerical implementation in Python.
- Give students the operator relations they need rather than assuming recall —
  the existing sheets restate the `J_\pm\ket{J,J_z}` formulas inline.
- Where a problem asks for code, point to the relevant program in
  `doc/BookManybody/BookMaterial/Programs/` or `doc/Programs/Python/`, and use
  the `manybody-code` skill for the implementation.
- Tie each exercise to the chapter that covers it (`chapter2.tex` for second
  quantization, etc.).

## Verifying

```bash
cd doc/BookManybody/Exercises && pdflatex -interaction=nonstopmode ExercisesWeekNN.tex
```

Check that every macro used is defined in that file's own preamble — these
files get no help from `book.tex`.
