# CLAUDE.md

> **Note:** Core project architecture and coding standards are defined in [AGENTS.md](./AGENTS.md).

@AGENTS.md

This file holds what is specific to Claude Code: how to work in this
repository. Commands and everything about the code itself are in
AGENTS.md.

## Workflow

- Run single test files while iterating; run `pixi run test` before
  calling a change done.
- Port a geometry-central file with its tests: the GoogleTest cases
  that cover it, and reference values from geometry-central itself
  where the test has none (see `docs/plan.md`, Testing).
- After changing doc strings, run `pixi run doc-check`; if it improved,
  `pixi run doc-check --update` lowers the baseline.
- A change molina needs (a solver option, a bug) goes into molina,
  with its tests there, not into a copy here; then move its `rev` in
  `pixi.toml`.
- Don't import mojodec: MojoDEC may depend on mogeometry, not the
  other way.

## Where to look

| Question                                   | Where                                   |
| ------------------------------------------ | --------------------------------------- |
| Commands, architecture, conventions, traps | AGENTS.md                               |
| What to port next, open design questions   | `docs/plan.md`                          |
| Doc string style                           | `docs/writing-docs.md`                  |
| The original C++                           | geometry-central, commit in `NOTICE`    |
| Sparse matrices and solvers                | `../molina/AGENTS.md`, Architecture     |
