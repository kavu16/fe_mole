# Architecture decision records

One record per choice that had real alternatives. Half a page each — if a
decision needs more than that, the extra belongs in `docs/theory/`.

Write one when a reasonable engineer would ask "why is it done this way?" and
the answer is not obvious from the code. Do not write one for a choice with no
alternative.

Records are immutable once merged. If a decision is reversed, add a new record
that supersedes it and mark the old one — the point of the log is that the
reasoning at the time stays legible, including the reasoning that turned out to
be wrong.

## Index

| # | Decision | Status |
|---|---|---|
| [0001](0001-particle-storage-layout.md) | Hand-rolled struct-of-arrays for particle storage | Accepted |
| [0002](0002-single-crate-layout.md) | Single crate rather than a Cargo workspace | Accepted |

## Template

```markdown
# NNNN — Title

**Status:** Proposed / Accepted / Superseded by NNNN
**Date:** YYYY-MM-DD
**Milestone:** MN

## Context

What forced a choice. Constraints that were actually binding.

## Options

Each real alternative, with the strongest argument *for* it — not a straw man.

## Decision

What was chosen, and the reason that actually decided it.

## Consequences

What this makes easy, what it makes hard, and what would trigger revisiting.
```
