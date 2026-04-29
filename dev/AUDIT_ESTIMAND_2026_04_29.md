# Estimand / identification audit — 2026-04-29

Adversarial methodological pass. Skills:
`separable-effects-ubiquitous-language`,
`causal-inference-ubiquitous-language`,
`survival-competing-risks-ubiquitous-language`, `grill-claude-code`.

8 sections. 3 deal-breakers, 3 followup-emails, 2 well-handled.

## 1. Estimand statement — DEAL-BREAKER

- The estimand is never named, never written in counterfactual notation,
  never present in `print()` / `summary()` / roxygen.
- User sees arm tuples (`arm_11`, `arm_10`, etc.) and learns the
  convention only by reading `dev/API_v1_LOCKED.md`.
- Stensrud's papers frame separable effects as a **new estimand**, not
  a procedure. Silence on the estimand reads as not having internalized
  this.
- Fix: implement the planned `assumptions()` accessor (TODO.md, deferred).
  Have it print the target counterfactual quantity in notation and the
  conditions under which the printed numbers identify it.

## 2. Identifying assumptions — DEAL-BREAKER

- No mention anywhere user-facing of: exchangeability, positivity (as
  such — see #6), consistency, generalized decomposition assumption,
  dismissible component conditions Δ1/Δ2, independent censoring.
- IPCW is described in roxygen as relaxing independent censoring, but
  this isn't surfaced in print/summary either.
- Fix: roxygen `@section Identification` on `separable_effects()`;
  `assumptions()` accessor block; `summary()` references it.

## 3. Isolation regime — FOLLOWUP EMAIL

- Architecture is ready: w_d and w_y are computed and stored.
- No accessor reports w_d / w_y proximity to 1, no diagnostic message
  classifies the regime (full / partial-D / partial-Y / no-iso+Z_k /
  not-identifiable per Stensrud & Young 2021 Lemma 4).
- TODO.md acknowledges and defers, but for v1 release with Stensrud's
  name attached, a regime read-out is the minimum.
- Fix: in `diagnostic(fit)` add summary statistics on w_d, w_y
  (mean, quantiles, max distance from 1) and a one-line regime hint.

## 4. Rep 1 vs Rep 2 — FOLLOWUP EMAIL

- Both are computed. Methodological status (alternative estimators of
  the same estimand under full isolation, divergent under partial /
  misspecification) is not stated to the user.
- No unit test pins down the agreement-under-full-isolation property
  on synthetic data.
- Fix: vignette section on Rep 1 vs Rep 2 interpretation; unit test on
  synthetic full-isolation DGP.

## 5. Decomposition A vs B — WELL-HANDLED (clarity gap)

- Both computed; `decomp` column labels A and B in the contrasts table.
- Print output shows the four arms but does not flag that A and B are
  different estimands (different reference for the unintervened
  component).
- Minor fix: print/summary footer notes the reference arm for each
  decomposition.

## 6. Positivity diagnostics — FOLLOWUP EMAIL

- `check_fitted_positivity()` (`R/hazards.R:143–200`) flags fitted P
  near 0/1 — but the warning isn't framed as a *positivity* diagnostic
  to the user.
- Rep 2 hazard-zero floor (`eps ≈ 1.5e-8` in `separable_swap.R:92,97,99,
  106,108`) silently absorbs positivity violations. TODO.md plans
  row-level flagging; not implemented.
- Fix: rename the warning category to `positivity`; add row-level
  flagging when Rep 2 hits the floor; report in `diagnostic()`.

## 7. Notation consistency with Stensrud — WELL-HANDLED

- Arm labels match `(a_Y, a_D)` tuples.
- Swap-weight construction (`R/separable_swap.R:42–57,89–115`) matches
  Stensrud 2020 Appendix.
- Lag trick: Y uses lagged cumprods, D uses contemporaneous — correct
  per the paper, documented in TODO.md.
- Minor drift: code uses `A_y`/`A_d` (lower) instead of paper's
  `A_Y`/`A_D` (subscript). Cosmetic; fix in roxygen pass.

## 8. What Stensrud would ask first

1. **What is the estimand?** — currently no good answer in package text.
2. **Are the dismissible component conditions named?** — no.
3. **Which isolation regime are we in?** — not surfaced.
4. **Is the lag trick correct in Rep 2?** — yes, verified.
5. **What happens at positivity boundary?** — partially handled, Rep 2
   floor swallows the signal.

## Bottom line

Math is correct. Documentation of the estimand and identifying
assumptions is the gap. Fixes are documentation + diagnostics, not
re-derivations.
