# Refactor Plan — `CausalCompetingRisks` v1

Approved 2026-04-27 after 4-agent adversarial review (`grill-claude-code` × API drift,
`simplify` × kill list, `improve-codebase-architecture` × pipeline, test-suite review)
and user grilling on file restructure / classical-vs-separable separation.

Companions: `dev/API_v1_LOCKED.md` (user-facing surface), `dev/ARCHITECTURE.md`
(internal file layout + dep graph).

10 stages. Each stage = small set of commits. Tests stay green between stages.

---

## Stage 0 — Pin contrast formulas (doc fix + regression test)

Initial adversarial review thought `R/contrasts.R` had wrong formulas.
Methodological re-derivation against Stensrud et al. (2020) §3 found the
**spec doc was wrong, not the code**. The code is methodologically faithful.

| # | Change | File |
|---|---|---|
| 0.1 | Correct "Pre-defined contrasts" section (SDE/SIE labels were swapped) | `dev/API_v1_LOCKED.md` ✓ done 2026-04-27 |
| 0.2 | Add regression test pinning the four formulas to fixture values | `tests/testthat/test-contrasts-formulas.R` |
| 0.3 | Update roxygen on `contrast()` to cite Stensrud 2020 §3 | `R/contrasts.R` |

Why first: pin the correct formulas in a test before any refactor moves them.
No code-formula change required.

---

## Stage 1 — File restructure (mechanical, no behavior change)

| # | Change |
|---|---|
| 1.1 | Create empty new files per `dev/ARCHITECTURE.md` |
| 1.2 | Move `fit_hazard_models()` + `predict_hazards()` → `hazards.R` |
| 1.3 | Extract `cif_from_hazards()` utility into `hazards.R` |
| 1.4 | Move `add_d_swap_weights` → `separable_swap.R` |
| 1.5 | Move classical IPW weight build → `weights.R` |
| 1.6 | Move `estimate_weighted_cum_inc_rep1/2` → `separable_ipw.R` |
| 1.7 | Move g-formula CIF math → `gformula_core.R` (classical) + `separable_gformula.R` (overlay) |
| 1.8 | Update `NAMESPACE` / `@export` tags |
| 1.9 | Run full test suite after every move |

Goal: final layout matches `dev/ARCHITECTURE.md`, behavior identical.

---

## Stage 2 — API rename + signature replace

| # | Change |
|---|---|
| 2.1 | Rename `causal_cr()` → `separable_effects()` |
| 2.2 | Replace signature with locked args |
| 2.3 | Drop dual-path branch + `pt_data`/`person_time` class |
| 2.4 | Drop `extreme_weight_*`; add `truncate = c(0.01, 0.99)` (default; `NULL` = no truncation) |
| 2.5 | Add `stabilize`, `K`; rename `censoring_weights` → `ipcw` |
| 2.6 | Drop warning-capture scaffolding (`withCallingHandlers`, `predict_with_warning`, `$warnings` slot); keep `warning()` calls and `$diagnostics` slot |

---

## Stage 3 — Fill IPW propensity gap (correctness)

Current code computes swap weights only. Locked spec requires `1/π(A|L)` (treatment weight).

| # | Change | File |
|---|---|---|
| 3.1 | Implement `fit_propensity()` + `predict_propensity()` | `R/propensity.R` |
| 3.2 | Wire `1/π(A|L)` into IPW pipeline | `R/weights.R` |
| 3.3 | Stabilize handling: `P(A) / π(A|L)` when `stabilize=TRUE` | `R/weights.R` |
| 3.4 | Coverage tests | `tests/testthat/test-propensity.R` |

---

## Stage 4 — Long-format `$risk` everywhere

| # | Change |
|---|---|
| 4.1 | `cif_from_hazards()` returns long format `(k, arm, cum_inc, method)` |
| 4.2 | Update `gformula_core.R`, `ipw_core.R`, separable estimators to emit long |
| 4.3 | Rewrite `risk()`, `contrast()` for long input (no-arg) |
| 4.4 | Update `print.R`, `plot.R` for long input |

---

## Stage 5 — Implementation defaults wired

| # | Change |
|---|---|
| 5.1 | `glm.control(maxit = 100)` in single funnel `fit_one_model()` |
| 5.2 | Truncation default `[0.01, 0.99]` applied in `weights.R` |
| 5.3 | IPW always runs Rep 1 + Rep 2 |
| 5.4 | Warning: Y/C formula RHS differ under gformula+ipcw=FALSE |
| 5.5 | Warning: method="ipw" + ipcw=FALSE (independent censoring) |

---

## Stage 6 — Bootstrap rewrite

| # | Change |
|---|---|
| 6.1 | Cluster-by-id resampling |
| 6.2 | Mutate fit (no separate class), `B`/`seed` args |
| 6.3 | Internal `.fit_separable_effects()` worker; bootstrap reuses (skips revalidation) |
| 6.4 | Store B × K × |arms| array on `$bootstrap` |

---

## Stage 7 — Accessors + S3 methods

| # | Change |
|---|---|
| 7.1 | `risk(fit)` no-arg, returns long df |
| 7.2 | `contrast(fit)` no-arg, returns long df, CIs from `$bootstrap` |
| 7.3 | New: `assumptions(fit)` accessor, prints assumption block |
| 7.4 | `print(fit)` minimal (single-line method) |
| 7.5 | `summary(fit)` pre/post-bootstrap contracts per locked doc |
| 7.6 | `plot(fit, arms, colors, labels, show_ribbons)` — returns ggplot, user-extensible via `+` |
| 7.7 | `risk_table(fit)` kept; slimmed |

---

## Stage 8 — Validation strict checklist

| # | Change |
|---|---|
| 8.1 | Implement strict locked validation (column existence, binary 0/1, monotone time, ≤1 event/subject, etc.) — `R/validate.R` |
| 8.2 | Implement conventional warnings (`min(time)≠0`, treatment varies, covariate NA on event row) |
| 8.3 | Drop `validate_subject_level()` entirely |
| 8.4 | Validate `formulas` keys ⊆ `{Y, D, C, A_model}` |

---

## Stage 9 — Test rebuild

Delete: `test-to-person-time.R`, `test-validate-subject-level.R`.

Build (target ~140 tests across 10 files):

| File | Tests | Tier |
|---|---|---|
| `test-validate.R` | ~30 | low |
| `test-covariate-quality.R` | ~8 | low |
| `test-fit-structure.R` | ~20 | medium |
| `test-method-ipcw-matrix.R` | ~8 | medium |
| `test-numerical-correctness.R` | ~10 | **HIGH** |
| `test-formulas.R` | ~10 | low |
| `test-accessors.R` | ~15 | medium |
| `test-bootstrap.R` | ~12 | medium/high |
| `test-edge-cases.R` | ~10 | low/medium |
| `test-reproducibility.R` | ~5 | low |

10 load-bearing tests live in `test-numerical-correctness.R` and `test-bootstrap.R`
(DGP recovery, Rep1≡Rep2 under isolation, CI coverage, cluster-by-id, seed reproducibility).

---

## Stage 10 — Documentation

| # | Change |
|---|---|
| 10.1 | Roxygen on every exported function (Hernán & Robins prose voice) |
| 10.2 | `README.md` with quick start |
| 10.3 | Vignette: basic use (Feynman + worked-example-first) |
| 10.4 | `NEWS.md` v1 changelog |

(Methodology walkthrough vignette already exists — not in scope here.)

---

## Deferred to post-v1

- **Benchmark vs `gfoRmula`** on shared DGP. Set up `dev/benchmark/`, compare numerical
  agreement, document any discrepancies.
- v2 features per `dev/API_v1_LOCKED.md` (Z_k generalization, time-varying covariates,
  composite, controlling direct effect, marginal net risk).

---

## Sequencing logic

- **Stage 0 first**: bugs are independent; ship correct numbers immediately.
- **Stages 1–8**: in order. Each unblocks the next (restructure → rename → propensity → long format → defaults → bootstrap → accessors → validation).
- **Stage 9 (tests)**: after structure + API stable, before docs. Tests written against rebuilt API, not legacy.
- **Stage 10**: post-implementation docs.

Each stage = a tracked GitHub issue (next step: `to-issues`).
