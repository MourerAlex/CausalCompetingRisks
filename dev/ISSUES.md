# GitHub Issues — `CausalCompetingRisks` v1 refactor

Copy-paste each block below into a new GitHub issue at
https://github.com/MourerAlex/CausalCompetingRisks/issues/new

11 issues total: Stage 0 + Stages 1–10 (per `dev/REFACTOR_PLAN.md`).

Suggested labels: `refactor`, `bugfix`, `tests`, `docs`, `correctness`.
Suggested milestone: `v0.1.0`.

---

## Issue 0 — Pin Decomposition A & B contrast formulas (doc fix + regression test)

**Labels:** `docs`, `tests`, `correctness`
**Depends on:** none
**Blocking:** all downstream stages (only because we want the regression test in place before refactoring)

### Context

Adversarial review (2026-04-27) initially flagged `R/contrasts.R` as having
wrong contrast formulas relative to `dev/API_v1_LOCKED.md`. Methodological
re-derivation against Stensrud et al. (2020) §3 found the **opposite**: the
**code is correct**; the **spec doc was wrong**.

The doc had SDE/SIE labels swapped at every position (it labelled the
indirect-effect arms as "direct" and vice versa). Fixed in the spec on the
same date.

Correct formulas (now in `dev/API_v1_LOCKED.md`):

|          | Decomp A (vary A_Y first; intermediate arm = arm_10) | Decomp B (vary A_D first; intermediate arm = arm_01) |
|----------|---------------------------|---------------------------|
| SDE      | `arm_10 − arm_00`         | `arm_11 − arm_01`         |
| SIE      | `arm_11 − arm_10`         | `arm_01 − arm_00`         |
| Total    | `arm_11 − arm_00`         | `arm_11 − arm_00`         |

Convention: arm `"xy"` = `(a_Y = x, a_D = y)`. `A_Y` is the non-D-mediated path
(its effect = SDE). `A_D` is the D-mediated path (its effect = SIE).

### Definition of done

- [x] `dev/API_v1_LOCKED.md` "Pre-defined contrasts" section corrected (2026-04-27).
- [ ] New test file `tests/testthat/test-contrasts-formulas.R` pins the four formulas with a fixture where each arm CIF is known. Example: `arm_11 = 0.4`, `arm_10 = 0.2`, `arm_01 = 0.3`, `arm_00 = 0.1` → SDE-A = 0.1, SIE-A = 0.2, SDE-B = 0.1, SIE-B = 0.2; Total = 0.3.
- [ ] Algebraic identity tested: `SDE + SIE = Total` under each decomposition.
- [ ] Roxygen on `contrast()` accessor cites Stensrud et al. (2020) §3.
- [ ] No code change to `R/contrasts.R` required (formulas already correct).

### Files

- `tests/testthat/test-contrasts-formulas.R` (new)
- `R/contrasts.R` — roxygen update only (citation; no formula change)
- `dev/API_v1_LOCKED.md` — already corrected

---

## Issue 1 — File restructure to architecture spec

**Labels:** `refactor`
**Depends on:** Issue 0
**Blocking:** Issues 2–10

### Context

Reorganize `R/` per `dev/ARCHITECTURE.md` so classical estimators are separated
from the separable-effects overlay. Mechanical move only — no behavior change.

### Final layout (15 files, all in `R/`, flat)

Core:
`hazards.R`, `propensity.R`, `weights.R`, `gformula_core.R`, `ipw_core.R`

Separable overlay:
`separable_arms.R`, `separable_swap.R`, `separable_gformula.R`, `separable_ipw.R`

Inference:
`bootstrap.R`, `contrasts.R`

Validation:
`validate.R`

API surface:
`separable_effects.R`, `accessors.R`, `risk_table.R`, `print.R`, `plot.R`

Removed: `data_prep.R`.

### Definition of done

- [ ] All 15 files exist; legacy `causal_cr.R`, `data_prep.R`, `hazard_models.R`, `ipw.R` removed.
- [ ] Dependency graph matches `dev/ARCHITECTURE.md` (verified by inspection — no cycles).
- [ ] `cif_from_hazards()` extracted as shared utility in `hazards.R`.
- [ ] `NAMESPACE` regenerated; `@export` tags on the same symbols as before.
- [ ] Full existing test suite passes after move (no behavior change).

---

## Issue 2 — API rename + signature replacement

**Labels:** `refactor`, `breaking-change`
**Depends on:** Issue 1

### Context

Replace `causal_cr()` with `separable_effects()` per `dev/API_v1_LOCKED.md`.

### Definition of done

- [ ] `causal_cr()` renamed to `separable_effects()`; old name removed (no alias — no users yet).
- [ ] Signature exactly matches `dev/API_v1_LOCKED.md`:
  ```r
  separable_effects(data, id, time_name, treatment, event_y, event_d,
                    covariates, method, ipcw, stabilize, formulas, K,
                    truncate)
  ```
- [ ] `truncate = c(0.01, 0.99)` default; `NULL` opts out (Cole & Hernán 2008).
- [ ] `stabilize = TRUE` default.
- [ ] `K = NULL` defaults to `max(data[[time_name]])`.
- [ ] `censoring_weights` arg renamed to `ipcw`.
- [ ] `extreme_weight_adjust` and `extreme_weight_threshold` removed.
- [ ] Dual-path branch (`pt_data`/`causal_cr_pt`) removed; long-format data.frame only.
- [ ] Warning-capture scaffolding removed (`withCallingHandlers`, `predict_with_warning`, `$warnings` slot). `warning()` calls and `$diagnostics` slot retained.

### Files

- `R/separable_effects.R` (new entry point, formerly `causal_cr.R`)

---

## Issue 3 — Implement IPW propensity model `π(A|L)`

**Labels:** `correctness`, `feature`
**Depends on:** Issue 2

### Context

Locked API requires `1/π(A|L)` treatment weights for IPW. Current code computes
swap weights only — IPW estimator is not identified per locked spec.

### Definition of done

- [ ] `R/propensity.R` implements `fit_propensity()` and `predict_propensity()`.
- [ ] Default formula `A_model ~ L1 + L2` (auto-built from `covariates`); user override via `formulas$A_model`.
- [ ] `R/weights.R` combines treatment weights with IPCW (and swap weights when separable).
- [ ] `stabilize = TRUE` applies `P(A) / π(A|L)`; `stabilize = FALSE` applies `1/π(A|L)`.
- [ ] `glm.control(maxit = 100)` applied (see Issue 5).
- [ ] `tests/testthat/test-propensity.R`: covers default formula, override, stabilize on/off, propensity diagnostics in `$diagnostics$min_propensity`.

### Files

- `R/propensity.R` (new content)
- `R/weights.R`
- `tests/testthat/test-propensity.R` (new)

---

## Issue 4 — Long-format `$risk` everywhere

**Labels:** `refactor`, `breaking-change`
**Depends on:** Issue 3

### Context

Current return shape is wide per-method data.frames keyed by arm columns
(`arm_11`, `arm_00`, etc.). Locked API requires long format `(k, arm, cum_inc, method)`
with arm strings `"11"`, `"00"`, `"10"`, `"01"`.

### Definition of done

- [ ] `cif_from_hazards()` returns long format directly.
- [ ] All estimators (gformula, IPW Rep1, IPW Rep2) emit long-format CIF.
- [ ] `fit$risk` is a single long data.frame (not a list of per-method wide dfs).
- [ ] `risk(fit)` returns the long df unchanged (no-arg signature).
- [ ] `contrast(fit)` consumes long input; returns long contrasts df.
- [ ] `print.causal_cr_fit` and `plot.causal_cr_fit` consume long input.
- [ ] All structural tests updated for long shape.

### Files

- `R/hazards.R`
- `R/gformula_core.R`, `R/ipw_core.R`, `R/separable_gformula.R`, `R/separable_ipw.R`
- `R/accessors.R`, `R/print.R`, `R/plot.R`
- All tests touching `$cumulative_incidence` or `arm_11`-style columns

---

## Issue 5 — Implementation defaults + warnings

**Labels:** `correctness`, `polish`
**Depends on:** Issue 4

### Definition of done

- [ ] `glm.control(maxit = 100)` wired into single funnel `fit_one_model()` (used by `hazards.R` and `propensity.R`).
- [ ] Truncation `[0.01, 0.99]` default applied in `weights.R`; `truncate = NULL` opts out.
- [ ] IPW always runs both Rep 1 and Rep 2 (no user choice; sensitivity comparison is built-in).
- [ ] Warning emitted when `Y` and `C` formula RHS differ under `method = "gformula"` and `ipcw = FALSE`. Message references Bang & Robins 2005 / Kawahara et al. 2020.
- [ ] Warning emitted when `method = "ipw"` and `ipcw = FALSE` (independent censoring assumption).
- [ ] `tests/testthat/test-formulas.R`: both warnings fire on the documented inputs.

### Files

- `R/hazards.R`, `R/propensity.R`, `R/weights.R`, `R/separable_effects.R`
- `tests/testthat/test-formulas.R` (new)

---

## Issue 6 — Bootstrap rewrite

**Labels:** `refactor`, `feature`
**Depends on:** Issue 5

### Definition of done

- [ ] `bootstrap(fit, B = 500, seed = NULL)` cluster-resamples by `id`.
- [ ] Returns the same `causal_cr_fit` object with `$bootstrap` slot populated (no separate `causal_cr_bootstrap` class).
- [ ] `$bootstrap$replicates` is a B × K × |arms| array; `$bootstrap$ci_summary` provides percentile CIs.
- [ ] Internal `.fit_separable_effects()` worker introduced; `bootstrap()` calls it directly per replicate (skips revalidation).
- [ ] `tests/testthat/test-bootstrap.R`: cluster-by-id resampling, seed reproducibility, B × K × |arms| shape, CIs propagate to `summary()`.

### Files

- `R/bootstrap.R`
- `R/separable_effects.R` (worker exposed internally)
- `tests/testthat/test-bootstrap.R`

---

## Issue 7 — Accessors + S3 methods

**Labels:** `refactor`, `feature`
**Depends on:** Issue 6

### Definition of done

- [ ] `risk(fit)` no-arg, returns `fit$risk` (long df).
- [ ] `contrast(fit)` no-arg, returns long df of all 5 predefined contrasts (Total, SDE-A, SIE-A, SDE-B, SIE-B); `lower`/`upper` populated iff `fit$bootstrap` present, else `NA`.
- [ ] `assumptions(fit)` new accessor; prints identifying assumptions block (per `dev/API_v1_LOCKED.md`).
- [ ] `print(fit)` minimal: `<causal_cr_fit>`, single-line `Method: gformula` or `Method: ipw`, plus accessor hint.
- [ ] `summary(fit)` full breakdown; pre-bootstrap shows arms only with note "run bootstrap()"; post-bootstrap shows contrasts with 95% CIs.
- [ ] `plot(fit, arms = NULL, colors = NULL, labels = NULL, show_ribbons = TRUE)` returns ggplot; ribbons iff `fit$bootstrap` present; user-extensible via `+`.
- [ ] `risk_table(fit)` slimmed; returns long df of (k, arm, n_at_risk, n_events_y, n_events_d, n_censored).

### Files

- `R/accessors.R`, `R/print.R`, `R/plot.R`, `R/risk_table.R`
- `tests/testthat/test-accessors.R`

---

## Issue 8 — Validation strict checklist

**Labels:** `feature`
**Depends on:** Issue 7

### Definition of done

`R/validate.R` enforces the strict locked checklist (`dev/API_v1_LOCKED.md` §
Validation):

**Errors:**
- [ ] `data` is data.frame-like.
- [ ] All named columns exist.
- [ ] `treatment`, `event_y`, `event_d` ∈ {0, 1} only.
- [ ] Monotone increasing `time_name` within each subject.
- [ ] At most one event per subject (no rows after Y=1 or D=1).
- [ ] `method` ∈ `{"gformula", "ipw"}`.
- [ ] `formulas` keys ⊆ `{Y, D, C, A_model}`.

**Warnings:**
- [ ] `min(time_name) ≠ 0`.
- [ ] `treatment` varies within subject.
- [ ] Covariates missing on event rows.

- [ ] `validate_subject_level()` and `validate_input_shape` (where redundant) removed.
- [ ] `tests/testthat/test-validate.R` covers each rule with a pass and a fail case.

### Files

- `R/validate.R`
- `tests/testthat/test-validate.R`

---

## Issue 9 — Test suite rebuild

**Labels:** `tests`
**Depends on:** Issues 0–8 (some sub-tasks already drafted as part of those issues; this issue closes any remaining gaps)

### Context

Existing suite is ~80% out-of-scope under the locked API. Target ~140 tests
across 10 files. 10 load-bearing tests anchor the math and inference.

### Definition of done

**Delete** (out of scope):
- [ ] `tests/testthat/test-to-person-time.R`
- [ ] `tests/testthat/test-validate-subject-level.R`

**Build / verify** (some already created in earlier issues):

| File | Tests | Tier | Status |
|---|---|---|---|
| `test-validate.R` | ~30 | low | Issue 8 |
| `test-covariate-quality.R` | ~8 | low | new |
| `test-fit-structure.R` | ~20 | medium | new |
| `test-method-ipcw-matrix.R` | ~8 | medium | new |
| `test-numerical-correctness.R` | ~10 | **HIGH** | new |
| `test-formulas.R` | ~10 | low | Issue 5 |
| `test-accessors.R` | ~15 | medium | Issue 7 |
| `test-bootstrap.R` | ~12 | medium/high | Issue 6 |
| `test-edge-cases.R` | ~10 | low/medium | new |
| `test-reproducibility.R` | ~5 | low | new |

**Load-bearing 10** (must pass):
1. gformula recovers true CIF (4 arms) on hand-built DGP within MC tolerance.
2. IPW Rep1 recovers true CIF.
3. SDE-A point estimate within tolerance + sum identity (SDE-A + SIE-A = Total).
4. IPW Rep1 ≡ Rep2 under full-isolation DGP.
5. Bootstrap CIs cover true SDE-A at nominal rate (small-MC).
6. Bootstrap reproducibility under fixed seed.
7. Cluster-by-id resampling unit (not row-level).
8. `$risk` long-format shape (`k, arm, cum_inc, method`).
9. `contrast(fit)` emits the 5 named decomposition contrasts with correct lower/upper.
10. `method="ipw"+ipcw=FALSE` warning AND Y/C formula mismatch warning fire.

- [ ] All 10 load-bearing tests pass.
- [ ] Total test count ≥ 140.

### Files

- `tests/testthat/*.R`

---

## Issue 10 — Documentation

**Labels:** `docs`
**Depends on:** Issue 9

### Definition of done

- [ ] Roxygen on every exported function (Hernán & Robins prose voice — named-example-first, explicit assumptions, no promotional language).
- [ ] `README.md` with quick start (5–10 line code example).
- [ ] Vignette: basic use (Feynman + worked-example-first).
- [ ] `NEWS.md` with v1 changelog.
- [ ] `R CMD check` clean (0 errors, 0 warnings, 0 notes — or NOTES documented).

(Methodology walkthrough vignette already exists separately.)

### Files

- `R/*.R` (roxygen)
- `README.md`
- `vignettes/basic-use.Rmd` (new)
- `NEWS.md`

---

## Deferred (post-v1, not part of this milestone)

- Benchmark vs `gfoRmula` on shared DGP. `dev/benchmark/`.
- v2 features per `dev/API_v1_LOCKED.md` (Z_k generalization, time-varying covariates, composite, controlling direct effect, marginal net risk).
