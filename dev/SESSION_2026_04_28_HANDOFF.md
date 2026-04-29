# Refactor session 2026-04-28 — handoff

## Status

**All refactor commits applied and tested. 218/218 PASS throughout.**

Working tree dirty (no commits yet). Tomorrow: review the diff and commit in
sensible batches.

---

## What landed

### File restructure (C5–C11)

Final `R/` layout (15 files, grouped):

```
1 · CORE (classical, reusable)
   hazards.R          predict_hazard_under, cumprod_survival, fit_logistic,
                      fit_hazard_models
   propensity.R       fit_propensity (incl. stabilization numerator model)
   weights.R          ipw, ipw_cens, ipw_static_trt, ipw_time_varying_trt (v2),
                      apply_weight_truncation, summarize_weights
   gformula_core.R    classical g-formula CIF
   ipw_core.R         cum_inc_from_weighted, weighted_hazard_by_k
                      (vector signatures, framework-agnostic)

2 · SEPARABLE OVERLAY
   separable_arms.R   separable_arm_hazards (per-arm hazard prep + Y-lag)
   separable_swap.R   swap_d_weights, swap_y_weights
   separable_ipw.R    ipw_estimate (orchestrator), estimate_weighted_cum_inc,
                      weighted_arm_cum_inc

3 · INFERENCE
   bootstrap.R, contrasts.R

4 · VALIDATION & DATA PREP
   validate.R, data_prep.R

5 · API SURFACE
   separable_effects.R (entry), accessors.R, risk_table.R, print.R, plot.R
```

Deleted: `R/ipw.R` (old monolith, ~636 LOC dismantled).

### Renames

- `causal_cr()` → `separable_effects()`
- `fit_causal_cr()` → `fit_separable_effects()`
- S3 class `causal_cr` → `separable_effects`
- `fit_one_model()` → `fit_logistic()`
- `ipcw()` → `ipw_cens()` (matches `ipw_*` suffix-by-mechanism convention)
- `add_d_swap_weights` → `swap_d_weights` (+ new `swap_y_weights`)
- `add_d_hazards`/`add_y_hazards` → unified into `separable_arm_hazards()` in
  separable_arms.R; primitives `predict_hazard_under()` + `cumprod_survival()`
  in core hazards.R
- Test file: `test-causal-cr-structure.R` → `test-separable-effects-structure.R`

### W_A — propensity weight (closes Stensrud 2020 IPW gap)

Previously missing. Now end-to-end:

- `fit_propensity()` in `R/propensity.R` — fits both full (`A ~ L`) and
  numerator (`A ~ 1`) models when `stabilize = TRUE`. Reuses `fit_logistic()`.
- Wired through `fit_separable_effects()` only when `"ipw"` is requested
  (option iii — saves the propensity fit when only g-formula runs).
- `ipw_static_trt()` in `weights.R` — point-treatment IPTW. Uses per-subject
  scalar broadcast (no cumprod-times-1 trick). Stabilized via `model_num`.
- W_A is now a **column** on pt_data (`w_a_raw` / `w_a`), so it goes
  through `apply_weight_truncation()` and `summarize_weights()` like other
  source weights.
- `weighted_arm_cum_inc()` reads `d$w_a` from the column (symmetric with
  `w_cens`) instead of computing inline post-subset.

### Public API — single `truncate` arg (Issue 2 simplification)

```r
separable_effects(pt_data, ...,
                  truncate = c(0.01, 0.99),   # default; NULL opts out
                  censoring_weights = TRUE)
```

**Replaces** the previous `extreme_weight_adjust` + `extreme_weight_threshold`
pair. Behavior changes:

- **Symmetric truncation** (Cole & Hernán 2008 AJE), not upper-tail-only.
- **Trim removed.** Truncate is the only adjustment.
- Old default `extreme_weight_threshold = 0.999` gone — new default
  `c(0.01, 0.99)`.

`fit$truncate` slot replaces `fit$extreme_weight_adjust` /
`fit$extreme_weight_threshold` slots. Bootstrap, print, accessors all
updated.

### Internal architecture — three-layer weights

```
Layer 1  Primitives          predict_hazard_under, cumprod_survival
Layer 2  Weight constructors ipw (core), ipw_cens, ipw_static_trt,
                              ipw_time_varying_trt (v2 stub),
                              swap_d_weights, swap_y_weights
Layer 3  Post-processors     apply_weight_truncation, summarize_weights
```

`ipw()` core: per-row inverse of cumulative probability of observed
history. Stabilization is structural (two-cumprod ratio), not
multiplicative.

### Truncation/summary scope

Both `apply_weight_truncation()` and `summarize_weights()` now operate on
**all six source weight columns**:

```r
c("w_cens", "w_a", "w_d_arm_10", "w_d_arm_01", "w_y_arm_10", "w_y_arm_01")
```

Previously missing: `w_a` (didn't exist as column), `w_y_arm_10/01` (added
in C7+ but never wired into truncation/summary lists).

---

## Remaining work (post-refactor)

Per `dev/REFACTOR_PLAN.md` and the locked plan:

### Issue 4 — long-format `$risk` slot
Cumulative incidence currently per-method, wide-format data.frames.
Locked API says long format with `(k, arm, cum_inc, method)` columns.
Touches gformula, ipw, contrasts, accessors, plot.

### Issue 5 — bootstrap rewrite
Cluster-by-id (already done), but signature simplification + accessor
contracts (CIs only shown when bootstrap attached) still pending.

### Issue 6 — accessor + S3 surface
`risk()`, `contrast()`, `plot()`, `summary()`, `print()`, `assumptions()`,
`risk_table()` per locked API contracts.
- `summary()` shows contrasts only when bootstrap attached
- `assumptions()` prints identifying assumptions block (new accessor)

### Issue 7 — strict validation checklist
Per locked plan, full strict validation in `validate.R`.

### Issue 8 — test rebuild
- Delete obsolete (already done for some)
- Build numerical-correctness DGP tests (10, high-leverage)
- Method × ipcw matrix (8)
- IPW Rep 1 ≡ Rep 2 under full isolation (sanity)
- Bootstrap CI coverage on simulated data
- Target: 128–160 tests across 10 files

### Issue 9 — roxygen polish + README + vignette
- Hernán-style roxygen on all exports
- README with quick start
- One user-facing vignette (Stensrud 2020 walkthrough already exists in blog)

### Stage 11 (deferred) — benchmark vs gfoRmula
On shared DGP, post-implementation. Document any differences.

### v2 deferrals (out of scope)
- Z_k generalization (Stensrud 2021)
- Time-varying covariates / treatment
- Composite, controlling-direct, marginal-net estimands
- ICE g-formula
- Causal subgroup analysis (TODO note: same machinery as separable mediation,
  not yet formalized in literature — easy add when separable code stabilizes)

---

## Tomorrow — start here

1. **Review the diff** — `git diff` and `git status` to inspect everything
   from today's session.
2. **Commit in sensible batches** — suggested groupings:
   - Phase A renames (3 file moves, no content change)
   - Phase B placeholders (7 empty stubs)
   - Hazard primitives + separable_arm_hazards (C5)
   - IPCW + ipw() core (C6)
   - Swap weights (C7)
   - IPW core math relocation (C8)
   - Rep 1/Rep 2 dispatcher merge (C9)
   - fit_logistic rename
   - causal_cr → separable_effects rename
   - W_A end-to-end (propensity model + wiring)
   - ipw_cens / ipw_static_trt / ipw_time_varying_trt
   - C10 (truncation/summary move + extensions + symmetric truncation +
     drop trim)
   - Public API: `truncate` arg
   - C11 (final relocation + delete ipw.R)
3. **Pick next issue** — recommend Issue 4 (long-format `$risk`) since it
   ripples through every accessor and unlocks Issue 6.

State to remember: `tests/manual/test_session_refactor.R` exists but isn't
run by `test_dir()`. May need pruning.
