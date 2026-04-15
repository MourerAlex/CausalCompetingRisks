# SPEC_CC.md — `CausalCompetingRisks` v1: Consolidated Specification

**Date:** 2026-04-15
**Status:** Agreed architecture from interactive review session.
**Supersedes:** SPEC.md for decisions where they differ.

---

## 1. API

Single entry point: `causal_cr()`. Returns S3 object of class `"causal_cr"`.

```r
fit <- causal_cr(
  data,
  id         = "id",
  time       = "time",
  event      = "event_type",
  treatment  = "A",
  covariates = c(...),
  event_y    = NULL,           # NULL = 0/1/2 convention (integer only)
  event_d    = NULL,           # character event -> must specify all three
  event_c    = NULL,
  method     = "both",         # "both", "gformula", "ipw1", "ipw2"
  times      = NULL,           # NULL = 12 equally-spaced intervals
  eval_times = NULL,           # time points for contrasts table
  formulas   = NULL,           # named list: override model formulas
  time_varying = NULL,         # not implemented v1 (error if used)
  extreme_weight_adjust = "truncate",   # "truncate", "trim", "none"
  extreme_weight_threshold = 0.999,
  bootstrap  = FALSE,
  n_boot     = 500,
  alpha      = 0.05
)
```

### Methods

```r
# S3 methods
print(fit)
summary(fit)
confint(fit)                  # bootstrap CIs at eval_times

# Accessor functions (return sub-objects with their own plot/print methods)
risk(fit)                     # cumulative incidence curves
contrast(fit)                 # effect estimates over time
diagnostic(fit)               # weight distributions, isolation check

# Plot methods (all ggplot2-based)
plot(risk(fit))               # cumulative incidence curves (3 arms)
plot(contrast(fit))           # effect-over-time curves
plot(diagnostic(fit))         # weight distribution plots
```

---

## 2. Input Validation

### Hard errors (stop):
- Column doesn't exist in data
- Treatment not binary (0/1)
- Event column not exactly 3 unique values
- Character event without event_y/event_d/event_c all specified
- Unknown method
- Duplicate rows
- Duplicate (id, time) pairs in person-time data
- time < 0
- time_varying not NULL (v1 placeholder)

### Warnings (continue):
- Covariates with NAs
- Covariate with > 20 categories (likely miscoded continuous)
- Few events (< 5 for Y or D)
- time == 0 events -> shift with tstart = -1e-7
- Continuous (non-integer) times -> auto-discretize with message

### Data format detection:
- `n_distinct(id) == nrow(data)` -> subject-level (wide), expand internally via survSplit
- `n_distinct(id) < nrow(data)`, no duplicate (id, time) -> person-time (long), use as-is

### Event specification:
- Integer 0/1/2 + event_y/d/c all NULL -> convention (0=cens, 1=Y, 2=D), print informative message
- Character event -> must specify event_y, event_d, event_c explicitly, else error
- Always enforce exactly 3 unique values

---

## 3. Time Handling

- `times = NULL` -> 12 equally-spaced intervals from 0 to max event time
- `times = integer scalar` -> that many equally-spaced intervals
- `times = numeric vector` -> use as cut points directly
- time == 0 edge case: set tstart = -1e-7 with warning
- Continuous data: auto-discretize using `times` argument as cut points

### Internal time convention:
- `tstart = 0` for everyone (month 0 = enrollment)
- Events shifted +1 so month 1 = first interval at risk
- `k = tstart` as interval index
- `survSplit` with `event_indicator = 1L` for ALL subjects (including censored)

---

## 4. Person-Time Construction

Single `survSplit()` call with `event_indicator = 1L` for all subjects. This is critical:
censored subjects need their terminal row flagged to derive `c_event`.

After expansion, derive indicators from `event_type`:
```r
pt_df <- pt_df %>%
  mutate(
    y_event = ifelse(event == 1 & event_type == 1, 1, 0),
    d_event = ifelse(event == 1 & event_type == 2, 1, 0),
    c_event = ifelse(event == 1 & event_type == 0, 1, 0)
  )
```

### Temporal ordering (C_k, D_k, Y_k):
```r
pt_df <- pt_df %>%
  mutate(
    d_event = ifelse(c_event == 1, NA, d_event),
    y_event = ifelse(c_event == 1, NA, y_event),
    y_event = ifelse(d_event == 1, NA, y_event)
  )
```

Censored -> both D and Y become NA. D occurred -> Y becomes NA.
NAs define risk sets for `glm()` via `na.action`.

### Treatment copies:
```r
pt_df$A_y <- pt_df$A
pt_df$A_d <- pt_df$A
```
Equal in observed data; diverge in cloned datasets for cross-arm prediction.

---

## 5. Hazard Models

### Default formula:
```r
event ~ A + k + I(k^2) + I(k^3) + covariates
```
Additive (no interaction). User can override via `formulas` argument.

### Which models for which method:
| Method | Models fitted |
|--------|--------------|
| gformula | fit_y_haz (uses A_y), fit_d_haz (uses A_d) |
| ipw1 | fit_d_haz (uses A_d), fit_cens |
| ipw2 | fit_y_haz (uses A_y), fit_cens |
| both | all three: fit_y_haz, fit_d_haz, fit_cens |

`method = "both"` is the default because all sensitivity analysis checks require IPW representations.

---

## 6. G-Formula Estimation

1. Fit `fit_y_haz` (with `A_y` as treatment) and `fit_d_haz` (with `A_d` as treatment) on observed person-time data
2. Create cloned datasets for each arm via `make_clone()`:
   - arm (1,1): `A_y = 1, A_d = 1` (observed treated)
   - arm (0,0): `A_y = 0, A_d = 0` (observed control)
   - arm (1,0): `A_y = 1, A_d = 0` (separable — cross-arm)
3. `predict_hazards()`: predict `haz_y` from `fit_y_haz` (sees `A_y`) and `haz_d` from `fit_d_haz` (sees `A_d`)
4. `compute_cum_inc()`: per-subject cumulative incidence, then average across subjects

### Cumulative incidence formula:
```
F_Y(t) = (1/n) * sum_i sum_{s<=t} h_Y(s|i) * (1 - h_D(s|i)) * S(s-1|i)
```
where `S(s|i) = prod_{r<=s} (1 - h_Y(r|i)) * (1 - h_D(r|i))`

**Critical**: Per-subject cumprod then average (NOT average then cumprod). The cumprod is nonlinear and does not commute with the mean. The code order (mean by k, then cumsum) is equivalent to (cumsum per subject, then mean) by commutativity of the double sum — but the cumprod must happen per-subject first. See MEMO_part5_math_traceback.md.

### Censoring in g-formula:
No separate censoring model. Censoring is handled by conditioning — include predictors of censoring as covariates in the Y and D hazard models.

---

## 7. IPW Estimation

### IPW Representation 1:
- Stand in `A = a_Y` arm
- `W_D` reweights D-hazard distribution across arms
- `W_C` corrects for censoring
- Needs: D model + censoring model (NOT Y model)

### IPW Representation 2:
- Stand in `A = a_D` arm
- `W_Y` reweights Y-hazard distribution
- Needs: Y model + censoring model (NOT D model)

### Variable naming:
- `d_free_k_ad` / `d_free_k_ay`: interval-level D-free probabilities (1 - h_D at each k)
- `cs_surv_d_ad` / `cs_surv_d_ay`: cumulative cause-specific survival (cumprod of D-free probs)
  - These are NOT marginal survival — they are P(D hasn't occurred by k | event-free up to k)
- `w_d`: ratio of cause-specific survivals across arms
- `w_cens`: inverse probability of censoring weight
- `w_total`: combined weight

### Censoring model:
Uses **observed treatment A** (not overridden to a_Y), because censoring weights correct for the censoring mechanism as it actually operated.

---

## 8. Weight Handling

- `extreme_weight_adjust = "truncate"` (default): cap at percentile threshold, keep everyone
- `extreme_weight_adjust = "trim"`: remove subjects beyond threshold
- `extreme_weight_adjust = "none"`: no adjustment
- One or the other, not both
- Flagged/affected IDs stored in `$diagnostics`
- Stabilized weights: deferred to v1.1

---

## 9. Return Object

S3 class `"causal_cr"`:

```r
list(
  cumulative_incidence = data.frame(
    k, arm_11, arm_00, arm_10     # F_Y(t; a_Y, a_D) at each time point
  ),
  contrasts = data.frame(
    k,
    total_rd, total_rr,
    sep_direct_rd, sep_direct_rr,       # separable direct (A_Y): arm(1,0) - arm(0,0)
    sep_indirect_rd, sep_indirect_rr    # separable indirect (A_D): arm(1,1) - arm(1,0)
  ),
  models = list(
    model_y = <glm>,           # Y-hazard model
    model_d = <glm>,           # D-hazard model
    model_c = <glm or NULL>    # censoring model (IPW only)
  ),
  person_time = <data.frame>,  # expanded person-time data
  diagnostics = list(
    weight_summary = <data.frame>,    # mean, median, min, max, percentiles of weights
    truncated_ids  = <integer>,       # IDs with truncated/trimmed weights
    isolation_check = <list>          # W_D, W_Y departure from 1
  ),
  bootstrap = NULL | list(
    replicates = <array>,      # n_boot x n_arms x n_times (full draws)
    ci_curves  = <data.frame>, # CI bands on cumulative incidence
    ci_contrasts = <data.frame> # CI on RD/RR at eval_times
  ),
  call     = <call>,
  method   = <character>,
  times    = <numeric vector>,
  n        = <integer>
)
```

### Contrasts (v1):
- Risk differences (additive) AND risk ratios (multiplicative)
- For: total, separable direct (A_Y), separable indirect (A_D)
- At eval_times (default: all time points)

### Deferred contrasts (v1.1):
- Direct effect ("zombie", h_D = 0)
- Indirect effect ("killer", total minus zombie)
- Composite endpoint

---

## 10. Accessor Functions

### `risk(fit)`
Returns a `"causal_cr_risk"` object with cumulative incidence data.
`plot()` method produces ggplot2 cumulative incidence curves.

### `contrast(fit)`
Returns a `"causal_cr_contrast"` object with effect estimates.
`plot()` method produces effect-over-time curves.

### `diagnostic(fit)`
Returns a `"causal_cr_diagnostic"` object with weight summaries and isolation checks.
`plot()` method produces weight distribution plots.

---

## 11. Plots

All plots use ggplot2. Colorblind-safe Okabe-Ito palette.

### Cumulative incidence plot — `plot(risk(fit))`:
- 3 curves: arm(1,1), arm(0,0), arm(1,0)
- Risk table below plot: at risk, events, censored by arm (unweighted)
- Colors: black for reference arms, distinguished for cross-arm

### Effect-over-time plot — `plot(contrast(fit))`:
- Colors:
  - Total effect: #000000 (black)
  - Separable direct (A_Y): #009E73 (green, Okabe-Ito)
  - Separable indirect (A_D): #D55E00 (vermillion, Okabe-Ito)
- Reference line at 0 (for RD) or 1 (for RR)
- CI ribbons when bootstrap available

### Diagnostic plots — `plot(diagnostic(fit))`:
- Weight distributions (histogram/density)
- W_D and W_Y departure from 1 over time (isolation check)

### Plot customization:
- Title, subtitle, legend handled via ggplot2 (users add `+ labs(...)` etc.)
- Risk table below main plot (at risk, censored, events by arm — unweighted for v1)

---

## 12. Bootstrap

- **Resampling**: Subject-level (resample IDs with replacement, pull all person-time rows)
- **Per replicate**: Re-fit all models, re-run estimation, store cumulative incidence + contrasts
- **Storage**: 3D array `n_boot x n_arms x n_times` for cumulative incidence; contrasts matrix `n_boot x n_contrasts x n_eval_times`
- **CIs**: Percentile method (quantiles of bootstrap distribution)
- **Defaults**: `n_boot = 500`, `alpha = 0.05`
- **Implementation**: Simple `lapply` loop (no parallelization in v1)
- **Progress**: `message()` with counter (e.g., "Bootstrap replicate 50/500")
- **Duplicate ID handling**: Append replicate suffix to create unique IDs for resampled subjects
- **Full replicates stored** in `$bootstrap$replicates` for custom analysis
- **Integration**: `confint()` prints CI table; `plot(risk(fit))` and `plot(contrast(fit))` auto-add CI ribbons when bootstrap exists

---

## 13. Sensitivity Analysis

| Type | What | How | Version |
|------|------|-----|---------|
| Model SA | Models correctly specified? | G-formula vs IPW1 vs IPW2 comparison | v1 |
| Weight diagnostics | Weights well-behaved? | Inspect W_D, W_Y, W_C distributions | v1 |
| Isolation check | Which isolation level? | W_D, W_Y departure from 1 -> full/partial/none | v1 |
| t_k sensitivity | Robust to dismissible violations? | Re-estimate under grid of t_k values | v1 |
| Z_k partition check | Proposed covariate partition valid? | Compute W_LD, W_LY with user partition | v2 |
| Isolation adjustment | Time-varying common causes? | Compare with/without W_LD, W_LY | v2 |

### Isolation hierarchy (from 2021 paper):
1. **Full isolation** (W_D ~ 1, W_Y ~ 1): strongest, separable effects cleanly identified
2. **Partial D-isolation** (W_D ~ 1, W_Y != 1): D-model correct under both arms
3. **Partial Y-isolation** (W_Y ~ 1, W_D != 1): Y-model correct under both arms
4. **No isolation + Z_k partition**: still identifiable with additional W_LD, W_LY weights
5. **No isolation, no partition**: not identifiable (Lemma 4)

All SA checks require IPW representations -> reason `method = "both"` is the default.

### Positivity / overlap:
- G-formula hides positivity problems (extrapolates silently)
- IPW makes them explicit (extreme weights)
- Running both methods: IPW weights serve as diagnostic for g-formula reliability
- No independent g-formula positivity diagnostic exists

---

## 14. Bundled Data

Ship the Byar & Green (1980) prostate cancer dataset from `Hmisc::byar`. Public domain clinical trial data — no license restriction.

```r
data(prostate, package = "Hmisc")
# Clean, subset to DES + placebo, save as prostate_data
usethis::use_data(prostate_data)
```

Also provide internal simulated data generator for unit tests (controlled scenarios with known true effects).

---

## 15. Dependencies

### Imports (hard):
- `survival` (survSplit)
- `stats` (glm, predict, binomial)
- `ggplot2` (plot methods — core feature, not optional)

### Suggests:
- `Hmisc` (for rebuilding bundled data from source)

### Explicitly NOT dependencies:
- `dplyr`: package uses base R internally
- `boot`: custom bootstrap loop
- `future.apply`: deferred to v1.1
- `progressr`: deferred to v1.1

---

## 16. File Structure

```
R/
  causal_cr.R           # Main entry point
  validate.R            # validate_input()
  data_prep.R           # to_person_time() — survSplit + indicator derivation
  hazard_models.R       # fit_hazard_y(), fit_hazard_d(), fit_censoring()
  gformula.R            # gformula_estimate(), make_clone(), predict_hazards(), compute_cum_inc()
  ipw.R                 # ipw_estimate(), compute_weights()
  contrasts.R           # compute_contrasts()
  bootstrap.R           # bootstrap_causal_cr()
  accessors.R           # risk(), contrast(), diagnostic()
  plot.R                # plot methods for accessor objects
  print.R               # print/summary/confint methods
  utils.R               # safe_cumprod(), internal helpers
  data.R                # dataset documentation
data/
  prostate_data.rda
man/                    # roxygen2-generated
vignettes/
tests/
  testthat/
DESCRIPTION
NAMESPACE
```

---

## 17. Version Roadmap

### v1 (now):
- G-formula + IPW1 + IPW2 for separable effects
- Total, separable direct (A_Y), separable indirect (A_D) contrasts
- Bootstrap (sequential, percentile CIs)
- ggplot2 plots with Okabe-Ito palette and risk tables
- Weight diagnostics and isolation check
- Bundled prostate data

### v1.1:
- Direct effect (zombie, h_D = 0) and indirect effect (killer)
- Stabilized weights
- Parallelized bootstrap (future.apply)
- Composite endpoint

### v2:
- time_varying covariates
- Z_k partition check (no isolation case)
- IPCW + g-formula hybrid (censor_weights)
- Arm-stratified models option

### Future:
- Continuous-time (Cox) estimation
- Doubly robust / TMLE estimators
