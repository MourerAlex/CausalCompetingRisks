# SPEC.md — `CausalCompetingRisks`: Discrete-Time Separable Effects for Competing Events

## Overview

`CausalCompetingRisks` is an R package implementing discrete-time pooled logistic regression estimators for total, direct, and separable effects in the presence of competing events, following Stensrud et al. (2020, JASA) and Stensrud et al. (2021, Lifetime Data Analysis).

**Scope:** Discrete-time estimation only. No existing R package implements separable effects estimation in any form (discrete or continuous time). Martinussen & Stensrud (2023, Biometrics) develop continuous-time estimators theoretically, but these have not been packaged. This package fills the discrete-time gap; continuous-time Cox-based estimation is a potential future extension.

**Target users:** Epidemiologists and biostatisticians working with time-to-event data where competing events are present.

**Relation to future packages:** A companion package, `CausalSurvivalMediation`, will extend the separable effects framework to mediation analysis in survival settings (illness-death models), building on the theory of Didelez (2019) and extended by Breum et al. (2024). The two packages share underlying methodology (pooled logistic models, cloning, cumulative incidence computation, bootstrap). For now, `CausalCompetingRisks` is built as a self-contained package with clean internal modularity. Once the mediation package is underway, genuinely shared functions will be extracted into a lightweight core package. This avoids premature abstraction.

---

## Theoretical Background

### The Problem

In survival analysis with competing events, a treatment `A` may affect both the event of interest `Y` (e.g., prostate cancer death) and a competing event `D` (e.g., death from other causes). Classical approaches (cause-specific hazards, subdistribution hazards) conflate these pathways.

### Separable Effects Framework

The key insight is to decompose a single treatment `A` into two hypothetical components:

- `A_Y`: the component affecting the event of interest `Y`
- `A_D`: the component affecting the competing event `D`

This decomposition defines a hypothetical 4-arm trial with regimes `(a_Y, a_D) ∈ {(0,0), (0,1), (1,0), (1,1)}`. The **separable direct effect** compares arms that differ only in `a_Y` while holding `a_D` fixed.

### Estimands

Given cumulative incidence function `F_Y(t; a_Y, a_D)`:

1. **Total effect:** `F_Y(t; 1, 1) - F_Y(t; 0, 0)` — the overall effect of treatment on `Y`, allowing both pathways.
2. **Direct effect ("zombie" / eliminative):** `F_Y(t; 1, elim) - F_Y(t; 0, elim)` — effect on `Y` under hypothetical elimination of `D` (set `h_D = 0`). Called the "zombie" effect because subjects who would have died from `D` are kept alive (as zombies) to remain at risk for `Y`.
3. **Indirect effect ("killer"):** Total minus direct — the component of the total effect that operates through `D`. Treatment "kills off" subjects via the competing event, removing them from the `Y` risk pool.
4. **Separable direct effect (A_Y):** `F_Y(t; 1, 0) - F_Y(t; 0, 0)` — effect through the `A_Y` pathway only, while `A_D = 0`.
5. **Separable indirect effect (A_D):** `F_Y(t; 1, 1) - F_Y(t; 1, 0)` — effect through the `A_D` pathway only. Equals total minus separable direct.
6. **Composite endpoint:** `F_{Y∨D}(t; 1, 1) - F_{Y∨D}(t; 0, 0)` — total effect on the composite of `Y` or `D` (whichever occurs first). Useful as a clinical summary when both events are adverse.

The cross-arm regime `(a_Y = 1, a_D = 0)` is the distinctive contribution of the separable framework: it requires predicting `Y`-hazards under one treatment assignment while using `D`-hazards from another.

### Two Estimation Methods

#### G-formula (parametric)
- Fit pooled logistic models for `h_Y` and `h_D` on person-time data
- Create cloned datasets for each arm configuration `(a_Y, a_D)`
- Predict hazards in each clone, compute cumulative incidence by averaging over subjects
- Censoring is handled by conditioning: covariates that predict censoring are included in the hazard models; no separate censoring model is needed

#### IPW (inverse probability weighting)
- Fit pooled logistic models for the competing event hazard and censoring hazard
- Construct weights that reweight observed data to mimic the target intervention
- For the separable effect: `w_sep = w_D / w_C` where `w_D` reweights the `D`-hazard across arms and `w_C` adjusts for censoring
- Requires an explicit censoring model (unlike g-formula)
- Supports both stabilized and unstabilized weights

---

## Package Design

### User-Facing API

```r
# Primary fitting function
fit <- causal_cr(
  data,                          # data.frame: person-level (wide) data
  id = "id",                     # character: subject identifier column
  time = "time",                 # character: event/censoring time column
  event_y = "y_event",           # character: primary event indicator
  event_d = "d_event",           # character: competing event indicator
  censor = "c_event",            # character: censoring indicator (optional for g-formula)
  treatment = "A",               # character: treatment column
  covariates = c(...),           # character vector: baseline covariates
  time_varying = c(...),         # character vector: time-varying covariates (optional)
  method = "gformula",           # "gformula" or "ipw"
  times = 1:K,                   # integer vector: discrete time points
  formulas = NULL,               # optional: named list of user-specified model formulas
  stabilized = TRUE,             # logical: use stabilized weights (IPW only)
  bootstrap = FALSE,             # logical: compute bootstrap CIs
  n_boot = 500                   # integer: number of bootstrap samples
)

# Methods
summary(fit)                     # print effect estimates at final time point
print(fit)                       # print object
plot(fit)                        # cumulative incidence curves for all arms
confint(fit)                     # bootstrap confidence intervals

# Extract components
coef(fit)                        # named vector of effect estimates
cumulative_incidence(fit)        # data.frame of cumulative incidence by arm and time
hazards(fit)                     # fitted hazard objects
contrasts(fit)                   # risk differences and risk ratios
```

### Return Object

`causal_cr()` returns an S3 object of class `"causal_cr"` containing:

```r
list(
  cumulative_incidence = data.frame(
    time, arm_11, arm_00, arm_10, arm_01  # F_Y(t; a_Y, a_D) for each regime
  ),
  contrasts = data.frame(
    time, total_rd, total_rr,
    direct_rd, direct_rr,             # "zombie" / eliminative
    indirect_rd, indirect_rr,         # "killer" (total minus direct)
    separable_ay_rd, separable_ay_rr,
    separable_ad_rd, separable_ad_rr,
    composite_rd, composite_rr        # composite endpoint (Y or D)
  ),
  models = list(
    model_y = <glm>,             # fitted Y-hazard model (predicts h_y)
    model_d = <glm>,             # fitted D-hazard model (predicts h_d)
    model_c = <glm or NULL>,     # fitted censoring model (IPW only)
    model_c_num = <glm or NULL>  # numerator censoring model (stabilized IPW)
  ),
  person_time = <data.frame>,    # the expanded person-time dataset used
  bootstrap = list(              # NULL if bootstrap = FALSE
    estimates = <array>,         # (n_times x n_estimands x n_boot) array
    ci_lower = <data.frame>,
    ci_upper = <data.frame>
  ),
  call = <call>,
  method = "gformula" | "ipw",
  times = <integer vector>,
  n = <integer>                  # number of subjects
)
```

### Internal Function Architecture

```
R/
├── causal_cr.R             # Main entry point: causal_cr()
├── data_prep.R             # to_person_time(), validate_input()
│                           #   - survSplit-based expansion from wide to long
│                           #   - event indicator construction (Y, D, C columns)
│                           #   - NA handling for censored/competing rows
├── hazard_models.R         # fit_hazard_y(), fit_hazard_d(), fit_censoring()
│                           #   - Pooled logistic (glm, family=binomial)
│                           #   - Default formulas: A * (time + time^2 + time^3) + covariates
│                           #   - User-overridable formula interface
├── gformula.R              # gformula_estimate()
│                           #   - Clone datasets for each arm (a_Y, a_D)
│                           #   - Predict hazards in each clone
│                           #   - Compute per-subject cumulative survival
│                           #   - Average cumulative incidence across subjects
├── ipw.R                   # ipw_estimate()
│                           #   - Compute censoring weights (stabilized/unstabilized)
│                           #   - Compute competing event weights for separable effects
│                           #   - Weighted non-parametric cumulative hazard
│                           #   - Weighted Aalen-Johansen cumulative incidence
├── cumulative_risk.R       # compute_cumincidence_parametric() [g-formula]
│                           #   compute_cumincidence_nonparametric() [IPW]
│                           #   - Iterative form: F(t) = sum_{s<=t} h_Y(s) * (1-h_D(s)) * S(s-1)
│                           #   - Competing event variant: F_D(t) = sum_{s<=t} h_D(s) * S(s-1)
├── contrasts.R             # compute_contrasts()
│                           #   - Total effect: arm(1,1) - arm(0,0)
│                           #   - Direct effect ("zombie"): arm(1,elim) - arm(0,elim) where h_D=0
│                           #   - Indirect effect ("killer"): total - direct
│                           #   - Separable A_Y effect: arm(1,0) - arm(0,0) [or arm(1,1) - arm(0,1)]
│                           #   - Separable A_D effect: total - separable A_Y
│                           #   - Composite endpoint: F_{Y∨D}(t;1,1) - F_{Y∨D}(t;0,0)
│                           #   - Risk differences and risk ratios for all
├── bootstrap.R             # bootstrap_causal_cr()
│                           #   - Subject-level resampling (resample IDs, expand to person-time)
│                           #   - Percentile confidence intervals
├── methods.R               # S3 methods: print, summary, plot, coef, confint
├── utils.R                 # safe_cumprod(), validate_formulas(), etc.
└── data.R                  # Documentation for bundled example dataset
```

### Data Flow

```
User provides wide (person-level) data
         │
         ▼
    to_person_time()          # Expand to long (person-time) format
         │                    # Construct event indicators: y_event, d_event, c_event
         │                    # Handle NAs: censor rows get NA for Y and D
         │                    #             D-event rows get NA for Y
         ▼
    fit_hazard_y()            # glm(y_event ~ A_y * time_terms + covariates)
    fit_hazard_d()            # glm(d_event ~ A_d * time_terms + covariates)
    [fit_censoring()]         # IPW only: glm(c_event ~ A * time_terms + covariates)
         │
         ▼
    ┌─── method? ───┐
    │               │
  g-formula        IPW
    │               │
    ▼               ▼
 Clone datasets   Compute weights
 for each arm     (censoring + competing event)
    │               │
    ▼               ▼
 Predict hazards  Weighted non-parametric
 in each clone    cumulative hazards
    │               │
    ▼               ▼
 Per-subject      Aalen-Johansen
 cum. incidence   cum. incidence
    │               │
    └───────┬───────┘
            ▼
      compute_contrasts()     # RD, RR for total/direct/separable
            │
            ▼
      Return "causal_cr" object
```

---

## Key Implementation Details

### The Cross-Arm Prediction Trick (G-formula)

This is the core mechanism for separable effects. For arm `(a_Y=1, a_D=0)`:

1. Create a cloned dataset where `A_Y = 1` (for predicting Y-hazard) and `A_D = 0` (for predicting D-hazard)
2. In the reference code, this is achieved by having two copies of the treatment variable: `rx` (maps to `A_Y`) and `Orx` (maps to `A_D`)
3. Predict `hazardP` using the Y-model with `rx=1`, predict `hazardO` using the D-model with `Orx=0`
4. Compute survival: `s = (1 - hazardP) * (1 - hazardO)`

The package must support this by either:
- Automatically creating the duplicated treatment column, or
- Allowing users to specify which treatment variable maps to which hazard model

### The Cross-Arm IPW Weights

For IPW separable effects, the key weight for arm `(a_Y=0, a_D=1)`:

```
w_sep(t) = [∏_{s≤t} (1 - ĥ_D(s | A=0))] / [∏_{s≤t} (1 - ĥ_D(s | A=1))]
```

This reweights subjects in the treated arm (`A=1`) so their competing event hazard looks like the untreated arm. The 2021 paper also uses an *inverse* of this weight for the complementary separable effect.

Implementation note: The IPW censoring model uses observed `A` (not overridden to `a_Y`), because censoring weights apply uniformly within each observed arm.

### Cumulative Incidence Computation

**G-formula (parametric, per-subject):**
```
F_Y(t) = (1/n) * Σ_i Σ_{s≤t} ĥ_Y(s|i) * (1 - ĥ_D(s|i)) * Ŝ(s-1|i)
```
where `Ŝ(s|i) = ∏_{r≤s} (1 - ĥ_Y(r|i)) * (1 - ĥ_D(r|i))`

This is an *iterative* cumulative sum averaged over subjects.

**IPW (non-parametric, weighted Hájek/Horvitz-Thompson):**
```
ĥ_Y(t) = Σ w_i * I(Y_i=1, T_i=t) / Σ w_i * I(T_i≥t, not censored, D not yet occurred)
```
Then plug into the Aalen-Johansen formula.

### Person-Time Data Construction

The package must handle:
1. `survSplit()`-based expansion from wide to person-time format
2. Proper event indicator construction:
   - Censored rows: `y_event = NA`, `d_event = NA`
   - Competing event rows: `y_event = NA` (removed from Y risk set)
   - Normal rows at risk: `y_event = 0/1`, `d_event = 0/1`
3. `group_by(id)` for all `cumprod()` operations to prevent cross-subject contamination

### Handling Censoring

**G-formula:** No separate censoring model. Censoring is handled by conditioning — include covariates that predict censoring in the Y and D hazard models. This works because the g-formula conditions on the full covariate history.

**IPW:** Requires an explicit censoring model: `glm(c_event ~ A * time_terms + covariates)`. The censoring weight is `w_C(t) = 1 / ∏_{s≤t} (1 - ĥ_C(s))`.

### Bootstrap

- Subject-level resampling: sample subject IDs with replacement, then expand to person-time
- The 2021 supplementary code shows the pattern: create `new_id` for bootstrap individuals, refit all models, recompute all estimates
- Store results in a 3D array: `(n_times × n_estimands × n_boot)`
- Percentile confidence intervals at each time point

---

## Bundled Data

Include a cleaned version of the prostate cancer dataset used in the reference code (from `byar` data in the `Hmisc` package). This is publicly available and well-documented. The package vignette should reproduce the reference results.

---

## Dependencies

**Hard dependencies (Imports):**
- `survival` (for `survSplit()`)
- `stats` (for `glm()`, `predict()`, `binomial()`)

**Suggested:**
- `ggplot2` (for plot method)
- `boot` (alternative bootstrap interface)

Minimize dependencies. The reference code uses `Hmisc`, `splines`, `dplyr`, `zoo` — most of this can be replaced with base R.

---

## Testing Strategy

1. **Reproduce reference results:** The primary validation is that `causal_cr()` with the prostate data matches the published results from Rojas-Saunero et al. and Stensrud et al.
2. **Unit tests per function:** `to_person_time()` produces correct dimensions and indicators; `compute_cumincidence()` matches hand calculation on toy data; weights are positive and finite.
3. **Edge cases:** Single time point, no competing events observed, no censoring, all subjects censored at some time.

---

## Future Extensions (Out of Scope for v1)

- Continuous-time (Cox) estimation (theory in Martinussen & Stensrud 2023, Biometrics)
- Time-varying treatments
- Mediation analysis / illness-death models (separate package: `CausalSurvivalMediation`, theory of Didelez 2019, extended by Breum et al. 2024)
- Doubly robust / TMLE / one-step estimators (influence-function-based, per Martinussen & Stensrud 2023)
- Sensitivity analysis for separability conditions
