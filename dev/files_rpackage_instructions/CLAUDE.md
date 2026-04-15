# CLAUDE.md — Context for Building the `CausalCompetingRisks` R Package

This document summarizes key implementation knowledge, design decisions, and hard-won clarifications developed over an extended series of conversations about causal inference with competing events. It is intended as context for Claude Code when working on this package.

---

## 1. What This Package Does

`CausalCompetingRisks` implements discrete-time pooled logistic regression estimators for **separable effects** in survival analysis with competing events. The theoretical framework comes from:

- **Stensrud et al. (2020), JASA** — introduces separable effects for competing events
- **Stensrud et al. (2021), Lifetime Data Analysis** — extends to time-varying covariates and the SPRINT trial application
- **Young et al. (2020), Statistics in Medicine** — provides the g-formula and IPW implementations for total and direct effects

The package fills a genuine gap: **no existing R package implements separable effects estimation** in any form. Martinussen & Stensrud (2023, Biometrics) develop continuous-time estimators theoretically with accompanying simulation code, but this has not been packaged. The `mets` package (Scheike & Holst) handles general competing risks modeling but does not implement separable effects. The `gfoRmula` package cannot accommodate the cross-arm prediction trick required for separable effects.

### Package Architecture

`CausalCompetingRisks` is the first of two planned user-facing packages:

- **`CausalCompetingRisks`** — separable effects for competing events (this package)
- **`CausalSurvivalMediation`** — separable effects for mediation in survival settings (illness-death models, building on the theory of Didelez 2019, extended by Breum et al. 2024)

Both share underlying methodology (pooled logistic models, cloning, cumulative incidence, bootstrap). For now, `CausalCompetingRisks` is built self-contained with clean internal modularity. Once the mediation package is underway, genuinely shared functions will be extracted into a lightweight core/engine package. This avoids premature abstraction — we refactor once we have two concrete consumers, not before.

See `SPEC.md` for the full API specification and architecture.

---

## 2. The Core Conceptual Framework

### Treatment Decomposition

The central idea: decompose a single treatment `A` into two hypothetical components `A_Y` (affects event of interest) and `A_D` (affects competing event). This creates a hypothetical 4-arm trial:

| Arm | A_Y | A_D | Interpretation |
|-----|-----|-----|----------------|
| (1,1) | 1 | 1 | Full treatment (observed treated arm) |
| (0,0) | 0 | 0 | Full control (observed control arm) |
| (1,0) | 1 | 0 | Treatment's Y-pathway only |
| (0,1) | 0 | 1 | Treatment's D-pathway only |

The arms `(1,1)` and `(0,0)` are observed. The arms `(1,0)` and `(0,1)` are **counterfactual** — they require the cross-arm prediction trick.

### Three Effect Types (plus composite)

- **Total effect:** `F(t; 1,1) - F(t; 0,0)` — uses observed arms only
- **Direct effect ("zombie" / eliminative):** Set `h_D = 0` in both arms, compare `F(t; 1, elim) - F(t; 0, elim)` — hypothetical elimination of competing event. Called "zombie" because subjects who would have died from D are kept alive to remain at risk for Y.
- **Indirect effect ("killer"):** Total minus direct — the component operating through D. Treatment "kills off" subjects via the competing event, removing them from the Y risk pool.
- **Separable direct effect (A_Y):** `F(t; 1,0) - F(t; 0,0)` — the cross-arm quantity. Requires stronger assumptions (separability) but answers a more clinically meaningful question.
- **Separable indirect effect (A_D):** `total - separable_A_Y` = `F(t; 1,1) - F(t; 1,0)` — effect through the competing event pathway.
- **Composite endpoint:** `F_{Y∨D}(t; 1,1) - F_{Y∨D}(t; 0,0)` — total effect on whichever event occurs first. Useful when both events are adverse outcomes.

### Separability Conditions

The separable effects are identified under conditions that formalize the treatment decomposition. In brief:
1. The treatment can be split into components that affect Y and D through separate pathways
2. The component affecting D does not affect Y except through D
3. Standard causal conditions (consistency, exchangeability, positivity) hold for the decomposed treatment

These are **untestable** assumptions. The package should document them clearly but not attempt to verify them programmatically.

---

## 3. Reference Implementations

Three reference code files are available in the project:

### `gform.txt` (Rojas-Saunero workshop code)
- Prostate cancer dataset from `Hmisc::byar`
- G-formula implementation for total, direct, and separable effects
- Key function: `calculateCumInc()` — iterative per-subject cumulative incidence
- Treatment variable: `rx` (DES=1, placebo=0)
- Competing event: `otherDeath`; event of interest: `prostateDeath`
- For separable effects: duplicated treatment as `rx` (→ A_Y) and `Orx` (→ A_D)

### `ipw.txt` (Rojas-Saunero workshop code)
- Same prostate dataset
- IPW implementation with stabilized and unstabilized weights
- Key functions: `nonParametricCumHaz()`, `nonParametricCumInc()`, `discrete_cuminc_prost()`
- Censoring model fitted only on `dtime > 50` (data-specific: no censoring before month 50)
- For separable effects: cross-arm weights computed from predicted D-hazards under both arms

### `10985_2021_9530_MOESM1_ESM.r` (Stensrud et al. 2021 supplementary)
- SPRINT trial data (not publicly available, but code structure is the reference)
- IPW-only implementation with time-varying covariates (blood pressure)
- Key function: `bootSeparableEffectsLK()` — wraps the full estimation pipeline for bootstrap
- Arm-stratified models: separate Y-hazard models for treated and control
- Bootstrap: subject-level resampling with ID reassignment
- Uses `zoo::na.locf()` for carrying forward time-varying covariates

### `utility_functions.R`
- Contains all shared utility functions used across the workshop code
- `nonParametricCumHaz()`: weighted hazard estimation for IPW
- `nonParametricCumInc()`: Aalen-Johansen cumulative incidence from cause-specific hazards
- `calculateCumInc()`: parametric g-formula cumulative incidence (per-subject, matrix-based)
- `discrete_cuminc_prost()`: simple weighted cumulative incidence (Horvitz-Thompson style)

### Dudukina tidyverse translation (external reference)
- URL: https://www.elenadudukina.com/post/competing-events/2022-11-24-competing-events/
- Tidyverse translation of the Rojas-Saunero workshop code
- Key improvements over the base R version:
  - `calculateCumInc()` parameterizes hazard column names (`hazardO = "hazard_odeath"`, `hazardP = "hazard_pdeath"`) instead of hardcoding — closer to what the package API needs
  - Uses `group_by(patno) %>% mutate(cumprod(...))` instead of `aggregate()` — the pattern the package should follow
  - Names the duplicated treatment `odeath_rx` instead of `Orx` — more self-documenting
  - **Fixes a bug in the `calculateCumInc()` initialization**: uses `hazardP * (1 - hazardO)` at time 0 instead of the incorrect `hazardP * hazardO` in the original
- Covers IPW (direct + total effect) and g-formula (separable effects) sections
- Does NOT cover the IPW separable effects section

---

## 4. Critical Implementation Details

### 4.1 The Cross-Arm Prediction Trick (G-formula)

This is the hardest part to get right. For the separable effect arm `(a_Y=1, a_D=0)`:

```r
# Create cloned dataset
treatAy <- baseline[rep(1:n, each = length(cutTimes)), ]
treatAy$dtime <- rep(cutTimes, n)
treatAy$rx <- 1     # A_Y = 1: used by Y-hazard model
treatAy$Orx <- 0    # A_D = 0: used by D-hazard model

# Predict hazards using DIFFERENT treatment values for each model
treatAy$hazardP <- predict(plrFitP, newdata = treatAy, type = 'response')  # uses rx=1
treatAy$hazardO <- predict(plrFitO, newdata = treatAy, type = 'response')  # uses Orx=0

# Survival combines both hazards
treatAy$s <- (1 - treatAy$hazardP) * (1 - treatAy$hazardO)
```

**The `rx`/`Orx` pattern**: The Y-hazard model includes `rx` as the treatment variable. The D-hazard model includes `Orx`. Both are fit on the same observed data where `rx == Orx`. But in the counterfactual clones, they can diverge.

In the package, this should be abstracted: the user specifies one treatment column, and the package internally creates the duplicated column and manages which model sees which version.

### 4.2 Cumulative Incidence Computation (G-formula)

The `calculateCumInc()` function uses an iterative form:

```
At time 1:  CI[1,i] = h_Y(0|i) * (1 - h_D(0|i))   [non-competing; S(-1) = 1]
            CI[1,i] = h_D(0|i)                       [competing; S(-1) = 1]
At time t:  CI[t,i] = h_Y(t-1|i) * (1-h_D(t-1|i)) * S(t-2|i)   for non-competing
            CI[t,i] = h_D(t-1|i) * S(t-2|i)                      for competing
Final:      F_Y(t) = mean over subjects of cumsum(CI[,i])
```

**Important**: The `competing = FALSE` branch computes `subInputDataO * subInputDataP * survivalProb` where `subInputDataO = (1 - hazardO)`. This is the cause-specific subdensity: `h_Y(t) * (1-h_D(t)) * S(t-1)`.

**Known bug in Rojas-Saunero `calculateCumInc()`**: The original workshop code initializes the first time point as `hazardP * hazardO` instead of the correct `hazardP * (1 - hazardO)`. This is inconsistent with the formula used at subsequent time points. The Dudukina tidyverse translation (see below) corrects this to `hazardP * (1 - hazardO)` for the non-competing case and `hazardO` for the competing case. The bug has negligible impact in the prostate example because baseline `hazardO` is near zero, but the package MUST use the corrected form.

### 4.3 IPW Weight Construction

**Censoring weights:**
```r
predC <- 1 - predict(plrFitC, newdata = data, type = 'response')  # P(not censored)
cumPredC <- cumprod(predC)  # grouped by subject!
w_cens <- 1 / cumPredC      # unstabilized
w_cens_stab <- cumPredCnum / cumPredC  # stabilized
```

**Competing event weights (for direct effect):**
```r
predO <- 1 - predict(plrFitO, newdata = data, type = 'response')  # P(no competing event)
cumPredO <- cumprod(predO)  # grouped by subject!
w_compev <- 1 / cumPredO
w_direct <- w_cens * w_compev  # combined weight for direct effect
```

**Separable effect weights (2020 paper, cross-arm):**
```r
# Predict D-survival under BOTH arms for ALL subjects
predOa0 <- 1 - predict(plrFitO, newdata = clone_a0, type = 'response')
predOa1 <- 1 - predict(plrFitO, newdata = clone_a1, type = 'response')
cum_pred_O_0 <- cumprod(predOa0)  # grouped by subject
cum_pred_O_1 <- cumprod(predOa1)  # grouped by subject
ipw_d <- cum_pred_O_1 / cum_pred_O_0
ipw_sep_eff <- ipw_d / cum_pred_C
```

**Separable effect weights (2021 paper, arm-stratified):**
```r
# Predict Y-survival under TREAT model and STD model for ALL subjects
pred_AKI_TREAT <- 1 - predict(plr_fit_AKI_TREAT, newdata = data, type = 'response')
pred_AKI_STD   <- 1 - predict(plr_fit_AKI_STD,   newdata = data, type = 'response')
cum_TREAT <- cumprod(pred_AKI_TREAT)  # grouped by subject
cum_STD   <- cumprod(pred_AKI_STD)    # grouped by subject
# Weight for (a_Y=0, a_D=1): reweight Y-hazard from std to treat
ipw_y <- cum_STD / cum_TREAT   # simplified; see 2021 code for t=0 edge case
ipw_sep_eff <- ipw_y / cum_pred_CENS
```

### 4.4 Censoring: G-formula vs. IPW

This is a critical distinction that must be communicated clearly in documentation:

**G-formula**: Censoring is handled by **conditioning**. Include covariates that predict censoring (e.g., age, health status) in the hazard models. No separate censoring model is needed. The g-formula integrates over the full covariate history, which accounts for informative censoring as long as the relevant predictors are in the models.

**IPW**: Censoring requires an **explicit model**: `glm(c_event ~ covariates)`. The IPW estimator operates on observed individuals and needs to reweight them to correct for differential loss to follow-up.

**Key detail for IPW**: The censoring model uses the **observed treatment** `A`, not the overridden `a_Y`. This is because censoring weights correct for the censoring mechanism as it actually operated in the observed data, not in the hypothetical regime.

### 4.5 `group_by(id)` for `cumprod()`

All cumulative product operations MUST be grouped by subject. The reference code uses `aggregate(predC ~ patno, FUN = cumprod)` which handles this correctly. In tidyverse:

```r
data %>% group_by(id) %>% mutate(cum_surv = cumprod(1 - hazard))
```

Without grouping, `cumprod()` bleeds across subjects and produces nonsensical results. This is a common bug source.

### 4.6 NA Handling in Person-Time Data

The event indicators must be set to NA in specific situations to correctly define risk sets:

```r
# Censored individuals: remove from both risk sets at censoring time
longData$y_event[longData$c_event == 1] <- NA
longData$d_event[longData$c_event == 1] <- NA

# Competing event: remove from Y risk set at competing event time
longData$y_event[longData$d_event == 1] <- NA
```

This ensures that `glm()` drops the correct rows when fitting each model (via `na.action`).

### 4.7 Direct Effect: Setting h_D = 0

For the direct (eliminative) effect, the competing event hazard is set to exactly zero:

```r
treated$hazardO <- 0                                    # elimination
treated$s <- (1 - treated$hazardP) * (1 - treated$hazardO)  # simplifies to (1 - hazardP)
```

This is conceptually clean: what would the risk of Y be if D could never occur?

For IPW, the direct effect uses `w_direct = w_cens * w_compev`, which effectively weights the data as if the competing event never happened.

---

## 5. Data Format Conventions

### Input (Person-Level / Wide)

```
id | time | y_event | d_event | c_event | treatment | age | sex | ...
1  | 24   | 1       | 0       | 0       | 1         | 65  | M   | ...
2  | 36   | 0       | 1       | 0       | 0         | 72  | F   | ...
3  | 60   | 0       | 0       | 1       | 1         | 58  | M   | ...
```

Each row is one subject. Exactly one of `y_event`, `d_event`, `c_event` is 1 (or all 0 if administratively censored at end of follow-up). `time` is the event/censoring time.

### Internal (Person-Time / Long)

```
id | dtime | y_event | d_event | c_event | treatment | A_y | A_d | age | ...
1  | 0     | 0       | 0       | 0       | 1         | 1   | 1   | 65  | ...
1  | 1     | 0       | 0       | 0       | 1         | 1   | 1   | 65  | ...
...
1  | 23    | 1       | 0       | 0       | 1         | 1   | 1   | 65  | ...
```

`A_y` and `A_d` are the duplicated treatment columns for the cross-arm trick. They equal `treatment` in the observed data but can diverge in cloned datasets.

---

## 6. Naming Conventions (from Reference Code)

The reference implementations use opaque variable names (`hazardP`, `hazardO`, `Orx`) that obscure the treatment decomposition. The package must use names that make the `A_Y`/`A_D` structure explicit. Mapping:

| Reference (prostate) | Reference (SPRINT) | Dudukina | Package internal |
|----------------------|---------------------|----------|------------------|
| `patno` | `MASKID` | `patno` | `id` |
| `dtime` | `AKI_SAE_MONTH_C` | `dtime` | `time` |
| `rx` | `INTENSIVE` | `rx` | `treatment` / `a_y` |
| `Orx` | (arm-stratified) | `odeath_rx` | `a_d` |
| `prostateDeath` | `AKI_SAE_EVNT` | `prostateDeath` | `y_event` |
| `otherDeath` | `DEATH_BF_AKI_EVENT` | `otherDeath` | `d_event` |
| `eventCens` | `CENSORING_EVENT` | `eventCens` | `c_event` |
| `hazardP` | — | `hazard_pdeath` | `h_y` |
| `hazardO` | — | `hazard_odeath` | `h_d` |
| `s` | — | `s` | `survival` |
| `plrFitP` | `plr_fit_AKI_*` | `pdeath_fit` | `model_y` |
| `plrFitO` | — | `odeath_fit` | `model_d` |
| `plrFitC` | `plr_fit_CENS` | `pcens_fit` | `model_c` |

The `h_y` / `h_d` naming ties each hazard to the event it models (Y or D). The treatment component that *governs* each hazard (`A_Y` or `A_D`) is a property of the model specification and the cloned dataset, not the hazard output itself.

---

## 7. Modeling Choices and Defaults

### Default Hazard Model Formulas

The reference code uses polynomial time with treatment interaction:

```r
# Y-hazard (prostate)
prostateDeath ~ rx * (dtime + I(dtime^2) + I(dtime^3)) + normalAct + ageCat + hx + hgBinary

# D-hazard (prostate)
otherDeath ~ Orx * (dtime + I(dtime^2) + I(dtime^3)) + normalAct + ageCat + hx + hgBinary
```

Default formula pattern for the package:
```r
event ~ treatment * (time + I(time^2) + I(time^3)) + covariates
```

Users should be able to override this with custom formulas (e.g., splines, different interaction structures).

### Arm-Stratified vs. Pooled Models

The 2020 code uses **pooled models** (single model with treatment as covariate). The 2021 code uses **arm-stratified models** (separate models for treated and control). Both are valid:

- Pooled: more efficient, relies on correct specification of treatment-covariate interactions
- Stratified: more flexible, avoids interaction specification, but less efficient with small samples

The package should support both via an argument like `stratified = FALSE` (default pooled).

---

## 8. Testing and Validation

### Primary Validation

The package must reproduce the results from the prostate cancer analysis in the workshop code. Known results to match:
- Total effect RD at 60 months
- Direct effect RD at 60 months  
- Separable effect cumulative incidence curves for arms (1,1), (0,0), (1,0)

### Edge Cases to Test

- No censoring in the data (censoring model should be skipped or trivial)
- No competing events observed (separable effects collapse to total effect)
- Single time point
- Very late censoring (e.g., all censoring after month 50, as in the prostate data)
- Time-varying covariates with LOCF imputation

---

## 9. Documentation and Vignettes

A blog series at alexmourer.com covers this material pedagogically across multiple parts:

- **Part 1-3**: Background on competing events, causal frameworks, estimands
- **Part 4a**: The 2020 paper (basic separable effects)
- **Part 4b**: The 2021 paper (extensions, conditions, G-notation)
- **Part 5**: Implementation (g-formula and IPW code walkthroughs)

The existing HTML output files (`part4b-technical-2021.html`, `part5-implementation.html`) contain detailed explanations that can be adapted into package vignettes. Key pedagogical principles established:

- Show raw person-level data before person-time data
- Show full target formulas in equation boxes before each code section
- Explain g-formula and IPW censoring differences carefully and distinctly
- Use collapsible sections for content repeated across sections

---

## 10. Scope Boundaries

**In scope for v1:**
- Discrete-time pooled logistic estimation
- G-formula and IPW methods
- Total, direct, and separable effects
- Stabilized and unstabilized IPW weights
- Bootstrap confidence intervals
- Prostate cancer example dataset and vignette

**Out of scope for v1:**
- Continuous-time (Cox) estimation → Martinussen & Stensrud (2023) provide the theory; could be a future extension
- Time-varying treatments
- Mediation analysis / illness-death models → future `CausalSurvivalMediation` package (theory of Didelez 2019, extended by Breum et al. 2024)
- Doubly robust / TMLE estimators (one-step estimators from Martinussen & Stensrud 2023 are a potential v2 addition)
- Sensitivity analysis for separability conditions
