# API v1 — LOCKED

**Source of truth for `CausalCompetingRisks` v1 API.**
Produced in session 2026-04-24 via grill-claude-code + grill-me on separable effects scope.
Superseded decisions live in git history; do not re-debate without updating this file.

---

## Package

- **Package name**: `CausalCompetingRisks` (kept — discoverability on search engines).
- **Scope v1**: Stensrud et al. (2020) separable effects under **full isolation** with **baseline covariates** only.
- Deferred to v2: Stensrud et al. (2021) generalization (`Z_k`, partial isolation, MC g-computation), time-varying covariates (covariate models), composite endpoint, controlling direct effect, marginal net risk.

---

## Primary entry point

```r
separable_effects(data,
                  id         = "id",
                  time_name  = "k",
                  treatment  = "A",
                  event_y    = "Y",
                  event_d    = "D",
                  covariates = c("L1", "L2"),
                  method     = c("gformula", "ipw"),   # default "gformula"
                  ipcw       = FALSE,
                  stabilize  = TRUE,
                  formulas   = NULL,
                  K          = NULL)
```

### Arg semantics

| Arg          | Type                       | Default             | Notes                                                                                   |
|--------------|----------------------------|---------------------|-----------------------------------------------------------------------------------------|
| `data`       | long-format data.frame     | —                   | One row per (id, time). gfoRmula-style. No custom S3 class, no dual-path subject-level. |
| `id`         | character(1) column name   | `"id"`              | Subject identifier.                                                                     |
| `time_name`  | character(1) column name   | `"k"`               | Discrete time index; convention starts at 0. Warn if `min(time) != 0`.                  |
| `treatment`  | character(1) column name   | `"A"`               | Binary {0, 1}. Baseline in v1 (constant per subject).                                   |
| `event_y`    | character(1) column name   | `"Y"`               | Binary {0, 1} — primary event at this row.                                              |
| `event_d`    | character(1) column name   | `"D"`               | Binary {0, 1} — competing event at this row.                                            |
| `covariates` | character vector of cols   | `c("L1","L2")`      | Baseline L (time-fixed in v1). Time-varying → v2.                                       |
| `method`     | `c("gformula", "ipw")`     | `"gformula"`        | Pick one; IPW runs both Rep 1 and Rep 2 for sensitivity.                                |
| `ipcw`       | logical                    | `FALSE`             | See method × ipcw table below.                                                          |
| `stabilize`  | logical                    | `TRUE`              | IPW only. Stabilized weights `P(A) / π(A|L)`.                                           |
| `formulas`   | NULL or named list         | `NULL`              | Partial override allowed. Keys ⊆ `{Y, D, C, A_model}`. RHS-only (e.g. `Y = ~ A + L1`).  |
| `K`          | integer or NULL            | `NULL`              | Horizon for CIF. NULL → `max(data[[time_name]])`.                                       |

### No `...`. All args explicit.

### No censoring column arg
Censoring is inferred from the long-format structure: subject whose final row has `Y=0, D=0` and `time < K` is censored at that row. Standard gfoRmula / gfoRmulaICE convention.

### Method × ipcw table

| `method`    | `ipcw`  | Behavior                                                                                         |
|-------------|---------|--------------------------------------------------------------------------------------------------|
| `gformula`  | `FALSE` | Plain parametric g-formula. Censoring handled via conditioning on time-varying L in Y, D hazards. |
| `gformula`  | `TRUE`  | DR-for-censoring (Bang & Robins 2005; Kawahara et al. 2020). IPCW-weighted g-formula.            |
| `ipw`       | `TRUE`  | Standard IPW + IPCW. Both Rep 1 and Rep 2 computed.                                              |
| `ipw`       | `FALSE` | Warn loudly: assumes independent censoring. Proceeds.                                            |

### IPW defaults (fixed, not user-tunable)
- Weights truncated at percentiles **[0.01, 0.99]** both tails (Cole & Hernán 2008, AJE 168(6):656–664).
- Always computes **both Rep 1 and Rep 2** — consistency comparison as built-in sensitivity.
- Propensity `π(A|L)` fitted once on full data (no per-timepoint models; baseline A).

### Formula escape hatch

**Default** (when `formulas = NULL`):
```r
Y       ~ A + L1 + L2 + k       # primary event hazard
D       ~ A + L1 + L2 + k       # competing event hazard
C       ~ A + L1 + L2 + k       # censoring hazard (only if ipcw = TRUE)
A_model ~     L1 + L2           # propensity, baseline covariates only (no k)
```

**Override convention**:
- Flat named list with keys ⊆ `{Y, D, C, A_model}`. Partial override allowed; unlisted keys → default.
- **RHS only**: `list(Y = ~ A + L1 + splines::ns(k, 3))`. LHS inferred from key.
- **No auto-injection** when user supplies a formula — what they write is what gets fit.

**Warning**: if user supplies `formulas$Y` and `formulas$C` with differing RHS under `method = "gformula", ipcw = FALSE`, emit:
> "Y and C formulas differ. Conditioning may not suffice for dependent censoring. Consider `ipcw = TRUE` to add IPCW (DR-for-censoring, Bang & Robins 2005)."

---

## Return object

**Class**: `"causal_cr_fit"` (S3 list).

**Slots** (names indicative):
```
$call           # matched call
$method         # "gformula" or "ipw"
$ipcw           # logical
$K              # horizon
$data_summary   # N subjects, rows, event counts
$models         # fitted glm objects: $Y, $D, $C, $A_model (NULL if not fit)
$risk           # LONG data.frame: k, arm, cum_inc, method (see below)
$diagnostics    # positivity min(π), max(weight), truncation counts
$bootstrap      # NULL until bootstrap() attaches B samples of $risk
```

### Return shape of `$risk` (and `risk()` accessor)

**Long format, always**:

| k | arm  | cum_inc | method          |
|---|------|---------|-----------------|
| 0 | 11   | 0.000   | gformula        |
| 1 | 11   | 0.021   | gformula        |
| … | 00   | 0.018   | gformula        |
| … | 10   | …       | gformula        |
| … | 01   | …       | gformula        |

Under `method = "ipw"`, additional rows with `method ∈ {"ipw_rep1", "ipw_rep2"}`.

Arms are strings: `"11"`, `"00"`, `"10"`, `"01"` (for `(a_Y, a_D)`).

---

## Uncertainty

```r
fit_b <- bootstrap(fit, B = 500, seed = NULL)
```

- **Cluster-bootstrap by `id`** (always; no row-level option).
- Refits **everything** per replicate: hazards (`Y`, `D`, optionally `C`), propensity (`A_model` if IPW), weights, truncation, CIF computation.
- Returns the same `causal_cr_fit` object with `$bootstrap` slot populated: B × K × |arms| array + percentile summaries.
- Default B = 500.

**TODO** (deferred):
- Rare-events special case (parametric bootstrap fallback / stratified resampling). Note in roxygen as known limitation.

---

## Accessors (all no-arg, no filtering)

| Accessor          | Returns                                                                      |
|-------------------|------------------------------------------------------------------------------|
| `print(fit)`      | Minimal: call, method, N, "use summary/risk/contrast/plot/assumptions".      |
| `summary(fit)`    | Full breakdown; **contrasts shown only if `$bootstrap` populated**.          |
| `risk(fit)`       | Full long df of cumulative incidence — users filter with dplyr.              |
| `contrast(fit)`   | Long df of all predefined contrasts (Decomp A & B, all k). CIs iff bootstrap.|
| `plot(fit)`       | CIF curves by arm; ribbons iff bootstrap.                                    |
| `assumptions(fit)`| Prints identifying assumptions block for the fit's estimand.                 |

### Pre-defined contrasts (returned by `contrast()`)

Both decompositions share the same reference arm (`arm_00`) and the same total
(`arm_11 − arm_00`). They differ only in the **intermediate arm** — the path
order from `arm_00` to `arm_11`.

**Total** (identical under A and B): `arm_11 − arm_00`

**Decomposition A** — vary `A_Y` first, then `A_D` (intermediate arm = `arm_10`):
- **Direct (SDE-A)** = effect of `A_Y` at `a_D = 0` = `arm_10 − arm_00`
- **Indirect (SIE-A)** = effect of `A_D` at `a_Y = 1` = `arm_11 − arm_10`

**Decomposition B** — vary `A_D` first, then `A_Y` (intermediate arm = `arm_01`):
- **Indirect (SIE-B)** = effect of `A_D` at `a_Y = 0` = `arm_01 − arm_00`
- **Direct (SDE-B)** = effect of `A_Y` at `a_D = 1` = `arm_11 − arm_01`

Algebraic identity: `SDE + SIE = Total` under each decomposition.

Convention: arm tuple `"xy"` means `(a_Y = x, a_D = y)`. `A_Y` acts on the path
A → Y not through D (its effect is the *direct* effect, SDE). `A_D` acts on
the path A → D → Y (its effect is the *indirect* effect, SIE). See Stensrud
et al. (2020) §3.

All produced as point estimates; CIs populated on bootstrap.

### `summary()` output contracts
- Pre-bootstrap: arms shown with point estimates; contrasts section says *"Not computed. Effects require bootstrap for valid inference. Run: fit_b <- bootstrap(fit, B = 500)"*.
- Post-bootstrap: contrasts shown with percentile 95% CIs.
- Assumptions never in `print()`; always in `summary()` and `assumptions()`.

### `print()` convention
- No method-flavor parenthetical (`"(conditioning; no IPCW)"` dropped — cleaner).
- Single line for method: `Method: gformula` or `Method: ipw`.

---

## Validation (strict, all three groups active)

### Mandatory errors
- `data` is data.frame-like.
- Required columns exist.
- `treatment`, `event_y`, `event_d` ∈ {0, 1}.
- Monotone increasing `time_name` within each subject.
- At most one event per subject (no rows after Y=1 or D=1).
- `method` ∈ allowed set.
- `formulas` keys ⊆ `{Y, D, C, A_model}`.

### Conventional warnings
- `min(time_name) != 0` — warn.
- `treatment` varies within subject — warn (v1 assumes baseline).
- Covariates missing on event rows — warn.

### Hygiene notes (printed in summary, not warned)
- Count of subjects with covariate missingness.
- Positivity: min observed `π(A|L)`.
- IPW weight max and count truncated.

---

## Implementation fixed defaults

- **GLM fitting**: `glm.control(maxit = 100)` applied to every hazard and propensity fit. Note this in roxygen for each model-fit wrapper.
- **Weight truncation**: [0.01, 0.99] both tails, always on (no user knob).
- **Stabilization**: default `TRUE`, user-toggleable.
- **Horizon**: defaults to `max(data[[time_name]])`.

---

## Citations (required in roxygen)

### Primary methodology
- Stensrud MJ, Young JG, Didelez V, Robins JM, Hernán MA (2020). *Separable effects for causal inference in the presence of competing events.* JASA. doi:10.1080/01621459.2020.1765783
- Stensrud MJ, Hernán MA, Tchetgen Tchetgen EJ, Robins JM, Didelez V, Young JG (2021). *A generalized theory of separable effects in competing event settings.* Lifetime Data Anal. doi:10.1007/s10985-021-09530-8

### IPW / weights
- Cole SR, Hernán MA (2008). *Constructing inverse probability weights for marginal structural models.* AJE 168(6):656–664.

### DR for censoring
- Bang H, Robins JM (2005). *Doubly robust estimation in missing data and causal inference models.* Biometrics 61(4):962–973.
- Kawahara T, Shinozaki T, Matsuyama Y (2020). *Doubly robust estimator of risk in the presence of censoring dependent on time-varying covariates.* BMC Med Res Methodol 20:204.

### Related packages for benchmarking
- `gfoRmula` (Logan, Young, Hernán) — parametric g-formula via MC.
- `gfoRmulaICE` — iterative conditional expectation g-formula.
- `ltmle` (Petersen, van der Laan) — wide-format TMLE.

---

## Implementation TODOs (post-lock)

- [ ] Rename `causal_cr()` → `separable_effects()` across codebase.
- [ ] Remove custom `pt_data` class path. Accept long data.frame directly with column-name args.
- [ ] `glm.control(maxit = 100)` wired into every GLM fit; note in roxygen.
- [ ] Long-format output for `$risk` everywhere. Pivot internal wide → long.
- [ ] Warning: Y/C formulas differ under `ipcw = FALSE`.
- [ ] Warning: `method = "ipw"` + `ipcw = FALSE` (independent censoring assumption).
- [ ] Validation pass matching the strict checklist above.
- [ ] `summary()` pre/post-bootstrap behavior.
- [ ] `assumptions(fit)` accessor with text block (see draft in session).
- [ ] Benchmark harness vs `gfoRmula` on shared DGP.
- [ ] Rare-events bootstrap special case (research item).

## Deferred to v2

- Time-varying covariates (covariate evolution models).
- `Z_k` partition variable + partial isolation (Stensrud 2021).
- MC g-computation path.
- Composite endpoint estimand.
- Controlling direct effect estimand.
- Marginal net-risk estimand.
- User-supplied pre-fitted model objects (coxph / glm) via `formulas`.
- ICE-gformula (v2+, requires internal wide pivot).
