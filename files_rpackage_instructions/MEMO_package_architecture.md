# MEMO: Package Architecture Decisions (2026-04-15)

## API

Single entry point: `causal_cr()` with `method = "both"` (default), `"gformula"`, `"ipw1"`, `"ipw2"`.

```r
fit <- causal_cr(
  data,
  id         = "id",
  time       = "time",
  event      = "event_type",
  treatment  = "A",
  covariates = c(...),
  event_y    = NULL,          # NULL = 0/1/2 convention (integer only)
  event_d    = NULL,          # character event → must specify all three
  event_c    = NULL,
  method     = "both",        # "both", "gformula", "ipw1", "ipw2"
  times      = NULL,          # NULL = 12 equally-spaced intervals
  eval_times = NULL,          # where to print contrasts table
  formulas   = NULL,          # override model formulas
  time_varying = NULL,        # not implemented v1, error if used
  extreme_weight_adjust = "truncate",   # "truncate", "trim", "none"
  extreme_weight_threshold = 0.999,
  bootstrap  = FALSE,
  n_boot     = 500,
  alpha      = 0.05
)
```

## Input validation

### Hard errors (stop):
- Column doesn't exist
- treatment not binary
- event not exactly 3 unique values
- character event without event_y/event_d/event_c specified
- method unknown
- duplicate rows
- duplicate (id, time) pairs in person-time data
- time < 0
- time_varying not NULL (v1)

### Warnings (continue):
- Covariates with NAs
- Covariate with > 20 categories
- Few events (< 5)
- time == 0 events → set tstart = -1e-7
- Continuous (non-integer) times → discretize with message

### Data format detection:
- n_distinct(id) == nrow → subject-level, expand internally
- n_distinct(id) < nrow, no duplicate (id,time) → person-time, use as-is

### Event specification:
- Integer 0/1/2 + event_y/d/c all NULL → convention (0=cens, 1=Y, 2=D), print message
- Character → must specify event_y, event_d, event_c explicitly, error otherwise
- Always enforce exactly 3 unique values

## Time handling

- times = NULL → 12 equally-spaced intervals from 0 to max event time
- times = integer scalar → that many equally-spaced intervals
- times = vector → use as cut points
- time == 0 edge case: tstart = -1e-7 with warning
- Continuous data: auto-discretize based on times argument

## Outputs

### Return object (S3 class "causal_cr"):
```r
$cumulative_incidence   # ν(1,1,k), ν(0,0,k), ν(1,0,k) at all k
$contrasts              # RD and RR at eval_times
$models                 # fitted glm objects (model_y, model_d, model_c)
$person_time            # expanded person-time data
$diagnostics            # weight summaries, truncated/trimmed IDs
$bootstrap              # NULL or bootstrap results with CIs
$call, $method, $times, $n
```

### Contrasts:
- Risk differences (additive) AND risk ratios (multiplicative)
- For: total, separable direct (A_Y), separable indirect (A_D)
- At eval_times

### Plots (big part of the package):
- Cumulative incidence curves (3 arms)
- Effect-over-time (differences between curves, computed from $cumulative_incidence)
- TBD: more plot types

## Hazard models

- Default formula: additive A + k + k² + k³ + covariates
- User can override with formulas argument
- G-formula: fit_y_haz + fit_d_haz
- IPW1: fit_d_haz + fit_cens
- IPW2: fit_y_haz + fit_cens
- method = "both": all three models fitted

## Weight handling

- extreme_weight_adjust = "truncate" (default): cap at percentile, keep everyone
- extreme_weight_adjust = "trim": remove subjects beyond threshold
- extreme_weight_adjust = "none": no adjustment
- One or the other, not both
- Flagged IDs stored in $diagnostics

## Sensitivity analysis

| Type | What | How | Version |
|------|------|-----|---------|
| Model SA | Models correctly specified? | G-formula vs IPW1 vs IPW2 comparison | v1 |
| Weight diagnostics | Weights well-behaved? | Inspect W_D, W_Y, W_C distributions | v1 |
| Isolation check | Which isolation level? | W_D, W_Y departure from 1 → full/partial/none | v1 |
| t_k sensitivity | Robust to dismissible violations? | Re-estimate under grid of t_k | v1 |
| Z_k partition check | Proposed covariate partition valid? | Compute W_LD, W_LY with user partition | v2 |
| Isolation adjustment | Time-varying common causes matter? | Compare with/without W_LD, W_LY | v2 |

All SA checks require IPW representations → reason for method = "both" default.

## Positivity / overlap

- G-formula hides positivity problems (extrapolates silently)
- IPW makes them explicit (extreme weights)
- Running both: IPW weights serve as diagnostic for g-formula reliability
- No independent g-formula positivity diagnostic known

## Deferred to later versions

- time_varying covariates (v2): placeholder argument, error if used
- IPCW + g-formula hybrid: censor_weights argument, not implemented v1
- Direct effect (zombie, h_D = 0): v1.1
- Competing event cumulative incidence: already in arm_11/arm_00
- Stabilized weights: v1.1
- Z_k partition check: v2
- Continuous-time (Cox) estimation: future
