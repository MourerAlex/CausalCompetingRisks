# CausalCompetingRisks — TODO / Decisions log

Running log of decisions, caveats, and pending work.
Organized by topic rather than chronologically.

---

## Rep 2 IPW — positivity & zero-hazard safeguards

**Context.** Rep 2 constructs the Y-swap weight:

```
W_Y(s) = (h_Y^{a_y}(s) / h_Y^{a_d}(s)) * prod_{j<s} (1 - h_Y^{a_y}(j)) / (1 - h_Y^{a_d}(j))
```

If any `h_Y^{a}(s) = 0` (literally zero, not just small), the ratio becomes
`0/0 = NaN` or `x/0 = Inf`.

**What hazard = 0 means in practice:**

1. **Positivity violation** — the model says "no subject like this could ever
   have Y at this time under arm a". IPW's core assumption (non-zero prob of
   assignment in every covariate stratum) fails exactly here.
2. **Floating-point underflow** — logistic predictions like `1e-300` round
   to 0.
3. **Extreme covariate profile** — highly protective combination pushing
   the hazard arbitrarily close to 0.

All three are red flags. `check_fitted_positivity()` already warns when
`min(probs) < 1e-6` or `max(probs) > 1 - 1e-6`.

**Current implementation** (ipw.R, `weighted_arm_cum_inc_rep2`):
- Uses `eps <- .Machine$double.eps^0.5` floor (~1.5e-8) on both num and den
  of the hazard ratio and the lag-cumprod ratio.
- Keeps ratios finite; extreme ratios (>1e6) still propagate to `w_d` /
  `w_total` and trigger weight truncation downstream.

**TODO:**
- [ ] Confirm eps choice (1.5e-8) doesn't silently distort estimates when
      models are well-behaved. Add unit test on synthetic data.
- [ ] In diagnostic output, report **rows** where `haz_y_a{0,1} < eps` — these
      are the concrete positivity-violating subject-intervals.
- [x] User-facing documentation: `@section Positivity and zero-hazard
      safeguards` added to `causal_cr()` roxygen. Explains the three failure
      modes and what the user should do (inspect weights, simplify models,
      restrict target population, compare Rep 1 vs Rep 2).
- [ ] Also add similar section to `bootstrap()` and `reweight()` docs when
      those are implemented — mention that bootstrap replicates inherit the
      same positivity issues and reweight() can help explore thresholds.

---

## Rep 1 IPW — contemporaneous cumprod (not lagged)

**Status.** Current code is correct.

W_D = `cumprod_one_minus_hazd_{a_d} / cumprod_one_minus_hazd_{a_y}` — no lag.

**Why no lag?** D appears in the g-formula only as `(1 - h_D)` factors:
- `(1 - h_D(s))` at current interval (no-D-at-current)
- `prod_{j<s}(1 - h_D(j))` at previous intervals

Both are the same type of factor. They combine into a single `prod_{j<=s}(1 - h_D(j))` — what our `cumprod_one_minus_hazd_a*` column holds. No split needed.

Y is different (lag-trick required) because Y appears as both `h_Y(s)` and
`(1 - h_Y(j))` — two different types of ratios.

---

## `fit$cumulative_incidence` as per-method list

**Decision.** `fit$cumulative_incidence` is a named list with entries per
method run:
- `$gformula` — g-formula estimates (NULL if not run)
- `$ipw1` — IPW Rep 1 (NULL if not run)
- `$ipw2` — IPW Rep 2 (NULL if not run)

**TODO:** update downstream accessors/plot/print to iterate over list entries.

---

## `fit$contrasts` — removed

**Decision.** Do not store contrasts in the fit object. Contrasts without
confidence intervals are misleading — they look precise when they're not.

Instead: `contrast(fit, method = "X", ci = boot_obj)` computes on demand.
If no bootstrap is provided, shows point estimates with a warning about
missing uncertainty.

**TODO:** add analytic CI method (influence function / sandwich) as an
alternative to bootstrap for v1.x.

---

## Bootstrap — separate from `causal_cr()`

**Decision (Option B).** Bootstrap returns a standalone object, not
attached to fit:

```r
fit  <- causal_cr(pt)
boot <- bootstrap(fit, n_boot = 500)
plot(risk(fit), ci = boot)
contrast(fit, method = "gformula", ci = boot)
```

**Rationale.** Keeps `fit` pure (point estimates). Bootstrap is optional,
potentially slow, and users might want different settings (n_boot, alpha,
different method subset). Separate object makes this explicit.

**TODO:**
- [ ] Implement `bootstrap()` as exported function taking a `fit` object.
- [ ] Decide if `bootstrap()` re-fits on resamples (current `bootstrap_causal_cr`
      does this) or uses some Monte Carlo pool.
- [ ] Define `"causal_cr_bootstrap"` S3 class with `print()` / `plot()` methods.

---

## `fit$weights` slot (IPW runs only)

**Structure:**
```r
fit$weights = list(
  pt_data_weighted          = <pt_data with haz_*, cumprod_*, w_*_raw, w_*>,
  weight_summary            = <data.frame>,
  truncated_ids             = <vector>,
  extreme_weight_adjust     = "truncate" | "trim" | "none",
  extreme_weight_threshold  = 0.999
)
```

**Purpose:** enables `reweight()` without refitting models.

**TODO:**
- [ ] Implement `reweight(fit, extreme_weight_adjust, extreme_weight_threshold)`
      that reads `fit$weights$pt_data_weighted`, re-applies truncation to raw
      weights, re-runs estimators.

---

## `fit$gformula` slot — reserved

**v1:** NULL.
**v1.x:** per-arm clones or subdensity matrices for diagnostics.
**v2:** Monte Carlo g-formula trajectories.

---

## Rep 2 lagged cumprod — data column

**Column:** `lag_cumprod_one_minus_hazy_a{0,1}`

**Computation:**
```r
ave(
  cumprod_one_minus_hazy_a1,
  id,
  FUN = function(x) c(1, x[-length(x)])
)
```

Per subject: cumprod up to (not including) current interval. First interval
gets lag = 1 (empty product).

---

## Warning handling

**Mechanism.** `causal_cr()` wraps its body with `withCallingHandlers()`
that captures warning messages, muffles them during the run, then re-emits
them as a single group at the end. Storage in `fit$warnings`.

Also, `print(fit)` reports the count: `"Warnings: 3 - use fit$warnings to inspect"`.

**Categories (conceptual):**
- `validation` — covariate quality, rare events, time=0 (from `validate_*`)
- `models` — positivity warning per model (from `check_fitted_positivity()`)
- `weights` — truncation notice (from `ipw_estimate()`)

**TODO:** group warnings in the final `fit$warnings` by category instead of
flat character vector. Currently everything gets concatenated.

---

## Method semantics

- `"gformula"` — only g-formula
- `"ipw1"` — only IPW Rep 1
- `"ipw2"` — only IPW Rep 2
- `"all"` — runs all three (replaces the removed `"both"`)

---

## arm_01 sensitivity arm

**Context.** The main decomposition (A) uses arm (1, 0). Decomposition B
uses arm (0, 1). Both sum to the same total effect but differ whenever
there's an `A_Y x A_D` interaction on the additive scale.

**Status:**
- [x] g-formula returns arm_01 as part of the cumulative_incidence
      data.frame (always, no opt-in).
- [ ] IPW Rep 1 / Rep 2 arm dispatchers still compute only 3 arms.
      Extend to 4 for consistency.
- [ ] Bootstrap: `arm_names` now includes arm_01 but contrast CI array
      only computes Decomposition A (sep_direct / sep_indirect based on
      arm_10). Extend to also compute Decomposition B contrasts
      (direct_B = arm_11 - arm_01, indirect_B = arm_01 - arm_00) when
      arm_01 is populated.
- [ ] `contrast()` accessor: add `decomposition = c("A", "B")` argument
      or return both by default.
- [ ] Plot methods: show 4 lines when all arms present; Decomposition B
      effects in a paired color palette.

**Note.** Decomposition A (reference `a_D = 0` for direct, `a_Y = 1` for
indirect) matches the classical mediation "natural direct at control
mediator" convention. Decomposition B flips the reference. Interpretation
(user intuition from part 5 draft):
- A: "what if the active drug had a modified A_D component whose effect
  on D matched placebo's"
- B: "what if placebo carried an added A_D component whose effect on D
  matched the active drug's"

---

## Left-truncation

**Status: hard-rejected in validate_person_time.** Every subject must
have a row at k = 0. Not supported in v1; **will never be supported** in
future versions — structural assumption of the discrete-time pooled
logistic framework.

---

## Treatment coercion

**TODO.**

User-facing: `treatment` can be character, factor, or numeric. Always coerce
to 0/1 integer internally. If user doesn't specify `treatment_value`,
the "higher" level (factor levels[2], sorted max, etc.) maps to 1.

- [ ] Implement coercion in `to_person_time()` (not `causal_cr()`, so person-time
      input is already coded).
- [ ] Add `treatment_value` argument for explicit mapping.
- [ ] Store original labels as attribute for later display.

---

## Tests (per file)

- [ ] `test-validate.R` — update for split into `validate_subject_level()` and
      `validate_person_time()`
- [ ] `test-data_prep.R` — `to_person_time()` output structure, attributes,
      edge cases (event_time = 0, character events)
- [ ] `test-hazard_models.R` — model selection by method, positivity warning
- [ ] `test-gformula.R` — cumulative incidence computation, per-arm clones
- [ ] `test-ipw-rep1.R` — W_D construction, truncation, weighted hazards
- [ ] `test-ipw-rep2.R` — W_Y lag trick, hazard-zero safeguards, Rep1-vs-Rep2
      agreement on arm (1,1) and (0,0)
- [ ] `test-bootstrap.R` — resampling, percentile CI
- [ ] `test-accessors.R` — risk(), contrast(), diagnostic()
- [ ] `test-plot.R` — smoke tests (no error, returns ggplot)
- [ ] `test-print.R` — print methods for each class

---

## Deferred to v1.1 / later

- [ ] Time-varying covariates (blocked in v1 via `time_varying != NULL` error)
- [ ] Analytic confidence intervals (influence-function / sandwich)
- [ ] Bootstrap parallelization (`future.apply`)
- [ ] Bootstrap progress bar (`progressr`)
- [ ] Third non-intervened cause of death (beyond censoring)
- [ ] Direct effect (zombie / Geneletti) estimand
- [ ] Monte Carlo g-formula (gfoRmula-ICE style)
- [ ] Isolation hierarchy check (full / partial-D / partial-Y / no-iso+Z_k)
- [ ] `t_k` sensitivity analysis
- [ ] Risk table below plots (at-risk / events / censored by arm,
      adjustedCurves-style)
- [ ] Diagnostic plots: weight distributions, W_D/W_Y departure from 1
- [ ] Hybrid censoring weights with strata

---

## Accessor updates for new fit structure

**Context.** `fit$cumulative_incidence` is now a named list (`$gformula`,
`$ipw1`, `$ipw2`). Accessors need to adapt.

**TODO:**
- [ ] `risk(fit, method = NULL)` — default shows all methods stacked;
      `method = "gformula"` returns single-method view.
- [ ] `contrast(fit, method = "gformula", ci = NULL)` — requires method arg
      (no hidden default). If `ci` (a `causal_cr_bootstrap` object) passed,
      include CI columns. If not, warn about missing uncertainty.
- [ ] `diagnostic(fit)` — unchanged concept but reads from `fit$weights`
      (weight distributions, truncated IDs, positivity flags).

---

## Plot methods — new signatures

**Purpose.** Accept an optional `ci` argument (`causal_cr_bootstrap` object)
to show confidence bands.

**Color palette (Okabe-Ito, colorblind-safe):**
- Total effect / arm (1,1) line: black `#000000`
- Separable direct (A_Y) / arm (1,0): green `#009E73`
- Separable indirect (A_D): vermillion `#D55E00`
- Control / arm (0,0): blue-ish `#0072B2`
- Fourth arm (0,1) if added: yellow `#F0E442` or sky blue `#56B4E9`

**TODO:**
- [ ] `plot(risk(fit), ci = boot, method = "gformula")` — cumulative
      incidence lines + optional CI ribbons.
- [ ] `plot(contrast(fit, method, ci = boot), type = "rd")` — effect curves.
      Total black, sep direct green, sep indirect vermillion, dashed zero
      reference.
- [ ] `plot(diagnostic(fit))` — weight distributions (histograms), W_D / W_Y
      departure from 1 over time.
- [ ] Risk table below plots (adjustedCurves-style): rows = arms, columns
      = time points, cells = at-risk / events / censored counts.
      Unweighted counts in v1 (discussed).
- [ ] Allow customization: title, subtitle, legend position. (Confirmed:
      embedded in ggplot, users can `+ labs(...)` — no explicit args needed.)

---

## Print / summary / confint methods

**TODO:**
- [ ] `print.causal_cr()` — list cumulative incidence per method at final
      time + warning count. Already does method + n + time points. Needs
      update for list-valued `fit$cumulative_incidence`.
- [ ] `summary.causal_cr()` — table of point estimates per method; note
      that contrasts / CIs require `contrast(fit, ci = bootstrap(fit))`.
- [ ] Remove `confint.causal_cr()` — CIs belong to the bootstrap object now,
      use `confint(boot)` instead.
- [ ] Add `print.causal_cr_bootstrap()` — n_boot, alpha, per-method arm
      summaries at eval_times.

---

## Sensitivity analyses

Two distinct kinds discussed:

### Estimand sensitivity
- [ ] Arm (0, 1) — symmetric alternative to (1, 0). Disagreement → assumption
      failure or poor overlap.
- [ ] Compare Rep 1 vs Rep 2 estimates (already built-in via `method = "all"`).
      Disagreement → model misspecification.

### Model sensitivity
- [ ] Compare g-formula vs IPW under same covariate set.
- [ ] User-supplied `formulas` list to override defaults for Y / D / C.
- [ ] Under-fit vs over-fit stress: report estimates under minimal model vs
      full covariate adjustment.

### Isolation hierarchy check

From Stensrud 2021 Section 5 / Lemma 4:

- [ ] Full isolation — `w_d` and `w_y` close to 1 everywhere → strongest
      assumption, holds.
- [ ] Partial D-isolation — `w_d ≈ 1`, `w_y ≠ 1`.
- [ ] Partial Y-isolation — `w_y ≈ 1`, `w_d ≠ 1`.
- [ ] No isolation + `Z_k` partition — neither weight near 1, but a time-varying
      partition variable makes the decomposition identifiable.
- [ ] No isolation, no partition — not identifiable.

Implementation:
- [ ] Report summary statistics of `w_d`, `w_y` (mean, quantiles, distance
      from 1).
- [ ] Given these, suggest which isolation regime is plausible.
- [ ] Add `Z_k` argument for user to specify partition variable(s).

### t_k sensitivity

Stensrud 2021 Section 7:

- [ ] For user-supplied ranges of an unmeasured-confounding parameter `t_k`,
      compute bounds on the treatment effect.
- [ ] Output: lower/upper bound per time, per contrast.
- [ ] Plot: sensitivity curve with shaded region.

---

## Post-fit GLM checks

**Context.** Separate from the input validation — these are run AFTER
model fitting.

- [x] `check_fitted_positivity()` — warn when min/max predicted probability
      within `eps` of 0 or 1. (Implemented.)
- [ ] Check for collinearity (glm drops a term → `NA` coefficient).
- [ ] Check for convergence (`!model$converged`).
- [ ] Check for extreme standard errors (indicator of separation).
- [ ] Report these consolidated as `fit$model_diagnostics`.

---

## Review queue — skipped blocks

We need to revisit these written-without-walkthrough blocks:

1. `check_fitted_positivity()` in `hazard_models.R`
2. `w_cens_raw` / `w_d_raw` preservation logic in `ipw.R`
3. Truncation warning wording in `ipw.R`
4. `check_covariate_quality()` consolidation in `validate.R`
5. `withCallingHandlers` wrapper + warning capture in `causal_cr.R`
6. Warning count line in `print.causal_cr()`
7. `cumprod_one_minus_*` rename (final naming decision)
8. `estimate_weighted_cum_inc_rep1` dispatcher
9. `weighted_arm_cum_inc_rep1` full body
10. `cum_inc_from_weighted()` shared helper
11. `weighted_arm_cum_inc_rep2` full body (W_Y lag trick)

---

## End-to-end validation

- [ ] Run full pipeline on prostate data. Expected (from part5_run.R, n=252
      subset):
      - g-formula: arm_11 ≈ 0.23, arm_00 ≈ 0.33, arm_10 ≈ 0.25
      - IPW Rep 1: similar numbers
      - IPW Rep 2: similar numbers (cross-check with Rep 1)
- [ ] Snapshot test: save expected values to a test fixture so we detect
      regressions.

---

## Bundled data

- [x] `prostate_data.rda` in `data/` (prepared from `Hmisc::prostate`, subset
      to placebo + 5.0mg DES, labels stripped).
- [ ] Add documentation page (R/data.R already has it — review content).
- [ ] Internal `simulate_causal_cr_data()` helper for controlled unit tests
      (known true effects, edge cases). Not exported.

---

## Future: causal subgroup analysis (novel application)

**Idea.** The same separable-effects machinery (A_Y / A_D decomposition,
swap weights, per-arm Hajek CIF) can be adapted to **causal subgroup
analysis**: identify causal effects within subgroups by treating subgroup
membership as a decomposable component analogous to A_D.

**Status.** Methodologically not yet formalized in the literature. Likely
straightforward to implement given the v1 plumbing is in place.

**Why it matters.** Opens a novel methodological space using already-built
machinery. Potential publication / coauthor opportunity (Stensrud or others)
once v1 of CausalCompetingRisks ships and the causal-mediation companion
is established.

**Form.** Companion package or vignette demonstrating the subgroup-effect
decomposition. Reuses `ipw_longitudinal()`, `swap_d_weights()`,
`swap_y_weights()`, and the per-arm CIF dispatcher unchanged.

**Sequencing.** Post-v1 only. Do not block v1 release on this.

---

## Future: shared `causal_tools` / `causal_utils` package

**Context.** Once CausalCompetingRisks stabilizes and we start building a
companion causal-mediation package (and potentially others), some building
blocks will be shared across packages.

**Candidates to factor out** (in a separate low-level utilities package):

- `to_person_time()` / survSplit wrapper + temporal-ordering logic
- `validate_subject_level()` / `validate_person_time()` / `check_covariate_quality()`
- `fit_hazard_models()` / `check_fitted_positivity()` / `fit_one_model()`
  (pooled-logistic hazard fitting with model-check diagnostics)
- `weighted_hazard_by_k()` (Hajek-ratio weighted hazard)
- `cum_inc_from_weighted()` shape (discrete-time cumulative incidence)
- `%||%` and other base utilities

**What stays in CausalCompetingRisks:**

- The separable-effects specific math (`weighted_arm_cum_inc_rep1/rep2`,
  arm definitions, cross-arm prediction trick)
- The `causal_cr()` entry point wiring
- Plot / print methods specific to the separable-effects output

**What would stay in a mediation-specific package:**

- Mediation-specific estimators (natural direct / indirect effects,
  randomized interventional analogs)
- Mediation-specific sensitivity analyses

**When to do this.** Not for v1 of either package. Factor out on the
second time we need something — avoid premature abstraction. Once both
packages are stable enough that we can see the repeated patterns, make
a `causal_tools` (or similar) package and have both depend on it.

---

## Deferred to later version

- [ ] **`reweight(fit, extreme_weight_adjust, extreme_weight_threshold)`**
      (v1.x). Sensitivity-analysis helper: re-applies weight truncation /
      trimming to an existing IPW fit without re-running glm or prediction.
      Uses the stored raw weights (`w_cens_raw`, `w_d_arm_10_raw`,
      `w_d_arm_01_raw`) on `fit$weights$pt_data_weighted`. Fast — fit once,
      reweight over a grid of thresholds. Not blocking for v1 since all the
      raw data to enable it is already preserved.

---

## R internals to learn / teach

- [ ] **Condition system deep dive**: `withCallingHandlers` vs `tryCatch`,
      calling-handler scope and why `warning()` inside a handler recurses,
      `invokeRestart("muffleWarning")` semantics, the capture-and-replay
      pattern for labeled re-emission (see `predict_with_warning()` in
      `R/gformula.R`). Either I teach this, or we find / write an
      explainer — should be comfortable with it before extending the
      warning collection logic further.

---

## Build / release

- [ ] `roxygen2::roxygenise()` to populate `man/` and `NAMESPACE`
- [ ] `R CMD check` passing (no NOTEs, WARNINGs, ERRORs)
- [ ] Vignette(s) — at least one walkthrough on prostate data, one on
      sensitivity analyses.
- [ ] **Interpretation vignette** — explicit guidance on how to phrase
      conclusions from the numeric output. The target user is an
      epidemiologist unfamiliar with the separable-effects framework;
      they need ready-to-paste sentence templates. For each numeric
      quantity the package returns, show:
      - What the number literally means (counterfactual definition)
      - What it does NOT mean (common misreading)
      - A **"In your analysis, the correct phrasing would be:"** block
        with a worked example the user can copy and adapt
      - Examples to cover:
        * Total RD / RR
        * Separable direct effect (Decomposition A vs B) — emphasize the
          `a_D` reference value choice
        * Separable indirect effect (A vs B) — same
        * Rep 1 vs Rep 2 agreement as a sensitivity check
        * Arm-level cumulative incidence (not a contrast)
        * Confidence intervals from the bootstrap — percentile CIs, what
          the 95% interval means, when they're unreliable (positivity
          flagged)
      - Flag common mistakes to avoid in the phrasing: "treatment
        causes X" (should be "under intervention ...") / "the direct
        effect is Y" (should specify the reference arm) / overclaiming
        from wide CIs.
- [ ] `pkgdown` site (optional for v1)
- [ ] README.md with install + quick example
- [ ] NEWS.md to track version changes
- [ ] LICENSE.md (MIT template, already in DESCRIPTION)
- [ ] Consider CRAN / r-universe submission after v1 stabilizes

---

## Future structural extensions (v3+, design sketch)

These are model-structure extensions that go beyond v1's (Y, D, C) single-cause
competing-events setup. Kept here so we remember why the current architecture
was chosen and which seams to keep clean.

### Two censoring nodes (independent + dependent)

**Current.** One `c_flag` node with a single `censoring_weights = TRUE/FALSE`
switch.

**Future.** Split censoring into two mechanisms:
- `c_indep_flag`: administrative / truly independent (study end). Skip model;
  `w_cens_indep = 1`.
- `c_dep_flag`: informative (loss-to-follow-up dependent on covariates /
  treatment). Fit glm as today; `w_cens_dep = 1 / cumprod(1 - haz_c_dep)`.

**Combined weight.** `w_cens = w_cens_indep * w_cens_dep = w_cens_dep` (since
the independent piece is 1). User controls via a per-node flag.

**API sketch.** `censoring = list(indep = "c_admin_flag", dep = "c_ltfu_flag")`
or two separate columns with `censoring_weights = c("indep" = FALSE, "dep" = TRUE)`.

**Motivation.** Administrative censoring genuinely independent of outcomes; no
model needed and fitting one just adds noise. Informative censoring requires a
model. Current single-flag design conflates these.

### Two treatment nodes (per-protocol estimation)

**Caveat on theoretical grounding.** Stensrud et al. 2020/2021 develop the
separable-effects framework for **baseline (point) treatment**. Extending
the decomposition to time-varying treatment requires additional
identification conditions and is an active research problem, not a settled
framework. Any implementation here would need careful alignment with the
methodological literature at the time of release. The current v1 restriction
to baseline treatment is principled, matching the original theory.

**Motivation.** "Per-protocol" strategies need to distinguish
initial treatment assignment from subsequent adherence. A single binary `A`
cannot encode "assigned to treatment but dropped out at month 3".

**Future.** Two nodes:
- `A_init`: assigned treatment at baseline (static)
- `A_k`: observed treatment at interval k (possibly time-varying)

**Estimand.** Per-protocol effect = effect of sustained `A_init = a` over all
k (compare `A_k = a` throughout vs. `A_k = 1-a` throughout).

**Implementation.** Needs time-varying covariate support (see v1.1). Person-time
data already structured for this; just need `A_k` treated as a per-interval
variable in the hazard models and weight construction. The separable
decomposition would extend: `A_init_Y`, `A_init_D`, plus per-interval
adherence handling.

**Connection to existing framework.** Overlaps with the `gfoRmula` package's
handling of treatment policies. Could share the "intervention grammar" with
that package once we have both.

### Additional death event (multi-cause, possibly interventionable)

**Motivation.** Extend (Y, D, C) to (Y, D1, D2, C) where D1, D2 are two
distinct competing events. Examples:
- D1 = cardiovascular death, D2 = cancer death (both dependent on treatment)
- D1 = suicide (potentially independent of treatment), D2 = other-cause
  death (treatment-dependent)

**Two flavors:**

1. **Non-intervenable D2**: like censoring but from a biological mechanism.
   Weight construction extends similarly to dependent censoring (separate
   model, enters as a factor in `cumprod_one_minus_haz...`).

2. **Intervenable D2**: same role as D1 (competing event we could imagine
   eliminating or modifying via separable treatment components). Separable
   effects decomposition extends: `A_Y, A_D1, A_D2` — 2^3 = 8 arms instead
   of 4.

**Bounds-style analysis.** When D2 is intervenable but we cannot
plausibly specify `a_D2`, can still produce bounds on the Y-effect by
ranging `a_D2` over {0, 1} and reporting the widest / tightest
decomposition.

**Architecture implications.** The current code assumes a fixed "3 event
categories" structure. For v3 we'd need to generalize:
- `event_labels` becomes a richer object: `list(y = ..., d = list(D1 = ..., D2 = ...), c = ...)`.
- Flag columns generalize: `y_flag`, `d1_flag`, `d2_flag`, `c_flag`.
- G-formula and IPW estimators loop over all intervenable causes for the
  cross-arm construction.
- Arms indexed by tuple `(a_Y, a_D1, a_D2, ...)` — return data frames become
  wider or list-of-arms.

**Connection to existing code.** The generalization is intrusive but
follows the same pattern (per-arm clones + predictions for g-formula;
per-arm swap weights for IPW). Planning ahead: keep the arm-dispatch
pattern (`lapply(arms, ...)`) rather than hardcoding arm_11/arm_00/arm_10
in multiple places so future generalization is local.

### Summary of future extensions

| Extension | Severity of refactor | v target |
|---|---|---|
| Two censoring nodes | Medium (adds one weight factor) | v1.2 |
| Two treatment nodes (per-protocol) | Large (needs time-varying first) | v2 |
| Multi-cause (D1, D2, ...) | Large (arm space generalization) | v3 |
| Bounds on non-fixable intervenable causes | Moderate (builds on multi-cause) | v3 |

**Design principle to preserve now.** Keep arm-indexing logic abstracted
into the dispatcher layer (`estimate_weighted_cum_inc_rep*()`,
`gformula_estimate()`'s `arms` list). Avoid hardcoding arm names
downstream — anywhere we loop over `c("arm_11", "arm_00", "arm_10", "arm_01")`
today should eventually take the arm list as an argument.
