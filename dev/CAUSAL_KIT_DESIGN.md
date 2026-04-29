# CausalKit Design Document (DRAFT)

Plan for the future shared-foundation package in the `RCausalUniverse`
ecosystem. Not implemented yet — this document captures intent so that
when `CausalSeparableMediation` (or any second consumer) starts, the
extraction is guided by a concrete target.

Status: **DRAFT**. Nothing is built. `CausalCompetingRisks` continues
as a single self-contained package for now.

---

## 1. Ecosystem Overview

```
RCausalUniverse (ecosystem name, not a package)
│
├── CausalKit                   ← shared foundation (this doc)
│     Imports: ggplot2, stats, survival
│     Exports: hazards, IPW, g-formula, cuminc, plots, bootstrap, validate
│
├── CausalCompetingRisks          ← separable effects for competing events
│     Imports: CausalKit
│     Exports: causal_cr(), to_person_time() (thin wrappers), contrast(), ...
│
├── CausalSeparableMediation     ← mediation analysis (Didelez et al.)
│     Imports: CausalKit
│
├── CausalTrialEmulation         ← trial emulation design + analysis
│     Imports: CausalKit
│
└── future packages...
```

Strict layering: downstream packages depend on `CausalKit` only. No
cross-dependencies between downstream packages.

---

## 2. Naming — DECIDED: `CausalKit`

Selected name: **`CausalKit`**.

Rationale:
- Dual role: standalone usable (IPW + g-formula + plots can do complete
  analyses on their own) AND foundation for domain-specific packages.
- "Kit" reads naturally in both senses — users pick tools from it
  directly, and downstream packages compose from it.
- Avoids "Core" (implies internals), "Scikit" (Python baggage awkward
  in R), and "Foundation" (too formal).

Rejected: CausalCore, CausalFoundation, causalbase, CausalLib,
CausalEngine, CausalScikit, dtcausal.

---

## 3. Scope Principles

**Belongs in CausalKit** if:
- Used by ≥2 downstream packages (or clearly will be)
- Is framework-generic (not specific to separable-effects decomposition)
- Encapsulates a single, composable primitive
- Has a clean, stable API

**Stays in downstream package** if:
- Encodes a specific estimand (e.g., separable direct/indirect effect)
- Wraps/orchestrates generic primitives into a domain-specific workflow
- Is the user-facing entry point for a particular framework

---

## 4. Functional Inventory

Current `CausalCompetingRisks` functions, labeled by proposed destination.

### 4.1 — Data prep (→ `CausalKit`)

| Function                   | Role                                         |
|----------------------------|----------------------------------------------|
| `to_person_time()`         | subject-level → person-time (survSplit)      |
| `validate_input_shape()`   | NULL / data.frame / nrow checks              |
| `validate_subject_level()` | subject-level structure validation           |
| `validate_person_time()`   | person-time structure validation             |
| `check_covariate_quality()`| covariate quality warnings                   |

All generic. `to_person_time()` would move as-is; downstream packages
use the same conventions (y_flag, d_flag, c_flag, A_y, A_d, k).

### 4.2 — Hazard model fitting (→ `CausalKit`)

| Function                   | Role                                         |
|----------------------------|----------------------------------------------|
| `fit_hazard_models()`      | orchestrator (Y, D, C models)                |
| `fit_one_model()`          | single glm fit + diagnostics                 |
| `check_fitted_positivity()`| per-model diagnostics                        |

Generic: Y/D/C naming is conventional, but the actual mechanism (pooled
logistic regression with polynomial time) works for any discrete-time
hazard regardless of framework.

### 4.3 — G-formula primitives (→ `CausalKit`)

| Function                   | Role                                         |
|----------------------------|----------------------------------------------|
| `make_clone()`             | baseline × cut_times expansion per arm       |
| `predict_hazards()`        | hazards on cloned data                       |
| `predict_with_warning()`   | helper: predict + warning capture            |
| `compute_cum_inc()`        | per-subject cumprod → average cumulative inc |

### 4.4 — IPW primitives (→ `CausalKit`)

| Function                   | Role                                         |
|----------------------------|----------------------------------------------|
| `add_d_hazards()`          | predict D-hazards under both arms + cumprod |
| `add_y_hazards()`          | predict Y-hazards + cumprod + LAGGED cumprod |
| `add_censoring_weights()`  | w_cens, w_cens_raw                           |
| `apply_weight_truncation()`| truncate / trim weights                      |
| `weighted_hazard_by_k()`   | Hajek ratio                                  |
| `cum_inc_from_weighted()`  | weighted hazards → cumulative incidence      |
| `summarize_weights()`      | per-weight distribution stats                |

All generic. `add_d_swap_weights()` stays framework-specific (it
encodes separable-effects swap semantics) — see 4.8.

### 4.5 — Bootstrap (→ `CausalKit`)

The **infrastructure** is generic:
- subject-level resampling with unique synthetic IDs
- `pt_by_id` pre-split for O(1) lookup
- per-replicate `fit_causal_cr()` call
- progress messaging

But `bootstrap()` currently bakes in the separable-effects output shape
(methods × arms × times). Two design options:

**(a)** Keep a generic `resample_subjects()` in `CausalKit`, and keep
the full `bootstrap()` orchestration in `CausalCompetingRisks`.

**(b)** Design a more general `bootstrap()` in `CausalKit` that takes a
`fit_fn` argument — each downstream package passes its own estimator.
Return shape is parameterized.

**(a)** is simpler and what I'd propose first.

### 4.6 — Plotting (→ `CausalKit`)

| Component                         | Role                                |
|-----------------------------------|-------------------------------------|
| `build_risk_table_plot()`         | risk table ggplot primitive         |
| `risk_table_internal()`           | counts by arm at cut times          |
| Step-ribbon transform (helper)    | match step curves                   |
| Okabe-Ito palette (constant)      | default arm colors                  |
| Shared x-tick helper              | endpoint + pretty() interior        |

Full `plot.causal_cr_risk()` wraps these. The wrappers stay in
`CausalCompetingRisks` (framework-specific arm labels, legend text,
contrast annotations), but the tile primitives live in `CausalKit`.

### 4.7 — Contrast computation (→ `CausalKit`, partially)

| Function                   | Verdict                                      |
|----------------------------|----------------------------------------------|
| (generic) percentile CI from 4D array → long data frame | → CausalKit |
| `compute_contrasts()` with A/B decomposition specific to separable effects | stays in CausalCompetingRisks |

Proposal: `CausalKit::bootstrap_percentile_ci(replicates, alpha)` —
takes a generic array of replicates, returns long-format CIs at each
time. `CausalCompetingRisks::compute_contrasts()` calls it for each
contrast (total, direct_A, etc.).

### 4.8 — Separable effects specifics (stays in `CausalCompetingRisks`)

- `causal_cr()` + `fit_causal_cr()`
- `add_d_swap_weights()` (swap semantics specific to separable effects)
- `weighted_arm_cum_inc_rep1()` / `_rep2()`
- `estimate_weighted_cum_inc_rep1()` / `_rep2()`
- `ipw_estimate()` orchestrator
- `gformula_estimate()` orchestrator (the arm enumeration for (1,1), (0,0), (1,0), (0,1) is separable-effects specific)
- `compute_contrasts()` (decomposition A/B)
- `risk()`, `contrast()`, `diagnostic()` accessors
- `bootstrap()` orchestrator
- Full `plot.causal_cr_risk()` / `plot.separable_effects_contrast()`
- All `print.*` / `summary.*` / `confint.*` methods

---

## 5. `CausalKit` Public API Sketch

```r
# Data prep
to_person_time(data, ...)
validate_input_shape(x, name)
validate_subject_level(data, ...)
validate_person_time(pt_data, ...)

# Hazards
fit_hazard_model(formula, data, label)    # single model; generalized
check_fitted_positivity(model, label)

# G-formula primitives
make_clone(baseline, cut_times, overrides = list())
predict_hazards(clone, models)
compute_cum_inc(clone, cut_times, id_col)

# IPW primitives
add_hazard_and_cumprod(pt_data, model, arm_overrides, label, id_col)
add_lagged_cumprod(pt_data, cumprod_col, id_col)
add_ipcw_weights(pt_data, model_c, id_col)
apply_weight_truncation(pt_data, weight_cols, id_col, adjust, threshold)
weighted_hazard_by_k(d, flag_col, weight_col)
cum_inc_from_weighted(d, cut_times)
summarize_weights(pt_data, weight_cols)

# Bootstrap
resample_subjects(pt_data, id_col, seed = NULL)
bootstrap_percentile_ci(replicates, alpha)

# Plot primitives
plot_cum_inc(ci_data_long, ...)           # generic step curves
plot_risk_table(counts_data, ...)         # generic table
build_arm_colors(arms, palette = "okabe_ito")

# Shared constants
OKABE_ITO
DEFAULT_EPS_POSITIVITY
```

---

## 6. Dependency Graph

```
CausalKit
    │  Imports: ggplot2, stats, survival
    │  Suggests: patchwork
    ▼
CausalCompetingRisks
    │  Imports: CausalKit (+ inherits its deps transitively)
    ▼
(user-facing package)
```

S3 class design:
- `"person_time"` — defined in `CausalKit` (used by `to_person_time()`). Downstream packages respect the class + attributes.
- `"causal_cr"`, `"causal_cr_risk"`, `"separable_effects_contrast"`, `"separable_effects_bootstrap"`, `"separable_effects_diagnostic"` — defined in `CausalCompetingRisks`.
- Any methods for the `_pt` class would be defined in `CausalKit`; methods for the domain-specific classes live in the respective downstream package.

**S3 gotcha:** if a generic is defined in `CausalKit` (e.g., `plot_cum_inc` as generic) and methods for it are defined in `CausalCompetingRisks`, R's S3 dispatch works but the downstream package must register the methods explicitly via `@exportS3Method`. Roxygen handles this.

---

## 7. `CausalCompetingRisks` After the Split

After extraction, `CausalCompetingRisks` becomes much smaller:
- Separable-effects-specific estimation (rep1/rep2, arm_01, decomposition A/B)
- Orchestration: `causal_cr()`, `fit_causal_cr()`, `bootstrap()`
- Accessors: `risk()`, `contrast()`, `diagnostic()`
- Domain-specific plot/print/summary methods
- Documentation tying everything back to Stensrud et al.

Rough estimate:
- Current lines of R/: ~3,000
- After split → CausalCompetingRisks: ~1,200
- After split → CausalKit: ~1,800

---

## 8. Migration Plan (When We're Ready)

Trigger: starting work on `CausalSeparableMediation` (or another consumer).

Steps:
1. Create `CausalKit` package skeleton next to `CausalCompetingRisks`.
2. Move the inventory items from §4.1–4.7 into `CausalKit/R/`.
3. Add `Imports: CausalKit` to `CausalCompetingRisks`'s DESCRIPTION.
4. Update internal calls in `CausalCompetingRisks` to use `CausalKit::` prefix (or `@importFrom CausalKit ...`).
5. Move tests that cover pure-primitive functions into `CausalKit/tests/`.
6. Keep integration tests in `CausalCompetingRisks/tests/`.
7. Run `R CMD check` on both packages.
8. Version bump: `CausalKit 0.1.0`, `CausalCompetingRisks 0.2.0`.
9. Update vignettes.

Estimated effort: ~1 day of focused work if the API stays stable (no tweaks during the move).

---

## 9. Open Questions

1. **Naming.** Which of the candidates in §2?
2. **Bootstrap API.** Generic `fit_fn`-parameterized in core, or framework-specific orchestration in each downstream package? (Lean: framework-specific orchestration; only `resample_subjects()` and `bootstrap_percentile_ci()` in core.)
3. **Plot primitives granularity.** Expose `plot_cum_inc()` as a full generic, or just helper builders that downstream packages compose? (Lean: helper builders.)
4. **Positivity threshold.** `CausalKit` ships a default `eps = 1e-6`, but downstream packages might override. Parameterize cleanly.
5. **Person-time format.** `"person_time"` class name is framework-neutral (good). But columns `A_y`, `A_d` are separable-effects-specific. Options:
   - `CausalKit` knows nothing about `A_y`/`A_d`; downstream packages add those columns after `to_person_time()`.
   - `CausalKit` provides a helper `add_working_treatment_copies(pt, treatment, names = c("A_y", "A_d"))`.
   - Current package: adds them in `to_person_time()` — but that's separable-effects-coupled. Worth revisiting when we actually split.

---

## 10. Commit Strategy When We Split

Single multi-file commit (or PR):
- Creates `CausalKit/` subdirectory with skeleton (DESCRIPTION, NAMESPACE, R/, tests/testthat/)
- Moves files with `git mv` so history is preserved
- Updates `CausalCompetingRisks/DESCRIPTION` to import CausalKit
- Updates tests
- `R CMD check` both packages green

Avoid: piecemeal moves over multiple commits — intermediate states don't compile.
