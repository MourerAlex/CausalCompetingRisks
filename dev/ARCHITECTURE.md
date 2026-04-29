# Architecture — `CausalCompetingRisks` v1

Captured 2026-04-27 during adversarial re-grill (4 parallel agents).
Companion to `dev/API_v1_LOCKED.md` (which fixes the user-facing surface).
This file fixes the **internal file layout and dependency structure**.

---

## Design principle

Separable effects are an **overlay** on classical estimators.

| Layer                    | Classical                      | Separable add-on                                |
|--------------------------|--------------------------------|-------------------------------------------------|
| Hazards                  | fit `λ_Y(A,L)`, `λ_D(A,L)`    | shared — same hazards                           |
| Arms                     | `a ∈ {0,1}` (2 arms)           | `(a_Y, a_D) ∈ {0,1}²` (4 arms)                  |
| Counterfactual prediction| set `A := a` everywhere        | set `A := a_Y` in λ_Y, `A := a_D` in λ_D        |
| IPW weights              | `1/π(A|L) × IPCW`              | + swap weights `W_D`, `W_Y` (Stensrud appendix) |
| CIF computation          | pooled-hazard recursion        | same recursion, different arm spec              |

Classical = core machinery. Separable = arm-specification + swap-weight overlay.
Files reflect this separation, so classical pieces remain extractable / reusable.

---

## File layout (R/ stays flat — CRAN convention; groups are conceptual)

```
─── 1 · CORE ESTIMATION (classical, reusable) ───────────────────────────
   hazards.R              fit/predict λ_Y, λ_D, λ_C; + cif_from_hazards()
   propensity.R           fit/predict π(A|L)
   weights.R              classical IPW weights (1/π × IPCW, stabilize, truncate)
   gformula_core.R        classical g-formula CIF per arm
   ipw_core.R             classical IPW CIF per arm

─── 2 · SEPARABLE OVERLAY ───────────────────────────────────────────────
   separable_arms.R       arm specification (a_Y, a_D), shadow column trick
   separable_swap.R       W_D and W_Y construction (Stensrud 2020 appendix)
   separable_gformula.R   uses gformula_core + separable_arms
   separable_ipw.R        uses ipw_core + separable_swap

─── 3 · INFERENCE ───────────────────────────────────────────────────────
   bootstrap.R            cluster-by-id resample, refit, percentile CI
   contrasts.R            decomposition A & B contrast formulas

─── 4 · VALIDATION ──────────────────────────────────────────────────────
   validate.R             strict locked validation

─── 5 · API SURFACE ─────────────────────────────────────────────────────
   separable_effects.R    entry point
   accessors.R            risk(), contrast(), assumptions()
   risk_table.R           risk_table()
   print.R / plot.R       S3 methods
```

**No file prefix conventions.** Naming alone signals group: `separable_*.R` is the overlay; everything else is core, inference, validation, or API.

15 files total. `data_prep.R` removed (single-input long-format only).

---

## Dependency graph (strict DAG, verified cycle-free)

```
LEAF (no inter-file deps)
  hazards.R · propensity.R · separable_arms.R · validate.R · contrasts.R

LEVEL 2
  weights.R           ←  hazards.R, propensity.R
  separable_swap.R    ←  hazards.R, separable_arms.R
  gformula_core.R     ←  hazards.R
  ipw_core.R          ←  weights.R, hazards.R

LEVEL 3 (separable estimators)
  separable_gformula.R  ←  gformula_core.R, separable_arms.R, hazards.R
  separable_ipw.R       ←  ipw_core.R, separable_swap.R, separable_arms.R

LEVEL 4 (entry)
  separable_effects.R   ←  validate.R, separable_gformula.R, separable_ipw.R, contrasts.R

LEVEL 5 (post-fit)
  bootstrap.R           ←  .fit_separable_effects() worker (internal, in separable_effects.R)
  accessors.R           ←  contrasts.R (+ reads fit slots)
  risk_table.R          ←  reads fit slots only

LEVEL 6 (S3 methods)
  print.R · plot.R      ←  accessors.R
```

Topological order: leaf → 2 → 3 → 4 → 5 → 6. No cycles.

---

## Subtleties already resolved

### Shared CIF recursion
Both `gformula_core.R` and `ipw_core.R` use the pooled-hazard CIF recursion
`cum_inc = 1 − prod(1 − hazard)`. To avoid duplication, this lives as
`cif_from_hazards()` inside `hazards.R` (or its own tiny file if it grows).

### `hazards.R` stays classical
Hazards never know about separable arms. They fit and predict via whatever
formula is supplied. Separable awareness lives in three places:
- formula construction (uses `A_Y`, `A_D` as predictors)
- `separable_arms.R` arm dispatcher (clones data, writes counterfactual `A_y`/`A_d`)
- `separable_swap.R` swap-weight construction

### Internal `.fit_separable_effects()` worker
`separable_effects()` is a thin wrapper: validate → call `.fit_separable_effects()`
→ return. Bootstrap calls the worker directly per replicate, skipping
re-validation. Avoids `bootstrap.R → separable_effects.R` looking like a
post-fit-calls-entry-point smell.

### Truncation as a default, not hardcoded
`separable_effects(..., truncate = c(0.01, 0.99))` — default per Cole & Hernán
2008, but user-overridable. `NULL` = no truncation.

### Warning-list scaffolding dropped, warnings preserved
The `withCallingHandlers()` wrapper that captures warnings into `fit$warnings`
and re-prints them is deleted. Plain `warning()` calls and the `$diagnostics`
slot remain. Positivity / convergence / weight diagnostics stay user-visible.

### Plot extensibility
`plot(fit, arms = NULL, colors = NULL, labels = NULL, show_ribbons = TRUE)`
returns a `ggplot` so users can `+ theme_bw()`, `+ labs(title = "...")`, etc.
Last `+` wins per ggplot convention; defaults are fully overridable.

---

## What is NOT in this layout (and why)

- **Subdirectories under `R/`** — CRAN-historic friction. Stay flat.
- **Class-per-file (S4-style)** — package is functional with a single `causal_cr_fit` S3.
- **Dependency injection (`fitter`, `predictor` callables)** — over-engineering for v1; we only ever call `glm()`. Drop.
- **Custom person-time S3 class (`pt_data`, `causal_cr_pt`)** — locked API takes long-format data.frame directly.
- **`data_prep.R`** — `to_person_time()` removed; users with subject-level data convert externally (e.g. `survival::survSplit`).
