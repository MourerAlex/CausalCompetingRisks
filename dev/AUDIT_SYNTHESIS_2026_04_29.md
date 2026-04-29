# Audit synthesis — 2026-04-29

Cross-cutting findings from the three parallel audits
(`AUDIT_UL_2026_04_29.md`, `AUDIT_ARCHITECTURE_2026_04_29.md`,
`AUDIT_ESTIMAND_2026_04_29.md`). Goal: turn findings into a re-prioritized
v1 issue list before contacting Stensrud.

## Find — themes that appeared in 2+ audits

### T1. The package never states what it estimates
- Estimand audit §1, UL §1.
- Fixing this is one accessor (`assumptions()`) plus one roxygen
  section, but it cascades: print, summary, vignette all reference it.

### T2. The package never states under what conditions the estimate identifies the estimand
- Estimand audit §2, UL §1.
- Same accessor as T1. Add the dismissible component conditions Δ1/Δ2,
  generalized decomposition assumption, exchangeability, positivity,
  consistency, independent censoring.

### T3. Isolation regime is computed but not surfaced
- Estimand audit §3, partly UL §2 (effects called "direct"/"indirect"
  without isolation qualifier).
- Architecture is ready (w_d, w_y on pt_data). Need a read-out in
  `diagnostic()` and a regime classifier.

### T4. Arm space is hardcoded — blocks v3, also leaks into UL
- Architecture A1, B4; estimand audit §5 (decomposition A/B logic
  hardcodes which arm is the reference).
- Single fix: arms-as-metadata pattern; dispatchers iterate over
  `fit$arms`.

### T5. Method dispatch is string-matched — blocks v2, contributes to UL drift
- Architecture A3; UL §3 (`ipw1`/`ipw2` vs `Rep 1`/`Rep 2`).
- Helpers `.needs_propensity()`, `.needs_d_swap()`, `.needs_y_swap()`
  + collapse user token to `"ipw"` per locked API.

### T6. Renames from yesterday are half-done
- UL §4 (`causal_cr_*` classes), §6 (`n_boot` vs `B`), §11
  (`method = "all"` vs spec `"gformula"`); architecture mentions slot
  naming drift in C3.
- Mechanical sweep.

### T7. Positivity story is incomplete
- Estimand audit §6, UL §10 (`w_a` not glossed; weight columns
  unexplained).
- Add row-level flagging in Rep 2; rename diagnostic category to
  `positivity`; gloss weight columns.

## Target — re-prioritized v1 issue list

Replaces the order in `dev/SESSION_2026_04_28_HANDOFF.md` and the bullet
sequence I gave earlier today. Numbering is execution order, not
locked-API issue numbers.

### P0 — Stensrud-bar (do before any contact)
1. **Estimand + assumptions accessor**
   - Implement `assumptions(fit)` printing target counterfactual in
     notation + identifying assumptions (Δ1, Δ2, GDA, exchangeability,
     positivity, consistency, censoring).
   - Add `@section Identification` to `separable_effects()` roxygen.
   - Reference: estimand audit §1–2; UL #1.
2. **Arms-as-metadata refactor**
   - Single source of truth for arm space. Dispatchers iterate.
   - Reference: architecture A1, B4.
3. **Isolation diagnostic in `diagnostic(fit)`**
   - Summary stats on w_d, w_y; one-line regime classification.
   - Reference: estimand audit §3.
4. **Class rename + UL sweep**
   - All `causal_cr_*` → `separable_effects_*`.
   - Bare "direct"/"indirect" → "separable direct/indirect" everywhere.
   - `n_boot` → `B`.
   - `method = "all"` decision (recommend update locked spec, keep code).
   - Reference: UL #2, #4, #6, #11.

### P1 — Surface API completeness (locked-API alignment)
5. **`censoring_weights` → `ipcw` rename**
   - Sweep arg, roxygen, internal calls, slot.
6. **`method = c("gformula", "ipw")` collapse**
   - User-facing `"ipw"` runs both reps internally; output keys
     `ipw_rep1`/`ipw_rep2`.
   - Reference: architecture A3; UL #3.
7. **Long-format `$risk`** (Issue 4 in handoff)
   - Now safe to do — arms metadata in place.
8. **Method-dispatch helpers**
   - `.needs_propensity()`, `.needs_d_swap()`, `.needs_y_swap()`,
     `.needs_ipcw()`. Replace scattered string checks.
9. **Formula validation**
   - Warn on unused keys; error on missing required.
   - Reference: architecture A2.

### P2 — Output surface
10. **Accessors + S3 surface** (Issue 6 in handoff)
    - `risk()`, `contrast()`, `summary()`, `print()`, `assumptions()`,
      `risk_table()`, `diagnostic()`. Read from long-format `$risk`.
    - `contrast()` requires `method` arg (no hidden default).
    - `summary()` shows contrasts only when bootstrap attached.
11. **Bootstrap as standalone object** (Issue 5)
    - `bootstrap(fit, B, alpha)` returns
      `separable_effects_bootstrap`. Print method.
    - Inherits `truncate` from fit; reapplies per replicate.
12. **Positivity story**
    - Row-level flagging in Rep 2 when `haz_y < eps`.
    - Rename warning category to `positivity`.
    - Gloss `w_*` columns in roxygen.
    - Reference: estimand audit §6; UL #10.
13. **Plot methods**
    - `plot.separable_effects(risk(fit), ci = boot, method)`.
    - `plot.separable_effects_contrast(contrast(...))` with
      Okabe-Ito palette.
    - `plot.separable_effects_diagnostic(diagnostic(fit))`.
    - Risk table beneath plots.
    - Unify contrast labels with what `contrast()` returns
      (UL #8 cleanup).

### P3 — Numerical truth (deferred per user; end of v1)
14. **Controlled DGP harness**
    - `simulate_separable_effects_data()` with known true effects.
    - Unit tests: Rep 1 ≡ Rep 2 under full isolation; arm CIFs
      converge to truth as n grows; bootstrap CI coverage.
15. **Prostate end-to-end**
    - Snapshot fixture matching `part5_run.R` numbers.

### P4 — Build / vignette / Stensrud handoff
16. **Roxygen polish** (Hernán voice for prose).
17. **README** with quick-start.
18. **Walkthrough vignette** on prostate data.
19. **Interpretation vignette** with copy-pasteable phrasings.
20. **`R CMD check`** clean.
21. **Email Stensrud** with: README link, vignette link, issue list of
    open methodological choices for his input
    (decomposition A vs B default; isolation diagnostic threshold;
    Rep 1 vs Rep 2 reporting convention; v2 roadmap for Z_k partition).

## Propose — execution shape

Three observations on how to run this:

**(a) P0 is one coherent chunk.** Items 1–4 share the same files (entry
point, roxygen, accessors). Do them in one pass with one PR-equivalent
commit batch.

**(b) Defer items 7 (long-format `$risk`) and 13 (plot methods) until
after P0.** Both depend on the arms-metadata refactor (item 2) and the
class rename (item 4). Doing them first means redoing them.

**(c) Numerical truth (P3) at the end is the right call** if and only
if the controlled DGP confirms the math is sound. If it surfaces a
discrepancy (e.g., Rep 1 ≢ Rep 2 under exact full isolation on
synthetic data), that's a P0 escalation and unrolls the rest.

## Open questions for next session

1. `assumptions()` output format — plain text block, or structured
   list with `format()` method for re-use in vignette?
2. Isolation regime threshold — what counts as "w_d ≈ 1"? Quantile?
   Mean absolute deviation? User-tunable?
3. Decomposition A vs B in default print — both, or A with B opt-in?
4. Stensrud first-contact letter — short (vignette + README + ask)
   or long (vignette + audit synthesis + roadmap)?
