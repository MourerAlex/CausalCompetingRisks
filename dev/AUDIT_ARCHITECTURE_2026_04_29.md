# Architecture audit — 2026-04-29

Adversarial pass on yesterday's 15-file restructure. Skills:
`improve-codebase-architecture`, `simplify`, `grill-claude-code`.

10 findings. 3 blocking (would force rewrites at v2/v3), 4 medium, 3 polish.

## Blocking for v1

### A1. Hardcoded arm names in dispatchers
- Files: `R/separable_ipw.R`, `R/bootstrap.R`, `R/contrasts.R`, `R/accessors.R`
- Issue: arm vector `c("arm_11", "arm_00", "arm_10", "arm_01")` is
  hardcoded in 3+ places. v3 multi-cause extension (D1, D2, …) needs
  arm space `(a_Y, a_D1, a_D2, …)` → arbitrary tuple. v2 partial
  isolation may need arm subset.
- TODO.md "Future structural extensions" already flags this design
  principle: keep arm-indexing in a dispatcher layer.
- Why blocking now: every accessor / plot / contrast we ship today
  hardcodes 4 arms. Each one becomes a v3 migration point.
- Suggested: single `.arm_space()` helper returning the arm list;
  dispatchers and bootstrap iterate over it. Arms become metadata on
  the fit object (`fit$arms`).

### A2. Formulas container does not close over `active_methods`
- Files: `R/separable_effects.R`, `R/fit_separable_effects()`
- Issue: user can pass `formulas$A` while requesting `method = "gformula"`;
  it's silently ignored. Also reverse: pass `formulas$Y` only and the
  IPW path will silently fall back to defaults for D, C, A.
- Suggested: validate `formulas` against `active_methods` and `ipcw`;
  warn on unused keys, error on missing required ones.

### A3. Method dispatch logic scattered
- Files: `R/separable_effects.R:103–319`, `R/hazards.R:68–87`
- Issue: "should we run IPW Rep 2?" is reconstructed via string matching
  on `method` in 5+ places. Adding v2 methods (per-protocol, time-varying)
  requires editing each site.
- Suggested: module constants `.VALID_METHODS`, `.METHOD_GROUPS`;
  predicate helpers `.needs_propensity()`, `.needs_d_swap()`,
  `.needs_y_swap()`. Dispatcher reads these.

## Paying debt soon (v1.x)

### B1. `R/ipw_core.R` is shallow
- 2 functions (`weighted_hazard_by_k`, `cum_inc_from_weighted`), ~80 LOC.
- These are framework-agnostic primitives — natural fit either inside
  `R/hazards.R` (CIF computation belongs near the hazard primitives) or
  reserved for the future `causal_tools` package noted in TODO.md.
- Suggested: merge into `hazards.R` for now; lift to shared package later.

### B2. `separable_swap.R` reaches into column names produced by `separable_arm_hazards()`
- Implicit cross-file contract on names like `cumprod_one_minus_hazd_a0`.
- Suggested: `.cumprod_col_name(arm, event)` helper; `@section Requirements`
  on `separable_arm_hazards()` listing the columns it produces.

### B3. Three-layer weight model not enforced
- Layer 2 returns `_raw` weights; Layer 3 post-processes. Contract is
  conceptual, not codified.
- Suggested: `_raw` suffix as documented invariant; layer 2 functions
  forbid returning post-truncated columns; unit-test the boundary.

### B4. Contrasts arm schema implicit
- Dispatcher returns 4 arms; downstream defensively checks
  `if ("arm_01" %in% names(...))`. Brittle once arm space becomes
  dynamic.
- Suggested: dispatcher returns `(arms_metadata, results_per_arm)`;
  downstream loops over metadata.

## Polish (v1.1+)

### C1. `R/separable_arms.R` is single-function (~69 LOC)
- Could merge with `separable_swap.R` into one `separable_weights.R`.
  Optional — current split has pedagogical value.

### C2. `R/risk_table.R` is single-function (~89 LOC)
- Same observation. Defer.

### C3. Slot naming, time-polynomial, ipcw forward-compat
- `fit$cumulative_incidence` long; `fit$risk` is the locked target.
- Time polynomial in default formulas is unexplained anywhere.
- ipcw flag forward-compatibility for `method = "ipw"` running both reps
  with one ipcw model — confirm in code path.

## Architectural strengths

- Core/overlay separation is clean. Dependency DAG acyclic.
- Three-layer weight conceptual model is sound (just under-enforced).
- Validation layer comprehensive.
- Bootstrap correctness OK (cluster-by-id done; arms 4-up).

## Recommended ordering

1. A1 (arm dispatch) — touches every downstream, fix before Issue 4
   (long-format `$risk`).
2. A3 (method dispatch helpers) — same reason.
3. A2 (formula validation) — small surface, low risk, high UX value.
4. B1–B4 — during the next refactor pass.
5. C1–C3 — opportunistic.
