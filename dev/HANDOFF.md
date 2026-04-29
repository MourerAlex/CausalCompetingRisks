# HANDOFF: Causal Inference Blog Series — Complete State

## PROJECT OVERVIEW

Alex Mourer writes a pedagogical blog series at **alexmourer.com** on causal inference with competing events, targeting epidemiologists with bachelor-level math. Core focus: the **separable effects framework** from Stensrud et al. (2020, 2021). Adjacent deep dives cover related causal mediation methodology.

Alex is a researcher/writer. Strong familiarity with:
- The separable effects framework
- The g-formula (parametric, pooled logistic regression)
- Discrete-time survival methods
- NPSEM-IE vs FFRCISTG distinction (Hernán & Robins Ch. 23)
- Competing events methodology broadly

Dual goal: produce high-quality blog content AND potentially develop new R software to fill identified gaps in the ecosystem (the discrete-time pooled logistic separable effects approach is not in any existing package).

Location: Paris. Sometimes writes in French. English is fine.

---

## COMMUNICATION STYLE

- **Terse.** One-line messages common. Pushes back when explanations lack rigor.
- **Drives toward refined explanations.** When Claude gives an imprecise answer, Alex says so and expects a correction, not an apology loop.
- **Precision over hedging.** Correctness and clarity beat verbose qualifications.
- **Dislikes excessive apology or self-abasement.** Fix and move on.
- **Will ask "why" repeatedly** to drill into any claim. Welcome this — don't paper over with hand-waving.
- **Sometimes asks to stop implementing changes without confirming first.** When Alex says "stop updating html without asking" — respect that. Draft in text first, get approval, then apply.

---

## BLOG STYLE CONVENTIONS

### Tone
- **Impersonal third person** throughout. Never "you," "your," or "we" — use "the reader," passive voice, or just describe what's happening.
- Retain all mathematical superscripts including `\bar{c}=0`.
- No emojis. Em-dashes are fine in math/prose when functional (—).

### Structure
- Each post has a consistent HTML template (see structure below).
- Section numbers use `<h2 id="sN">`. Subsections use `<h3>`.
- Variable glossary in `<details>` collapsible at top of post.
- Contents box at top (green-left-bordered `#f7f6f3` background).
- Equations in `eq-box` class divs (green left border, `#2a5f4f`, light bg `#f7f6f3`).
- Collapsible `<details>` blocks for derivations, proofs, or code the reader might skip.

### Mobile formatting
- Plain notation preferred over nested LaTeX when possible.
- Long equations broken across multiple `$$...$$` lines.
- Tables use `table-layout: fixed` with explicit `<colgroup>` widths.

### Color coding in math
- `\color{teal}` for Y-pathway / $a_Y$
- `\color{orange}` for D-pathway / $a_D$

### MathJax pitfalls encountered
- `\texttt{}` **breaks MathJax**. Use `\text{}` instead, OR keep math pure and map to code variable names in prose below the equation.
- Underscores in `\text{var_name}` cause subscript parsing errors. If you need variable names in math, write them without underscores or use escaped `\_` — but even that can fail. **Safest: keep code names out of math entirely.**

---

## BLOG POST SERIES STRUCTURE

All completed and live at `/mnt/user-data/outputs/`. Files renamed over time:

### Part 1: Competing Events & Causal Inference [COMPLETED]
**File:** `competing-events-v3.html` (also exists as `competing-events-causal-inference.html`)
**Based on:** Young et al. (2020), Stat Med — estimands paper
**Sections:**
1. Risk and hazard without competing events
2. The competing event enters
3. Three hazards, zero causal effects
4. Is the competing event a censoring event?
5. Identification: from counterfactual to observed data
6. The limits of this framework

### Part 2: Separable Effects [COMPLETED]
**File:** `part2-separable-effects-2020.html`
**Based on:** Stensrud et al. (2020), JASA — the 2020 paper
**Sections:**
1. A different question
2. Why this solves the ambiguity
3. The decomposition assumption
4. The dismissible component conditions
5. How identification works
6. Limitations and forward pointers

### Part 3: Generalized Separable Effects [COMPLETED]
**File:** `part3-generalized-separable-effects-2021.html`
**Based on:** Stensrud et al. (2021), Lifetime Data Analysis — the 2021 generalization
**Sections:**
1. What changes from 2020
2. $Z_k$ vs $L_k$: what you need vs what you have
3. The isolation hierarchy
4. Generalized dismissible conditions
5. The generalized g-formula
6. Discussion

### Part 4a: Technical Companion — 2020 Mathematics [COMPLETED]
**File:** `part4a-technical-2020.html`
**Based on:** Stensrud et al. (2020) Appendix + Young et al. (2020) Appendix
**Sections:**
1. The problem
2. The target quantity
3. The assumptions
4. From the target to the g-formula
5. From the g-formula to IPW Representation 1: the swap
6. IPW Representation 2: the mirror
7. Summary: what each estimator needs
8. The deep proof: counterfactuals to observed data
9. Where the dismissible conditions come from
10. When the conditions hold and fail

### Part 4b: Technical Companion — 2021 Extensions [COMPLETED]
**File:** `part4b-technical-2021.html`
**Based on:** Stensrud et al. (2021) Appendix
**Sections:**
1. What changes from 2020
2. The new data structure
3. The four dismissible conditions
4. The generalized g-formula
5. From the g-formula to IPW Representation 1: the new weight $W_{L_D}$
6. IPW Representation 2: the mirror with $W_{L_Y}$
7. Summary: what each estimator needs
8. The isolation hierarchy and the partition
9. Connection to path-specific effects (TODO placeholder)

### Part 5: Implementation — From Formula to Code [COMPLETED THIS SESSION]
**File:** `part5-implementation.html`
**Companion R script:** `part5-separable-effects.R`
**Based on:** Rojas-Saunero (SER 2023 Workshop) + Stensrud SPRINT code
**Sections:**
1. The data
2. The person-time dataset
3. A note on censoring
4. The g-formula estimator
5. The IPW estimator
6. Causal contrasts
7. Summary

---

## PROJECT FILES IN `/mnt/project/`

Reference PDFs:
- `1901.09472v3.pdf` — Stensrud 2020 separable effects paper
- `s10985021095308.pdf` — Stensrud 2021 generalized paper
- `Conditional_Separable_Effects.pdf` — Stensrud 2023
- `s10985023095948.pdf` — Janvin et al. 2024 (recurrent events)
- `nihms1659396.pdf` — Young et al. 2020 (estimands)
- `didelez2018.pdf` — Didelez 2019 (mediation precursor)
- `2025-11-21-causal-inference-what-if-hernan.pdf` — Hernán & Robins What If book
- `9780429029684_webpdf.pdf` — Andersen & Ravn book
- `kwab029.pdf`, `martinussen2021.pdf`, `shpitser2013.pdf`, `annurevstatistics040522094556.pdf`, etc.

Reference code:
- `gform.txt` — Rojas-Saunero's DES Rmd file, g-formula version (pooled logistic with `rx`/`Orx` trick, NO censoring model — censoring handled by conditioning)
- `ipw.txt` — Rojas-Saunero's DES Rmd file, IPW version (uses `cumPredO_1`/`cumPredO_0` for W_D, `1/cumPredC` for W_C)
- `10985_2021_9530_MOESM1_ESM.r` — Stensrud's SPRINT code (separate models per arm, bootstrap, uses `ipw_y` Rep 2 with lag trick)
- `utility_functions.R` — DES IPW helper functions

---

## KEY TECHNICAL DECISIONS FOR PART 5

### Data naming conventions (CRITICAL — stick to these)
```r
# Columns in subject_df (one row per subject)
id             # subject identifier
A              # observed treatment (1 = DES, 0 = placebo) — NOT "trt" or "rx"
event_time     # month of event or censoring
event_type     # 0 = censored, 1 = prostate death (Y), 2 = other death (D)
normal_act     # baseline: normal daily activity (0/1)
age_cat        # baseline: age category
cv_hist        # baseline: cardiovascular history (0/1) — NOT "hx"
hemo_bin       # baseline: hemoglobin < 12 g/dL (0/1) — NOT "hgBinary"

# Columns added in pt_df
k              # time interval (0, 1, ..., 59) — NOT "dtime"
A_y            # copy of A for Y-hazard model — NOT "rx" or "trt"
A_d            # copy of A for D-hazard model — NOT "Orx" or "trt_d"
y_event        # Y at end of this interval (0/1/NA)
d_event        # D at end of this interval (0/1/NA)
c_event        # censored at end of this interval (0/1)

# Model objects
fit_y_haz, fit_d_haz, fit_cens

# Predictions
surv_d_ad      # Pr(D=0|..., A_d = a_d) per interval
surv_d_ay      # Pr(D=0|..., A_d = a_y) per interval
surv_cens      # Pr(C=0|...) per interval

# Cumulative products
cum_surv_d_ad, cum_surv_d_ay, cum_surv_cens

# Weights
w_d, w_cens, w_total

# G-formula clones
clone_11, clone_00, clone_10

# In clones
haz_y, haz_d, surv, inc, cum_inc

# In riskset_s / eventset_s (in estimate_sep_eff)
riskset_s      # risk set at time s — NOT "at_s"
eventset_s     # subjects with Y event at time s — NOT "y_events"

# Scalars
a_y <- 1       # A_Y intervention value
a_d <- 0       # A_D intervention value
n              # number of subjects
```

### Model formulas (SIMPLE ADDITIVE, no interactions)
After iteration, Alex chose the simple form:
```r
y_event ~ A_y + k + I(k^2) + I(k^3) + normal_act + age_cat + cv_hist + hemo_bin
```
NOT `A_y * (k + I(k^2) + I(k^3))`. The treatment effect $\theta_A$ is constant over time.

Corresponding math in eq-box:
$$\text{logit}(p) = \alpha_0 + \alpha_1 k + \alpha_2 k^2 + \alpha_3 k^3 + \theta_A \cdot A + \theta_L' L$$

One sentence notes interactions can be added (as in Rojas-Saunero) if desired.

### The `A_y`/`A_d` trick explained
In observed `pt_df`: `A_y = A_d = A` (all equal). Purpose is to give each hazard model its own "dial":
- `fit_y_haz` uses `A_y`
- `fit_d_haz` uses `A_d`
- `fit_cens` uses `A` (observed, no override needed since W_C only applied in a_Y arm)

At prediction time in clones (or via `mutate` for IPW), they diverge:
- `clone_10`: `A_y = 1`, `A_d = 0` — cross-arm prediction

### Censoring handling — CRITICAL asymmetry
- **G-formula:** conditioning (no censoring model). Censored individuals drop out of risk set naturally via `glm`'s NA handling. Under assumption E2 (censoring ignorable given covariates in model), this is correct. Reference: Hernán & Robins TP 21.10.
- **IPW:** weighting via W_C. Individuals disappear from the data when censored; W_C upweights uncensored to represent them.

**The censoring model in IPW uses observed `A`** — NOT `A_y`. Reasoning: W_C is only ever applied inside the `a_Y` arm (the `estimate_sep_eff` function filters there). In that arm, `A = a_Y` already, so predicting at observed `A` IS predicting at `a_Y`. No override needed. Matches Rojas-Saunero and Stensrud reference code.

### Iterative cumulative incidence formulation
Instead of the closed-form formula with big product, show iterative buildup:
$$S_i(0) = 1$$
$$S_i(s) = S_i(s-1) \times (1 - \hat{p}_Y(s)) \times (1 - \hat{p}_D(s))$$
$$\text{inc}_i(s+1) = \hat{p}_Y(s+1) \times (1 - \hat{p}_D(s+1)) \times S_i(s)$$
$$\hat{\nu} = \frac{1}{n}\sum_i \sum_s \text{inc}_i(s+1)$$

Maps directly to code: `cumprod(surv)` builds `S_i(s)`, `lag(cum_surv, default = 1)` gives `S_i(s-1)`.

### W_D vs W_Y asymmetry (important for future discussion)
- **W_D** is simple: D appears ONCE per time step in the g-formula (Block B). One ratio, one cumulative product.
- **W_Y** is complex: Y appears TWICE per time step (Block A hazard at $s+1$, Block C event-free at $j$). Because $\frac{p}{p'} \neq \frac{1-p}{1-p'}$, they can't be combined. W_Y = (current hazard ratio) × (cumulative product of past event-free ratios). This is the "lag trick" in Stensrud SPRINT code.

This asymmetry explains why Part 4a IPW Rep 2 has a front term outside the product — it's the current hazard ratio at $s+1$, while the product sweeps past event-free terms $j = 0, \ldots, s-1$.

### Causal contrasts section
After getting curves $\hat{\nu}(1,1)$, $\hat{\nu}(0,0)$, $\hat{\nu}(1,0)$:
- **Total effect** = $\hat{\nu}(1,1) - \hat{\nu}(0,0)$
- **Separable indirect effect (A_D)** = $\hat{\nu}(1,1) - \hat{\nu}(1,0)$
- **Separable direct effect (A_Y)** = Total − Indirect

---

## PART 5 CURRENT STATE (as of end of this session)

### HTML file: `part5-implementation.html` (~3,400 words)

**Current 7 sections:**
1. **The data** — raw one-row-per-subject table shown first with 3 patients (A, B, C). `event_type` 0/1/2 as integers with mapping explained in prose below the table (not in table).
2. **The person-time dataset** — `survSplit`, event type indicators with NA rules, person-time table showing same 3 patients, `A_y`/`A_d` copies defined, `baseline_df` extracted.
3. **A note on censoring** — eq-box showing target $\Pr(Y^{a_Y, a_D, \bar{c}=0}_{K+1} = 1)$, then paragraphs explaining:
   - G-formula: conditioning mechanism, no censoring model, hazard models must include covariates predicting both outcome and censoring
   - IPW: weighting mechanism, W_C upweights uncensored
   - Both require E2
4. **The g-formula estimator** — three-piece bridge (target → identification → estimation):
   - Target statement
   - Identification eq-box (full g-formula with Blocks A, B, C labeled)
   - Pooled logistic model eq-box (simple additive form)
   - Prose explaining pooled fitting, why polynomial in k works
   - 4.1: Fit models (simple formulas, no interactions)
   - 4.2: `make_clone` function
   - 4.3: Cross-arm prediction eq-box (pure math $A = a_Y$ vs $A = a_D$, code mapping in prose)
   - 4.4: Iterative cumulative incidence eq-box + `compute_cum_inc` function
5. **The IPW estimator** — motivated by "in the a_Y arm, Blocks A and C are already correct, only Block B needs swapping":
   - IPW Rep 1 eq-box (full weighted expectation + W_D, W_C definitions)
   - 5.1: Fit D model + censoring model
   - 5.2: Pure-math eq-box for W_D numerator/denominator, prose maps to `surv_d_ad`/`surv_d_ay`
   - 5.3: Cumulative products + weights (with `group_by(id)` required!)
   - 5.4: `estimate_sep_eff` function using `riskset_s`/`eventset_s`, including the `w_total = w_cens` override for observed arms
6. **Causal contrasts** — eq-box with total/indirect/direct definitions, code computing risk differences at month 60
7. **Summary** — table (G-formula vs IPW Rep 1) with columns: Models needed, Expectation over, Censoring handling. Plus forward pointer to time-varying covariates (W_LD) future post.

**Plots:** REMOVED from HTML (Alex said "delete the plots tabs"). Plots are in the R file only.

**IPW Rep 2:** REMOVED from HTML. Briefly mentioned in Summary as "mirror option" pointing to Part 4a Section 6.

### R file: `part5-separable-effects.R` (~410 lines)

Fully runnable script:
- Simulated DGP: 250 subjects, baseline covariates mimicking DES trial, competing events via discrete-time Bernoulli draws with temporal ordering (C, D, Y)
- All code from HTML in runnable form
- Three plots: g-formula curves, IPW curves, g-formula vs IPW comparison on `clone_10` / cross-arm
- Causal contrasts printed for both estimators
- Comments matching HTML sections

### Specific comment style used in Part 5 (IMPORTANT)
Two-line comments explaining W_D logic:
```r
# D survival at the rate we WANT: placebo cardiovascular death rate
# This is the numerator of W_D — the rate the formula prescribes for Block B
pt_df$surv_d_ad <- 1 - predict(
  fit_d_haz,
  newdata = pt_df %>% mutate(A_d = a_d),  # a_d = 0
  type = "response"
)
```

Scalars used (`A_d = a_d`, `A_d = a_y`) not literals (`A_d = 0`, `A_d = 1`) — defined at top of script.

---

## KEY REVISIONS APPLIED DURING PART 5 SESSION

1. Variable name evolution: `rx` → `trt` → `A`, `orx` → `trt_d` → `A_d`, added `A_y`. Also `y_k1` → `y_event`, `dtime` → `k`, `hx` → `cv_hist`, `hgBinary` → `hemo_bin`.
2. `transform(df, rx = a_y)` → `df %>% mutate(A_y = a_y)` throughout.
3. Removed separate Section 1 on pooled logistic — merged into the g-formula section 4 as part of the three-piece bridge.
4. Fixed MathJax issues: removed `\texttt{}` (breaks), removed code variable names with underscores from equations, kept math pure and mapped to code in prose.
5. Removed plot collapsible tabs from HTML (still in R file).
6. Simplified model formulas from `A_y * (k + I(k^2) + I(k^3))` (which expands to interactions) to `A_y + k + I(k^2) + I(k^3)` (simple additive). Math simplified from $(\theta_1 + \theta_2 k + \theta_3 k^2 + \theta_4 k^3) \cdot A$ to $\theta_A \cdot A$.
7. Raw person-level data shown first (Section 1), then person-time (Section 2). Earlier version had person-time first.
8. Added censoring note as its own Section 3 explaining the g-formula vs IPW asymmetry in handling.
9. Added Causal Contrasts section (Section 6) with total/direct/indirect risk differences.
10. Iterative cumulative incidence formulation instead of closed-form product.
11. `is.na()` checks removed from `estimate_sep_eff` (shown to be redundant given NA rules).
12. Renamed `at_s` → `riskset_s`, `y_events` → `eventset_s`.
13. Added `group_by(id)` comment emphasizing cumulative products must be per subject.
14. Explicit `return()` in all custom functions (`make_clone`, `predict_hazards`, `compute_cum_inc`, `estimate_sep_eff`).
15. Added two-line comments on W_D numerator/denominator explaining "rate we WANT" vs "rate we CANCEL" with reference to Block B.

---

## TECHNICAL DEEP DIVES FROM PRIOR SESSIONS (Parts 1-4 context)

### 2020 Paper Appendix (all covered in Part 4a)
- Appendix A: Simulation DGP, subpopulations $Q_k$ and $R_k$
- Appendix B: SWIG independencies → dismissible conditions
- Appendix C: G-formula proof with multiply-and-divide ratio trick
- Appendix D: IPW proof via chain rule — W'_C cancels censoring, W_D cancels wrong D arm's rate
- Appendix E: When conditions hold/fail

### 2021 Paper Appendix (all covered in Part 4b)
- Appendix A: Modified treatment assumption variants
- Appendix B: Proof of identifiability — 10 assumptions, Lemmas 1-2, Theorem 1
- Appendix C: $Z_k$ partition — Lemmas 3-8, isolation hierarchy levels
- Appendix D: Bayes' rule trick for $W_{L_D}$ and $W_{L_Y}$
- Appendix E: 6-step estimation algorithm for IPW Rep 1
- Appendix F: Sensitivity function $t_k$ for violations of Δ1

### Key principles reinforced over many sessions
- Censoring handled in g-formula through covariate conditioning (no explicit censoring model)
- IPW requires explicit censoring model because it operates on actually observed individuals
- The censoring model in IPW uses observed `A` — `W_C` is only applied within the `a_Y` arm
- `group_by(id)` MUST come before `cumprod` to prevent bleeding across subjects
- Randomization renders baseline confounders irrelevant for the dismissible conditions argument
- The Lin/Zheng survival mediational g-formula and Didelez's approach require careful distinction
- **Recanting twins (Vo et al. 2025/2026) = stochastic replacement of intermediate confounder values. Alex's view: "feels like cheating."**

---

## R PACKAGE LANDSCAPE (investigated during Part 5)

**What exists:**
- `mets` (Martinussen & Scheike): implements separable effects in **continuous time** via Cox models + influence functions. Based on Martinussen & Stensrud 2023 Biometrics. `mediatorSurv` function.
- `CMAverse`: natural direct/indirect effects, not separable effects. Cross-world counterfactuals.
- `intmed`: interventional effects for mediation (VanderWeele/Vansteelandt/Robins). Different estimand.
- `medflex`: natural effect models, no competing events.
- `gfoRmula`: standard g-formula for time-varying treatments. Does not support the cross-arm prediction trick.
- `crumble`, `medoutcon`: cross-sectional only.

**What's missing:**
- No package implements the **discrete-time pooled logistic separable effects approach** from Stensrud 2020/2021.
- No package implements the $W_D$ swap weight with the `A_y`/`A_d` trick.
- No package for the $W_{L_D}$ Bayes' rule trick from Stensrud 2021.
- `gfoRmula` could in principle be hacked but the cross-arm prediction isn't a standard use case.

**This is the gap:** a potential R package could be built on the Stensrud SPRINT code base implementing the discrete-time version.

---

## HERNÁN & ROBINS ON MODEL MISSPECIFICATION CHECKS

From What If book (retrieved this session):

**Primary diagnostic:** Compare g-formula and IPW estimates. If they differ beyond sampling variability (bootstrap the difference), at least one set of models is misspecified. Book quote: "An implication is that one should always estimate E[Y^ā] using both methods and, if the estimates differ substantially (according to some prespecified criterion), reexamine all the models and modify them where necessary."

**Caveat:** "There is no logical guarantee of no model misspecification even when the estimates from both parametric approaches are similar, as they may both be biased in the same direction."

**G-null paradox (TP 21.3):** With time-varying confounders, non-saturated parametric models for outcome and covariate distributions cannot all be correctly specified simultaneously. Theoretical concern; in practice bias appears small relative to sampling variability.

**Doubly robust estimators** are the next step — consistent if either outcome or weight models are correct. Stensrud 2023 conditional separable effects paper derives one.

---

## PENDING / NOT YET WRITTEN

1. **Part 6 (future):** Time-varying covariates implementation
   - Monte Carlo g-formula
   - ICE (iterative conditional expectation) g-formula
   - $W_{L_D}$ with Bayes' rule trick
   - Sensitivity analysis for violations of $\Delta 1$
   - Bootstrap confidence intervals
2. **TODO placeholders** to fill:
   - Part 4a Section 9 (where dismissible conditions come from — needs SWIG diagrams)
   - Part 4b Section 9 (path-specific effects)
3. **Potential R package** for discrete-time separable effects (identified gap)
4. **Potentially:** Claude Code project for real R implementation on DES dataset

---

## WORKFLOW NOTES

- **Outputs go in `/mnt/user-data/outputs/`.** Alex can download from there.
- **Never implement HTML changes without showing draft text first** (Alex said this explicitly).
- **`view` tool first, then `str_replace`.** Don't guess at exact text.
- **When in doubt, ask Alex to confirm** before large structural changes.
- **Keep responses terse.** Alex's messages are often one line; long responses from Claude feel excessive.

---

## CURRENT IMMEDIATE CONTEXT

Last thing discussed: the W_D vs W_Y asymmetry in Part 4a Section 6 (IPW Rep 2). Specifically why the front term is OUTSIDE the product:
- Y appears in two places in g-formula: hazard at $s+1$ (ONCE) + event-free at $j=0,...,s-1$ (PRODUCT)
- Front ratio = current hazard swap (one term)
- Product = past event-free swaps (cumulative product)
- Not because of W_Y construction choice — it reflects Y's actual structure in the g-formula

Conversation ended because it was getting laggy — lots of iteration history. All Part 5 work is complete and saved.
