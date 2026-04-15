# Memo: Part 5 Review Notes

Notes from interactive code review session. Use these to update `part5-implementation.html` and `part5_run.R`.

---

## 1. Naming & Style Changes

- **tstart / tstop**: rename `time_0` → `tstart`, `event_time` → `tstop` after survSplit. These are interval boundaries, not subject-level event times.
- **k = tstart**: interval index derived from tstart.
- **time_points**: use `time_points <- c(0, cut_times)` for clone datasets instead of the unclear `pred_times`.
- **dplyr throughout**: use `mutate()`, `case_when()`, `group_by()`, `filter()`, `rename()`, `n_distinct()` — not base R assignment. Package will use base R, but the blog/script should use dplyr for readability.
- **case_when for event_type**: replace nested `ifelse()` with `case_when()`.
- **Covariate renaming**: make explicit mapping from raw prostate names (`pf`, `hx`, `hg`, `age`) to clean names (`normal_act`, `cv_hist`, `hemo_bin`, `age_cat`).

## 2. survSplit Details

- **Single survSplit call**: no double expansion needed. Use `event_indicator = 1L` for all subjects, then derive `y_event`, `d_event`, `c_event` from `event_indicator + event_type`.
- **Why event_indicator = 1L for everyone (including censored)**: survSplit uses this column to flag each subject's terminal row. If set to 0 for censored subjects, their last row is NOT flagged, and `c_event` cannot be derived. Comment this clearly in the code.
- **Time convention**: `tstart = 0` for everyone, `event_time = dtime + 1`. Month 0 = enrollment, month 1 = first month at risk. This avoids the workshop's `Tstart = -0.01` hack and ensures `event_time > tstart` (survSplit requires start < end). 10 subjects have dtime=0 in raw data — the +1 shift handles them.
- **Formula vs named-argument form**: `survSplit(Surv(tstart, event_time, event_indicator) ~ ., ...)` and the named-argument form produce identical results. The two-argument `Surv(event_time, event_indicator)` also works but creates `tstart` under a different name. Either is fine.
- **Do NOT use `Surv(event_type > 0)` shortcut**: this consumes the `event_type` column from the output, which we need afterward to derive y/d/c events.

## 3. Temporal Ordering (NA assignment)

The base R version works but should be converted to dplyr `mutate()` with `ifelse()`:
```r
pt_df <- pt_df %>%
  mutate(
    d_event = ifelse(c_event == 1, NA, d_event),
    y_event = ifelse(c_event == 1, NA, y_event),
    y_event = ifelse(d_event == 1, NA, y_event)  # sequential within mutate
  )
```
Note: within a single `mutate()` call, assignments are sequential — the second `y_event` line sees the already-updated `d_event` (with NAs from censoring). This preserves the required ordering.

## 4. compute_cum_inc: Why Per-Subject Then Average

Add a comment explaining why we can't average hazards first then cumulate (Aalen-Johansen style). The parametric g-formula requires per-subject computation because:
- Each subject has their own hazard trajectory (covariates differ)
- `cumprod(mean(surv))` ≠ `mean(cumprod(surv))` — the cumulative product doesn't commute with averaging
- The nonparametric Aalen-Johansen estimator works at the population level (average first), but the parametric g-formula standardizes over L by computing conditional quantities per subject, then marginalizing

Suggested comment block:
```r
# The g-formula standardizes over covariates L:
#   ν = Σ_l f(l) P(L=l) ≈ (1/n) Σ_i f(l_i)
#
# This requires per-subject cumulative survival paths because each subject
# has different covariates → different hazards → different survival trajectory.
# We CANNOT average hazards across subjects first then cumulate:
#   cumprod(mean(surv)) ≠ mean(cumprod(surv))
# The nonlinear cumprod doesn't commute with averaging.
#
# This differs from the nonparametric Aalen-Johansen estimator which operates
# at the population level. Here we do: conditional (per-subject) → marginalize.
```

## 5. Math Traceback for the Cumulative Incidence Formula

The iterative computation traces back through:

1. **Stensrud et al. (2020) eq. (9)**: the full g-formula in the hypothetical 4-arm trial where both $A_Y$ and $A_D$ are randomized
2. **Stensrud et al. (2020) eq. (8)**: applying the dismissible component conditions ($\Delta 1$, $\Delta 2$) gives the cross-arm form — Y-hazard from $A = a_Y$ arm, D-hazard from $A = a_D$ arm
3. **Stensrud et al. (2020) eq. (10)**: adding censoring ($\bar{C} = 0$ conditioning)
4. **Young et al. (2020) eq. (23)**: the general g-formula for direct effects under elimination of competing events (same structure, different estimand)
5. **Beyersmann et al. (2012) eq. (4.8)**: the continuous-time foundation $F_h(t) = \int_0^t S(u) \alpha_h(u) du$

The per-interval incidence increment maps from eq. (8) as follows:

At interval $k$, for subject $i$:
- $S_i(k-1) = \prod_{j=0}^{k-1} (1 - h_{Y,i}(j))(1 - h_{D,i}(j))$ — survived all prior intervals
- At interval $k$: D is resolved first (temporal ordering), then Y can happen
- $\text{inc}_i(k) = h_{Y,i}(k) \times (1 - h_{D,i}(k)) \times S_i(k-1)$

This matches the code: `inc = haz_y * (1 - haz_d) * lag_cum_surv`.

The standardization $\sum_l \ldots \Pr(L=l)$ is replaced by the sample average $\frac{1}{n}\sum_i$ because each subject provides one draw from the covariate distribution (empirical distribution assigns weight $1/n$ to each $l_i$).

## 6. Model Specification Note

Add a note that our additive specification (`A_y + k + I(k^2) + I(k^3)`) assumes a constant treatment effect over time, while the workshop uses interactions (`rx * (dtime + ...)`) allowing the effect to vary. Both are valid. Results will differ from workshop due to this modelling choice, not bugs.

## 7. surv Column

`clone$surv = (1 - clone$haz_y) * (1 - clone$haz_d)` — this is the joint event-free probability WITHIN one interval, not cumulative survival. It assumes Y and D are independent within each interval conditional on covariates and being alive. In discrete time with small intervals this is standard (probability of both events in one interval is negligible).

## 8. IPW Risk Set Filtering

`estimate_sep_eff` should use explicit `!is.na()` checks rather than relying on `NA == 0` returning `NA` (which happens to exclude correctly but is fragile):
```r
riskset_s <- ay_df %>%
  filter(k == s, !is.na(y_event), !is.na(d_event),
         c_event == 0, d_event == 0)
```

## 9. Censoring Model Note

The censoring model has extreme coefficients (intercept ≈ -249) because there's almost no censoring before late in follow-up. Workshop handles this by fitting only on `dtime > 50`. Our model works but weights get large (max W_C ≈ 30). Worth noting in the post but not a bug.

## 10. dplyr Version Issue

If `dplyr` < 1.1 is installed alongside `vctrs` >= 0.7, `mutate()` crashes with `vec_is_vector() is defunct`. Fix: update dplyr (`install.packages("dplyr")`), then restart R session. The error is internal to dplyr, not in user code. Also need to restart R if old `lifecycle` package is still loaded (`lifecycle >= 1.0.5` required).
