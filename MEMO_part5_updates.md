# Memo: Changes to apply to Part 5 HTML (implementation post)

From the code review session on `part5_run.R`. The validated R script is the reference — the Part 5 HTML should be updated to match.

---

## Naming changes

| Old (in current HTML) | New | Why |
|---|---|---|
| `subject_df$any_event <- 1L` | `event_indicator = 1L` | Clearer name; must be 1 for ALL subjects including censored (survSplit needs it to flag terminal rows) |
| `time_0` / `event_time` after survSplit | `tstart` / `tstop` | These are interval boundaries [tstart, tstop), not the original event time |
| `surv_d_ad`, `surv_d_ay` | `d_free_k_ad`, `d_free_k_ay` | These are interval-level P(D=0), NOT cumulative survival. `1 - predict()` returns one-interval D-free probability |
| `cum_surv_d_ad` etc. | `cs_surv_d_ad` | "cs" = cause-specific. Cumulative D-free conditional on being event-free. Not overall survival |
| `surv_cens` | `uncens_k` | Interval-level P(not censored) |
| `cum_surv_cens` | `cs_surv_cens` | Cumulative uncensored probability |
| `pred_times` | `time_points <- c(0, cut_times)` | Clearer; reuses existing `cut_times` |

## Structural changes

1. **Time convention**: `tstart = 0` for everyone, `event_time = dtime + 1`. Month 0 = enrollment, month 1 = first month at risk. Needed because 10 subjects have `dtime = 0` (died in first month) and survSplit requires `start < end`.

2. **Single survSplit**: no double expansion needed. Use `event_indicator = 1L` for all subjects, derive `y_event`, `d_event`, `c_event` from `event_indicator + event_type` after expansion.

3. **dplyr throughout**: `case_when` for event_type, `mutate` for temporal ordering, `rename`, `filter`, `group_by`. Base R `$<-` assignments replaced.

4. **Temporal ordering in one mutate**: `d_event = ifelse(c_event==1, NA, d_event)` then `y_event = ifelse(d_event==1, NA, y_event)` — works because mutate is sequential within a call.

## Key comments to add in the HTML

### At `compute_cum_inc`:
- This is the **parametric g-formula**, not Aalen-Johansen. Must compute per-subject then average. `cumprod(mean(surv)) ≠ mean(cumprod(surv))` because each subject has different covariates → different hazard trajectories → different cumulative survival paths.
- Code averages by k first then cumsums. Equation writes cumsum per subject then average. Same result by commutativity of the double sum. Code order is more memory-efficient (n×K → K before cumulating).

### At `surv = (1 - haz_y) * (1 - haz_d)`:
- Joint event-free probability within ONE interval, not cumulative. Becomes cumulative only after `cumprod()`.

### At the standardization step:
- Σ_l f(l)P(L=l) is replaced by (1/n)Σ_i f(l_i). Each subject carries their own covariates; the empirical distribution assigns weight 1/n. The cloned datasets make this concrete.

### At IPW weights:
- `d_free_k` is per-interval. `cs_surv` is the cumulative product = cause-specific survival conditional on being event-free (Y-free AND D-free). The Y-free conditioning is implicit in the hazard model's risk set.
- W_D reweights D-free trajectories: if D more likely under a_Y → surviving D is "too rare" → W_D > 1 upweights.

### At `estimate_sep_eff`:
- Divides by n at BASELINE, not current risk set. This is Horvitz-Thompson: the weights already correct for attrition.
- No Y model needed. Just observed Y events, reweighted. Trade-off: needs D and censoring models correct instead.

## Math traceback for the cumulative incidence formula

The formula in Part 5 Section 4 box comes from:
1. **Beyersmann et al. (2012)** eq. 4.8: F_h(t) = ∫ S(u) α_h(u) du (continuous subdistribution)
2. **Young et al. (2020) SiM**: discrete-time g-formula for competing events (eq. 3–5)
3. **Stensrud et al. (2020) JASA** eq. (8)/(9): extends to cross-arm separable effects
4. **Rojas-Saunero (SER 2023)**: code implementation via `calculateCumInc()`
