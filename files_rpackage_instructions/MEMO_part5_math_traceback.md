# Memo: Math Traceback for compute_cum_inc

Add this derivation as a comment block or explanatory section near `compute_cum_inc()` in Part 5, showing exactly how the code maps to the paper.

---

## Source

**Stensrud et al. (2020), JASA**, Section 5.2, equations (8) and (9).

## The target

$$\Pr(Y_{k+1}^{a_Y, a_D} = 1)$$

The risk of Y by time $k+1$ under separable treatment $(a_Y, a_D)$.

## Equation (9): the g-formula in a hypothetical 4-arm trial

If both $A_Y$ and $A_D$ were randomized:

$$\Pr(Y_{k+1}^{a_Y, a_D} = 1) = \sum_l \Bigg[ \sum_{s=0}^{k} \Pr(Y_{s+1} = 1 \mid D_{s+1} = Y_s = 0, A_Y = a_Y, A_D = a_D, L = l)$$

$$\times \prod_{j=0}^{s} \Pr(D_{j+1} = 0 \mid D_j = Y_j = 0, A_Y = a_Y, A_D = a_D, L = l) \times \Pr(Y_j = 0 \mid D_j = Y_{j-1} = 0, A_Y = a_Y, A_D = a_D, L = l) \Bigg] \Pr(L = l)$$

## Equation (8): applying the dismissible component conditions

The dismissible component conditions ($\Delta 1$ and $\Delta 2$, p.13) allow us to decouple the treatment components:

- $\Delta 1$: the Y-hazard does not depend on $A_D$ → only needs $A = a_Y$
- $\Delta 2$: the D-hazard does not depend on $A_Y$ → only needs $A = a_D$

This gives the **cross-arm form** (eq. 8):

$$= \sum_l \Bigg[ \sum_{s=0}^{k} \underbrace{\Pr(Y_{s+1} = 1 \mid D_{s+1} = Y_s = 0, A = a_Y, L = l)}_{\text{Y-hazard from } a_Y \text{ arm}}$$

$$\times \prod_{j=0}^{s} \underbrace{\Pr(D_{j+1} = 0 \mid D_j = Y_j = 0, A = a_D, L = l)}_{\text{D event-free from } a_D \text{ arm}} \times \underbrace{\Pr(Y_j = 0 \mid \ldots, A = a_Y, L = l)}_{\text{Y event-free from } a_Y \text{ arm}} \Bigg] \Pr(L = l)$$

This is the key result: each conditional probability is identified from observed data under a **single** treatment value ($a_Y$ or $a_D$), not both. We only observe $A$, but the dismissible conditions let us read hazards at $A = a_Y$ for Y-related terms and $A = a_D$ for D-related terms separately.

## From equation to code

**Define shorthand** for subject $i$ with covariates $L = l_i$:

$$h_{Y,i}(k) = \Pr(Y_{k+1} = 1 \mid D_{k+1} = Y_k = 0, A = a_Y, L = l_i) \quad \text{(= haz\_y in code)}$$
$$h_{D,i}(k) = \Pr(D_{k+1} = 1 \mid D_k = Y_k = 0, A = a_D, L = l_i) \quad \text{(= haz\_d in code)}$$

**The product** $\prod_{j=0}^{s} (1 - h_{D,i}(j)) \times (1 - h_{Y,i}(j))$ is the joint event-free survival through interval $s$. In the paper, $j$ runs from 0 to $s$, so D-survival at $s+1$ is NOT in the product — it's in the conditioning of the Y-hazard ($D_{s+1} = 0$).

**Temporal ordering within interval $k$**: at interval $k$, D is resolved first, then Y. So the incidence increment at interval $k$ is:

1. Survive all prior intervals: $S_i(k-1) = \prod_{j=0}^{k-1} (1 - h_{Y,i}(j))(1 - h_{D,i}(j))$
2. At interval $k$: survive D (prob $1 - h_{D,i}(k)$), then experience Y (prob $h_{Y,i}(k)$)

$$\text{inc}_i(k) = h_{Y,i}(k) \times (1 - h_{D,i}(k)) \times S_i(k-1)$$

This maps directly to:
```r
inc = haz_y * (1 - haz_d) * lag_cum_surv
```
where `lag_cum_surv` = $S_i(k-1)$ via `lag(cumprod(surv), default = 1)`, and `surv` = $(1-h_Y)(1-h_D)$.

**The standardization** $\sum_l \ldots \Pr(L = l)$ becomes the sample average $\frac{1}{n}\sum_{i=1}^{n}$. Each subject $i$ provides one draw from the covariate distribution — the empirical distribution assigns weight $1/n$ to each observed $l_i$. This is why we compute per-subject first (each has a different hazard trajectory due to different covariates), then average across subjects:

$$\hat{\nu}_{a_Y, a_D}(k) = \frac{1}{n} \sum_{i=1}^{n} \sum_{s=0}^{k} \text{inc}_i(s)$$

In code: `mean(inc)` grouped by `k`, then `cumsum(mean_inc)`.

**Why not average first then cumulate?** The parametric g-formula requires per-subject computation because $\text{cumprod}(\text{mean}(\text{surv})) \neq \text{mean}(\text{cumprod}(\text{surv}))$ — the nonlinear cumulative product doesn't commute with averaging. This differs from the nonparametric Aalen-Johansen estimator, which operates at the population level using a single pooled survival estimate.

**Note on order of operations in code vs equation:** The equation writes $\frac{1}{n} \sum_i \sum_s \text{inc}_i(s)$ — cumsum per subject, then average. The code does it in reverse: `mean(inc)` by k first (average over subjects at each time), then `cumsum` (accumulate over time). Both are equivalent because the double sum commutes:

$$\frac{1}{n} \sum_{i=1}^{n} \sum_{s=0}^{k} \text{inc}_i(s) = \sum_{s=0}^{k} \frac{1}{n} \sum_{i=1}^{n} \text{inc}_i(s)$$

The code order is more memory-efficient: averaging first collapses $n \times K$ rows to $K$ rows before cumulating, so we never store the full per-subject cumulative incidence paths. The alternative (cumsum per subject inside `group_by(id)`, then average by k) is mathematically identical but holds $n \times K$ cumulative values in memory before averaging.
