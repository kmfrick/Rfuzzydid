# Rfuzzydid: Fuzzy Difference-in-Differences

## Title

**fuzzydid** — Estimation with Fuzzy Difference-in-Difference Designs

## Syntax

```r
fuzzydid(
  data,
  formula,
  group,
  time,
  group_forward = NULL,
  did = FALSE,
  tc = FALSE,
  cic = FALSE,
  lqte = FALSE,
  newcateg = NULL,
  numerator = FALSE,
  partial = FALSE,
  nose = FALSE,
  cluster = NULL,
  breps = 50,
  eqtest = FALSE,
  modelx = NULL,
  sieves = FALSE,
  sieveorder = NULL,
  tagobs = FALSE,
  backend = c("auto", "native"),
  seed = 1
)
```

## Description

`fuzzydid()` computes estimators of local average and quantile treatment effects in fuzzy DID designs, following de Chaisemartin and D'Haultfoeuille (2018a). It also computes their standard errors and confidence intervals.

**Arguments:**

- `data`: Data frame containing all variables.
- `formula`: A formula of the form `y ~ d` or `y ~ d + x1 + x2`, where `y` is the outcome variable, `d` is the treatment variable, and `x1`, `x2`, etc. are covariates.
- `group`: Name of the group variable (backward group for multi-period designs). See Section 4.2 of [de Chaisemartin et al. (2018b)](https://sites.google.com/site/clementdechaisemartin/statapaper_fuzzydid.pdf) for details on constructing this variable.
- `time`: Name of the time period variable.
- `group_forward`: Optional name of the forward group variable for multi-period designs.

A detailed introduction to the methodology is given in [de Chaisemartin et al. (2018b)](https://sites.google.com/site/clementdechaisemartin/statapaper_fuzzydid.pdf).

## Options

**Estimators:**

- `did`: Logical; computes the Wald-DID estimator.
- `tc`: Logical; computes the Wald-TC estimator.
- `cic`: Logical; computes the Wald-CIC estimator. Only available when no covariates are included.
- `lqte`: Logical; computes estimators of the LQTE for quantiles of order 5%, 10%, ..., 95%. Only available when D, G, and T are binary, and no covariates are included.

At least one of `did`, `tc`, `cic`, or `lqte` must be specified. If several are specified, all requested estimators are computed.

**Treatment categorization:**

- `newcateg`: Numeric vector of upper bounds to group treatment values together for Wald-TC and Wald-CIC. Useful when treatment takes many values. See Section 3.3 of [de Chaisemartin et al. (2018b)](https://sites.google.com/site/clementdechaisemartin/statapaper_fuzzydid.pdf).

**Numerators and bounds:**

- `numerator`: Logical; return only the numerators of Wald-DID, Wald-TC, and Wald-CIC estimators. Useful for placebo tests (see Section 3.3.3 of the supplement of de Chaisemartin and D'Haultfoeuille 2018a).
- `partial`: Logical; compute bounds on local average treatment effects in the absence of a "stable" control group. Only available without covariates.

**Inference:**

- `nose`: Logical; compute only point estimates, not standard errors.
- `cluster`: Name of cluster variable for block bootstrap. Only one clustering variable is allowed.
- `breps`: Number of bootstrap replications. Default is 50.
- `eqtest`: Logical; perform equality tests between estimands when at least two of `did`, `tc`, `cic` are specified.

**Covariates:**

- `modelx`: Character vector specifying parametric methods for estimating conditional expectations in Wald-DID and Wald-TC with covariates. Two entries required for binary treatments; three for ordered multi-valued treatments. Values must be `"ols"`, `"logit"`, or `"probit"`.
- `sieves`: Logical; use nonparametric sieve estimation for conditional expectations.
- `sieveorder`: Optional sieve order control when `sieves = TRUE`. Default `NULL` selects order by deterministic 5-fold CV. A scalar applies to both outcome and treatment sieve bases. A length-2 vector `(outcome, treatment)` is supported for backward compatibility. Values must be ≥ 2 and satisfy the basis cap `min(4800, floor(n/5))`.

When covariates are included and neither `modelx` nor `sieves` is specified, all conditional expectations are estimated by OLS by default.

**Other:**

- `tagobs`: Logical; return a logical mask of observations used by `fuzzydid()`.
- `backend`: One of `"auto"` or `"native"`. The `"stata"` backend is no longer supported.
- `seed`: Preserved for API compatibility. Bootstrap uses a fixed Stata-parity seed (1) when `nose = FALSE`.

## Returned Values

An object of class `"fuzzydid"` containing:

**Data frames:**

- `late`: LATE estimates with columns: `estimator`, `estimate`, `std.error`, `conf.low`, `conf.high`
- `eqtest`: Equality test results (if `eqtest = TRUE`)
- `lqte`: LQTE estimates at quantiles 0.05, 0.10, ..., 0.95 (if `lqte = TRUE`)

**Matrices (Stata-parity):**

- `matrices$b_LATE`: k × 1 matrix of requested estimators
- `matrices$se_LATE`: k × 1 matrix of bootstrap standard errors
- `matrices$ci_LATE`: k × 2 matrix of 95% percentile bootstrap confidence intervals
- `matrices$b_LQTE`: 19 × 1 matrix of LQTE estimates at quantiles 0.05–0.95
- `matrices$se_LQTE`: 19 × 1 matrix of LQTE bootstrap standard errors
- `matrices$ci_LQTE`: 19 × 2 matrix of LQTE 95% confidence intervals

**Counts:**

- `n`: Number of observations used
- `n11`, `n10`, `n01`, `n00`: Cell sizes for (G,T) combinations
- `n_reps`: Number of bootstrap replications requested
- `n_misreps`: Number of failed/degenerate bootstrap replications
- `share_failures`: Proportion of failed replications

## Examples

### Generate the dataset

```r

# Generate simulated data (saved to CSV for R/Stata parity verification)
set.seed(50321)
n_cell <- 80
df <- rbind(
  data.frame(y = rnorm(n_cell, 1 + 1.8 * rbinom(n_cell, 1, 0.20)), g = 0, t = 0, d = rbinom(n_cell, 1, 0.20)),
  data.frame(y = rnorm(n_cell, 1 + 0.5 + 1.8 * rbinom(n_cell, 1, 0.35)), g = 0, t = 1, d = rbinom(n_cell, 1, 0.35)),
  data.frame(y = rnorm(n_cell, 1 + 0.7 + 1.8 * rbinom(n_cell, 1, 0.30)), g = 1, t = 0, d = rbinom(n_cell, 1, 0.30)),
  data.frame(y = rnorm(n_cell, 1 + 0.7 + 0.5 + 1.8 * rbinom(n_cell, 1, 0.70)), g = 1, t = 1, d = rbinom(n_cell, 1, 0.70))
)

# Save for Stata comparison
write.csv(df, "fuzzydid_example.csv", row.names = FALSE)
```

### R

```r
library(Rfuzzydid)
df <- read.csv("fuzzydid_example.csv")

fit <- fuzzydid(
  data = df,
  formula = y ~ d,
  group = "g",
  time = "t",
  did = TRUE,
  tc = TRUE,
  cic = TRUE,
  breps = 50
)

summary(fit)
```

### Stata

```stata
import delimited "fuzzydid_example.csv", clear

fuzzydid y g t d, did tc cic breps(50)
```

**Note:** Point estimates from R and Stata will be identical, but bootstrap confidence intervals will differ due to RNG differences between the two platforms. Results remain comparable across implementations.

## References

- de Chaisemartin, C. and D'Haultfoeuille, X. 2018a. [Fuzzy Differences-in-Differences](https://sites.google.com/site/clementdechaisemartin/fuzzy_did.pdf). *Review of Economic Studies*, 85(2): 999-1028.

- de Chaisemartin, C., D'Haultfoeuille, X., and Guyonvarch, Y. 2018b. [Fuzzy Differences-in-Differences with Stata](https://sites.google.com/site/clementdechaisemartin/statapaper_fuzzydid.pdf). *Stata Journal*.


## License

AGPL-3.0
