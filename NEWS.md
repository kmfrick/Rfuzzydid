# Rfuzzydid 1.3.0

- Fixed bootstrap percentile confidence intervals to match Stata `fuzzydid`:
  failed replications (carrying a sentinel value) are dropped when computing the
  bootstrap standard error but retained when computing the percentile interval,
  so their extreme values widen the interval under weak support. The previous
  behavior discarded them for both, producing intervals that were too narrow
  whenever `share_failures > 0`.
- `W_CIC` is now reported as non-estimable (an informative error, rather than a
  finite estimate over the remaining strata) when a positive-mass treated
  baseline treatment stratum has no control-group support and thus an undefined
  quantile map. The previous behavior silently dropped the unsupported stratum,
  changing the estimand and reweighting the supported strata.
- Fixed multi-period aggregation of `W_DID`, `W_TC` and `W_CIC` to match Stata
  `fuzzydid`. Contributions are now weighted by the post-period switcher-cell
  count (`n11`, not `n11 + n10`) and the design sign of the switcher arm, and
  `W_TC`/`W_CIC` are weighted averages of the per-subdesign Wald ratios (anchored
  on the DID denominator) rather than ratios of summed numerators/denominators.
  Results now reproduce Stata to 1e-6 on unbalanced multi-period designs; the
  previous behavior differed on any unbalanced or multi-arm panel.
- Fixed the `sieves = TRUE` basis for two or more continuous covariates: it now
  uses a total-degree polynomial (including cross-products such as `x1 * x2`),
  matching Stata's Legendre sieve span. The previous additive-only basis omitted
  interactions and could not recover interaction effects.
- The covariate-adjusted stable-control (`special_case`) path now returns a
  missing `TC` estimate when the adjusted DID components are non-finite, instead
  of silently falling back to a treated-only before/after change.
- `numerator = TRUE` now matches Stata by erroring unless the design has exactly
  two time periods and two treatment groups (the only case in which the
  reduced-form numerator is defined).
- Added frozen Stata-parity regression tests for the unbalanced multi-period and
  two-covariate sieve designs.

# Rfuzzydid 1.2.1

- `broom` is now in `Depends`, so `library(Rfuzzydid)` attaches it
  automatically; `tidy()` and `glance()` work without a `broom::` prefix.

# Rfuzzydid 1.2

- Added SRR standards documentation and related package-review metadata.
- Improved the `plot()` method and examples for stable non-interactive checks.
- Expanded native API tests for type handling and extractor behavior.
- Refactored internal `fuzzydid()` helpers without changing the user-facing API.

# Rfuzzydid 1.1.1

- Fixed CRAN no-Suggests checks by making vignette rendering work when
  optional vignette-only packages are unavailable.
- Added rOpenSci-oriented package checks, coverage, pkgdown, and CI
  improvements.
- Added a package citation and clarified license metadata.
- Documented reproducible bootstrap seeds and added regression coverage.
- Simplified Stata parity maintenance by keeping frozen constants in the test
  suite and dropping committed generator tooling.

# Rfuzzydid 1.1.0

- Removed the vendored Stata package sources from the installed package.
- Kept Stata parity checks as frozen test goldens instead of a runtime test
  dependency.
- Added missing package metadata and user-facing contribution guidance.
- Expanded documentation for `fuzzydid()` and its S3 methods.
