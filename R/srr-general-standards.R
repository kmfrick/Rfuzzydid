#' srr_stats_general
#'
#' `Rfuzzydid` implements fuzzy difference-in-differences estimators from
#' de Chaisemartin and D'Haultfoeuille (2018) and the Stata `fuzzydid`
#' implementation described by de Chaisemartin, D'Haultfoeuille, and
#' Guyonvarch (2019). The package is an R implementation of existing causal
#' estimators with a formula-first interface and native R computation.
#'
#' @srrstatsVerbose TRUE
#'
#' @srrstats {G1.0} Primary academic references are listed in DESCRIPTION,
#' README, vignettes, and inst/CITATION.
#' @srrstats {G1.1} README and vignettes describe the package as an R port of
#' Stata `fuzzydid`.
#' @srrstats {G1.2} README documents the lifecycle and parity status.
#' @srrstats {G1.3} The methods vignette defines DID, TC, CIC, LQTE, first-stage,
#' switchers, partial bounds, and equality-test terminology.
#' @srrstats {G1.4} Exported functions use roxygen2 documentation.
#' @srrstats {G1.4a} Internal helpers are documented through scoped comments and
#' tests because they are implementation details in one source file.
#' @srrstats {G1.5} Reproducible vignettes include code for reported examples and
#' paper-replication claims.
#' @srrstats {G1.6} Parity tests and vignettes compare against frozen Stata
#' `fuzzydid` outputs, the relevant prior implementation.
#'
#' @srrstats {G2.0} Scalar, vector-length, and option-combination assertions are
#' implemented in `fuzzydid()` and validation helpers.
#' @srrstats {G2.0a} Function docs document scalar and vector expectations.
#' @srrstats {G2.1} Type assertions are implemented for data, formula, outcome,
#' treatment, group, time, covariates, options, and cluster variables.
#' @srrstats {G2.1a} Function docs document accepted vector types.
#' @srrstats {G2.2} Matrix/list columns are rejected for univariate roles.
#' @srrstats {G2.3} Character inputs are variable names or method names.
#' @srrstats {G2.3a} `match.arg()` and explicit allowed-value checks restrict
#' character arguments.
#' @srrstats {G2.3b} Method names are documented and tested as case-sensitive.
#' @srrstats {G2.4} Type conversion is explicit in formula parsing, covariate
#' expansion, treatment recoding, and output table construction.
#' @srrstats {G2.4a} Integer conversion is explicit for logical/group/time masks.
#' @srrstats {G2.4b} Numeric conversion is explicit for outcome, treatment, and
#' estimator computations.
#' @srrstats {G2.4c} Character conversion is explicit for qualitative covariates.
#' @srrstats {G2.4d} Factor inputs are consumed as qualitative covariates.
#' @srrstats {G2.4e} Factor covariates are converted through character labels.
#' @srrstats {G2.5} Factor covariates are documented as unordered qualitative
#' predictors.
#' @srrstats {G2.6} One-dimensional inputs are checked before computation.
#' @srrstats {G2.7} `data.frame` inputs are supported and documented.
#' @srrstats {G2.8} Pre-processing normalizes data to one internal list/data-frame
#' representation.
#' @srrstats {G2.9} The package rejects lossy inputs rather than silently
#' converting them.
#' @srrstats {G2.10} Column extraction uses explicit `drop = FALSE` where needed.
#' @srrstats {G2.11} Standard vector columns are accepted and tested.
#' @srrstats {G2.12} List and matrix columns in analysis roles trigger clear
#' errors.
#' @srrstats {G2.13} Missing data are checked before analytic routines.
#' @srrstats {G2.14} Complete-case missing-data handling is documented.
#' @srrstats {G2.14a} Error-on-missing is not offered because Stata `fuzzydid`
#' parity uses complete-case estimation.
#' @srrstats {G2.14b} Missing values are ignored through documented complete-case
#' filtering; `tagobs` exposes retained rows.
#' @srrstats {G2.14c} Imputation would change the causal estimand and is out of
#' scope.
#' @srrstats {G2.15} Mean, quantile, covariance, and bootstrap helpers filter or
#' reject missing/non-finite values before use.
#' @srrstats {G2.16} `Inf` and `-Inf` are rejected; `NA`/`NaN` are complete-case
#' filtered.
#'
#' @srrstats {G3.0} Floating-point comparisons are limited to support-membership
#' checks inherited from discrete treatment/group/time design variables.
#' @srrstats {G3.1} `vcov()` returns empirical bootstrap covariance for stored
#' estimates; user-selectable covariance algorithms are not part of these
#' bootstrap estimators.
#' @srrstats {G3.1a} See G3.1.
#' @srrstats {G4.0} Package functions do not write statistical outputs to local
#' files.
#' @noRd
NULL
