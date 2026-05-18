#' srr_stats_regression_diagnostics
#'
#' Regression standards addressed by runtime documentation, plotting, and
#' regression-specific tests.
#'
#' @srrstats {RE5.0} Runtime scaling is documented in the methods vignette.
#' @srrstats {RE6.0} `plot.fuzzydid()` provides a default plot method.
#' @srrstats {RE6.1} S3 plot dispatch is implemented.
#' @srrstats {RE6.2} Fitted-value plots are not meaningful for causal estimand
#' summaries; `plot()` displays estimates and intervals instead.
#' @srrstats {RE6.3} Forecast distinction is not applicable; no forecasts are
#' produced.
#' @srrstats {RE7.0} Tests cover noiseless predictor relationships.
#' @srrstats {RE7.0a} Perfectly noiseless input is not invalid for these
#' deterministic estimators.
#' @srrstats {RE7.1} Tests cover exact predictor-response relationships.
#' @srrstats {RE7.1a} Timing comparisons are not meaningful for these small
#' deterministic support checks.
#' @srrstats {RE7.2} Tests cover `tagobs`, `formula()`, and `nobs()` retention.
#' @srrstats {RE7.3} Tests cover `coef()`, `confint()`, `formula()`, `nobs()`,
#' `vcov()`, `print()`, and `plot()`.
#' @srrstats {RE7.4} Forecast intervals are not applicable.
#' @noRd
NULL
