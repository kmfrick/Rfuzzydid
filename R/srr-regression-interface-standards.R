#' srr_stats_regression_interface
#'
#' Regression standards addressed by the formula interface, preprocessing, and
#' S3 model-like methods.
#'
#' @srrstats {RE1.0} `fuzzydid()` uses a formula interface.
#' @srrstats {RE1.1} The methods vignette documents formula parsing and feature
#' construction.
#' @srrstats {RE1.2} Function docs document accepted/rejected predictor types.
#' @srrstats {RE1.3} `tagobs`, `formula()`, `nobs()`, and output row names retain
#' relevant submitted-data metadata.
#' @srrstats {RE1.3a} Documentation states that raw row names are not carried
#' into estimand tables because rows are estimands, not observations.
#' @srrstats {RE1.4} The methods vignette documents identifying assumptions and
#' implications of violations.
#' @srrstats {RE2.0} Treatment recoding, covariate expansion, and sieve
#' transformations are documented.
#' @srrstats {RE2.1} Missing/non-finite handling is documented and tested.
#' @srrstats {RE2.2} A fuzzy DID estimand uses one complete estimation sample;
#' fitting with missing responses for prediction is not applicable.
#' @srrstats {RE2.3} Centering/offsets are not part of the published fuzzy DID
#' estimators.
#' @srrstats {RE2.4} Degenerate and collinear designs are tested through
#' estimator support and nuisance-regression paths.
#' @srrstats {RE2.4a} Predictor collinearity is covered by covariate tests.
#' @srrstats {RE2.4b} Exact predictor-response relationships are covered by
#' deterministic fixture tests.
#' @srrstats {RE3.0} No iterative primary model is fit; nuisance `glm()` calls
#' are not exposed as fitted models.
#' @srrstats {RE3.1} See RE3.0.
#' @srrstats {RE3.2} See RE3.0.
#' @srrstats {RE3.3} See RE3.0.
#' @srrstats {RE4.0} `fuzzydid()` returns an S3 model-like object.
#' @srrstats {RE4.1} Empty model objects are not useful for these estimands.
#' @srrstats {RE4.2} `coef.fuzzydid()` extracts point estimates.
#' @srrstats {RE4.3} `confint.fuzzydid()` extracts stored intervals.
#' @srrstats {RE4.4} `formula.fuzzydid()` extracts the model formula.
#' @srrstats {RE4.5} `nobs.fuzzydid()` extracts sample size.
#' @srrstats {RE4.6} `vcov.fuzzydid()` extracts bootstrap covariance.
#' @srrstats {RE4.7} No primary convergence statistics exist.
#' @srrstats {RE4.8} Response variable role is retained in the formula.
#' @srrstats {RE4.9} Observation-level fitted responses are not defined for
#' causal estimand summaries.
#' @srrstats {RE4.10} Observation-level residuals are not defined for causal
#' estimand summaries.
#' @srrstats {RE4.11} `glance.fuzzydid()` returns sample and bootstrap summary
#' statistics.
#' @srrstats {RE4.12} Transformation behavior is documented in the methods
#' vignette.
#' @srrstats {RE4.13} Predictor roles are retained in the formula and options.
#' @srrstats {RE4.14} Forecast errors are not applicable; no `predict()` method
#' is provided.
#' @srrstats {RE4.15} See RE4.14.
#' @srrstats {RE4.16} See RE4.14.
#' @srrstats {RE4.17} `print.fuzzydid()` summarizes inputs and estimates.
#' @srrstats {RE4.18} `summary.fuzzydid()` prints estimator tables.
#' @noRd
NULL
