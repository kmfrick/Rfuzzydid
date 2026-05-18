#' @title fuzzydid
#' @description
#' Formula-first interface for fuzzy difference-in-differences estimators.
#' Estimation is fully native in R. The object returned by \code{fuzzydid()}
#' is an estimand summary rather than a predictive regression model: it stores
#' local average and local quantile treatment-effect estimates, bootstrap
#' uncertainty summaries, design cell counts, and metadata needed by extractor
#' methods.
#' @param data A \code{data.frame}.
#' @param formula Formula of the form `y ~ d + covariates`.
#' @param treatment Optional treatment variable name for multi-term formulas.
#'   If `NULL`, treatment is inferred from formula RHS when unambiguous.
#' @param group Name of the group variable (backward group for multi-period).
#' @param time Name of the time variable.
#' @param group_forward Optional name of the forward group variable for
#'   multi-period designs.
#' @param did Logical; compute the Wald-DID estimator.
#' @param tc Logical; compute the Wald-TC estimator.
#' @param cic Logical; compute the Wald-CIC estimator.
#' @param lqte Logical; compute local quantile treatment effects.
#' @param newcateg Optional numeric vector of upper bounds used to recategorize
#'   treatment values for TC/CIC.
#' @param numerator Logical; return estimator numerators for DID/TC/CIC.
#' @param partial Logical; request TC partial-identification bounds.
#' @param nose Logical; skip bootstrap standard errors and confidence intervals.
#' @param cluster Optional name of cluster variable for one-way clustered
#'   bootstrap resampling.
#' @param breps Integer number of bootstrap replications.
#' @param eqtest Logical; compute equality tests across requested LATE estimands.
#' @param modelx Optional native covariate-adjusted methods (`ols`, `logit`,
#'   `probit`). Two entries are required for binary treatments and three for
#'   ordered multi-valued treatments.
#' @param sieves Logical; use sieve expansion for continuous covariates.
#' @param sieveorder Optional sieve order control for `sieves = TRUE`.
#'   `NULL` (default) selects order by deterministic 5-fold CV.
#'   A scalar value applies to both outcome and treatment sieve bases.
#'   A length-2 vector is accepted for backward compatibility and interpreted
#'   as `(outcome_order, treatment_order)`.
#' @param tagobs Logical; return logical mask of observations used.
#' @param backend One of `"auto"` or `"native"`.
#' @param seed Optional integer seed used for bootstrap resampling when
#'   \code{nose = FALSE}. If \code{NULL} (default), bootstrap draws use the
#'   current RNG state. Supply a value to make bootstrap standard errors,
#'   confidence intervals, and diagnostics reproducible.
#' @details
#' \code{fuzzydid()} uses complete cases across the outcome, treatment, group,
#' time, optional forward-group, covariate, and cluster variables. Missing
#' \code{NA} and \code{NaN} values are dropped; non-finite numeric values such
#' as \code{Inf} and \code{-Inf} are rejected. The outcome and treatment must
#' be numeric vectors. Group and time identifiers must be numeric vectors; with
#' one group variable, group values must be in \code{{0, 1, NA}}. Covariates may
#' be numeric, factor, character, or logical vectors. Numeric covariates enter
#' as continuous predictors; factor, character, and logical covariates enter as
#' qualitative predictors expanded to indicator columns. When \code{sieves =
#' TRUE}, continuous covariates are expanded to polynomial sieve terms.
#'
#' Standard errors and confidence intervals are percentile bootstrap summaries.
#' Use \code{seed} to make bootstrap draws reproducible. If \code{tagobs =
#' TRUE}, the returned object includes a logical vector identifying the input
#' rows retained after complete-case filtering.
#'
#' @return An object of class \code{"fuzzydid"}. This is a list whose
#'   \code{late} component is a data frame of requested LATE-type estimators
#'   with columns \code{estimator}, \code{estimate}, \code{std.error},
#'   \code{conf.low}, and \code{conf.high}. \code{eqtest} is either
#'   \code{NULL} or an analogous data frame of pairwise equality contrasts,
#'   and \code{lqte} is either \code{NULL} or a data frame with columns
#'   \code{quantile}, \code{estimate}, \code{std.error}, \code{conf.low}, and
#'   \code{conf.high} for local quantile treatment effects. Additional
#'   components include \code{matrices}, a named list of Stata-style result
#'   matrices; \code{tagobs}, an optional logical mask of retained
#'   observations; sample-size diagnostics \code{n}, \code{n11}, \code{n10},
#'   \code{n01}, and \code{n00}; bootstrap diagnostics \code{n_reps},
#'   \code{n_misreps}, and \code{share_failures}; and metadata such as
#'   \code{backend}, \code{call}, and \code{options}. The estimate tables
#'   report point estimates and, unless \code{nose = TRUE}, bootstrap
#'   standard errors and percentile confidence limits.
#' @examples
#' make_example_cell <- function(g, t, ones, n_cell = 20L) {
#'   data.frame(
#'     g = rep.int(g, n_cell),
#'     t = rep.int(t, n_cell),
#'     d = c(rep.int(1L, ones), rep.int(0L, n_cell - ones))
#'   )
#' }
#'
#' df <- rbind(
#'   make_example_cell(0L, 0L, 4L),
#'   make_example_cell(0L, 1L, 8L),
#'   make_example_cell(1L, 0L, 6L),
#'   make_example_cell(1L, 1L, 16L)
#' )
#' df$id <- seq_len(nrow(df))
#' df$y <- 1 + 0.5 * df$g + 0.4 * df$t + 2 * df$d + sin(df$id / 7)
#' example_data <- df
#'
#' fit <- fuzzydid(
#'   data = example_data,
#'   formula = y ~ d,
#'   treatment = NULL,
#'   group = "g",
#'   time = "t",
#'   group_forward = NULL,
#'   did = TRUE,
#'   tc = TRUE,
#'   cic = TRUE,
#'   lqte = TRUE,
#'   newcateg = c(0, 1),
#'   cluster = NULL,
#'   modelx = NULL,
#'   sieveorder = NULL,
#'   seed = NULL,
#'   nose = TRUE
#' )
#'
#' fit$late
#' @export
fuzzydid <- function(
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
  seed = NULL,
  treatment = NULL
) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula like y ~ d + x1 + x2.", call. = FALSE)
  }

  backend <- match.arg(backend)
  logical_opts <- list(
    did = did,
    tc = tc,
    cic = cic,
    lqte = lqte,
    numerator = numerator,
    partial = partial,
    nose = nose,
    eqtest = eqtest,
    sieves = sieves,
    tagobs = tagobs
  )
  valid_logical <- vapply(
    logical_opts,
    function(x) is.logical(x) && length(x) == 1L && !is.na(x),
    logical(1)
  )
  if (!all(valid_logical)) {
    bad <- names(valid_logical)[!valid_logical][[1L]]
    stop(sprintf("`%s` must be TRUE or FALSE.", bad), call. = FALSE)
  }

  opts <- list(
    did = did,
    tc = tc,
    cic = cic,
    lqte = lqte,
    newcateg = newcateg,
    numerator = numerator,
    partial = partial,
    nose = nose,
    cluster = cluster,
    breps = breps,
    eqtest = eqtest,
    modelx = modelx,
    sieves = sieves,
    sieveorder = sieveorder,
    tagobs = tagobs,
    seed = seed
  )

  .validate_common_options(opts)
  parsed <- .parse_fuzzy_formula(formula, data, treatment = treatment)
  prepared <- .prepare_input_data(
    data = data,
    y_name = parsed$y_name,
    d_name = parsed$d_name,
    covariates = parsed$covariates,
    group = group,
    time = time,
    group_forward = group_forward,
    cluster = cluster
  )

  out <- .run_native_backend(
    prepared = prepared,
    opts = opts
  )

  out$backend <- "native"
  out$formula <- formula
  out$group <- group
  out$time <- time
  out$group_forward <- group_forward
  out$call <- match.call()
  out$options <- opts
  class(out) <- "fuzzydid"
  out
}

.validate_common_options <- function(opts) {
  if (!opts$did && !opts$tc && !opts$cic && !opts$lqte) {
    stop(
      "At least one estimator must be requested: did, tc, cic, or lqte.",
      call. = FALSE
    )
  }

  if (opts$nose && !is.null(opts$cluster)) {
    stop("`nose = TRUE` cannot be used with `cluster`.", call. = FALSE)
  }

  if (opts$eqtest && opts$numerator) {
    stop("`eqtest` cannot be used with `numerator`.", call. = FALSE)
  }

  if (opts$partial) {
    if (!opts$tc) {
      stop("`partial` can only be used with `tc = TRUE`.", call. = FALSE)
    }
    if (opts$numerator) {
      stop("`partial` cannot be combined with `numerator`.", call. = FALSE)
    }
  }

  if (!is.null(opts$newcateg)) {
    if (!is.numeric(opts$newcateg)) {
      stop("`newcateg` must be numeric.", call. = FALSE)
    }
    if (length(opts$newcateg) < 2L) {
      stop("`newcateg` must contain at least two values.", call. = FALSE)
    }
    if (is.unsorted(opts$newcateg, strictly = TRUE)) {
      stop("`newcateg` values must be strictly increasing.", call. = FALSE)
    }
  }

  if (!is.null(opts$modelx)) {
    if (!is.character(opts$modelx)) {
      stop("`modelx` must be a character vector.", call. = FALSE)
    }
    if (!(length(opts$modelx) %in% c(2L, 3L))) {
      stop("`modelx` must contain two or three entries.", call. = FALSE)
    }
    allowed <- c("ols", "logit", "probit")
    if (!all(opts$modelx %in% allowed)) {
      stop("`modelx` entries must be one of: ols, logit, probit.", call. = FALSE)
    }
  }

  if (!is.null(opts$sieveorder)) {
    if (!is.numeric(opts$sieveorder) || !(length(opts$sieveorder) %in% c(1L, 2L))) {
      stop("`sieveorder` must be numeric with length 1 or 2.", call. = FALSE)
    }
    if (!all(is.finite(opts$sieveorder))) {
      stop("`sieveorder` values must be finite.", call. = FALSE)
    }
    if (any(opts$sieveorder < 2)) {
      stop("`sieveorder` values must be >= 2.", call. = FALSE)
    }
    if (length(opts$sieveorder) == 2L) {
      warning(
        "Using length-2 `sieveorder` is deprecated; prefer a scalar order or `NULL` for CV.",
        call. = FALSE
      )
    }
  }

  if (!opts$sieves && !is.null(opts$sieveorder)) {
    stop("`sieveorder` cannot be used without `sieves = TRUE`.", call. = FALSE)
  }

  if (!is.null(opts$modelx) && (opts$sieves || !is.null(opts$sieveorder))) {
    stop("`modelx` cannot be used in combination with `sieves`/`sieveorder`.", call. = FALSE)
  }

  if (!is.null(opts$cluster)) {
    if (!is.character(opts$cluster) || length(opts$cluster) != 1L) {
      stop("`cluster` must be a single variable name.", call. = FALSE)
    }
  }

  if (!is.null(opts$seed)) {
    seed <- opts$seed
    is_integer_like <- is.numeric(seed) &&
      length(seed) == 1L &&
      !is.na(seed) &&
      is.finite(seed) &&
      (seed == as.integer(seed)) &&
      seed >= 0
    if (!is_integer_like) {
      stop("`seed` must be NULL or a non-negative integer scalar.", call. = FALSE)
    }
  }

  if (!opts$nose) {
    breps <- opts$breps
    is_integer_like <- is.numeric(breps) &&
      length(breps) == 1L &&
      !is.na(breps) &&
      is.finite(breps) &&
      (breps == as.integer(breps))
    if (!is_integer_like || breps < 2) {
      stop("`breps` must be an integer scalar >= 2.", call. = FALSE)
    }
  }
}

.parse_fuzzy_formula <- function(formula, data, treatment = NULL) {
  tr <- stats::terms(formula, data = data)
  rhs_terms <- attr(tr, "term.labels")
  if (length(rhs_terms) < 1L) {
    stop("Formula RHS must include treatment variable `d`.", call. = FALSE)
  }

  lhs_vars <- all.vars(formula[[2L]])
  if (length(lhs_vars) != 1L) {
    stop("Formula LHS must contain exactly one outcome variable.", call. = FALSE)
  }

  if (!is.null(treatment)) {
    if (!is.character(treatment) || length(treatment) != 1L || !nzchar(treatment)) {
      stop("`treatment` must be a single variable name.", call. = FALSE)
    }
    if (!(treatment %in% rhs_terms)) {
      stop("`treatment` must appear on the formula RHS.", call. = FALSE)
    }
    d_name <- treatment
  } else if (length(rhs_terms) == 1L) {
    d_name <- rhs_terms[[1L]]
  } else {
    d_matches <- which(rhs_terms == "d")
    if (length(d_matches) == 1L) {
      d_name <- "d"
    } else {
      stop(
        "Unable to infer treatment from formula RHS; provide `treatment = \"<var>\"`.",
        call. = FALSE
      )
    }
  }

  y_name <- lhs_vars[[1L]]
  covariates <- rhs_terms[rhs_terms != d_name]

  needed <- c(y_name, d_name, covariates)
  missing_vars <- setdiff(needed, names(data))
  if (length(missing_vars) > 0L) {
    stop(
      sprintf("Variables not found in `data`: %s", toString(missing_vars)),
      call. = FALSE
    )
  }

  list(
    y_name = y_name,
    d_name = d_name,
    covariates = covariates
  )
}

.infer_covariate_types <- function(data, covariates) {
  if (length(covariates) == 0L) {
    return(list(continuous = character(0), qualitative = character(0)))
  }

  is_qual <- vapply(
    covariates,
    function(v) {
      x <- data[[v]]
      is.factor(x) || is.character(x) || is.logical(x)
    },
    logical(1)
  )

  list(
    continuous = covariates[!is_qual],
    qualitative = covariates[is_qual]
  )
}

.is_plain_vector <- function(x) {
  is.atomic(x) && is.null(dim(x))
}

.validate_input_columns <- function(data, y_name, d_name, covariates, group, time, group_forward, cluster) {
  required_numeric <- unique(c(y_name, d_name, group, time, group_forward))
  required_numeric <- required_numeric[!is.na(required_numeric) & nzchar(required_numeric)]

  for (nm in required_numeric) {
    x <- data[[nm]]
    if (!.is_plain_vector(x) || !is.numeric(x)) {
      stop(sprintf("`%s` must be a numeric vector.", nm), call. = FALSE)
    }
    bad <- !is.na(x) & !is.finite(x)
    if (any(bad)) {
      stop(sprintf("`%s` must not contain Inf or -Inf.", nm), call. = FALSE)
    }
  }

  for (nm in covariates) {
    x <- data[[nm]]
    ok_type <- .is_plain_vector(x) &&
      (is.numeric(x) || is.factor(x) || is.character(x) || is.logical(x))
    if (!ok_type) {
      stop(
        sprintf("Covariate `%s` must be a numeric, factor, character, or logical vector.", nm),
        call. = FALSE
      )
    }
    if (is.numeric(x)) {
      bad <- !is.na(x) & !is.finite(x)
      if (any(bad)) {
        stop(sprintf("Covariate `%s` must not contain Inf or -Inf.", nm), call. = FALSE)
      }
    }
  }

  if (!is.null(cluster)) {
    x <- data[[cluster]]
    if (!.is_plain_vector(x)) {
      stop(sprintf("Cluster variable `%s` must be a vector.", cluster), call. = FALSE)
    }
  }

  invisible(NULL)
}

.prepare_input_data <- function(data, y_name, d_name, covariates, group, time, group_forward, cluster) {
  role_vars <- unique(c(y_name, d_name, group, time, group_forward, covariates, cluster))
  role_vars <- role_vars[!is.na(role_vars) & nzchar(role_vars)]

  missing_vars <- setdiff(role_vars, names(data))
  if (length(missing_vars) > 0L) {
    stop(
      sprintf("Required variable(s) not found in `data`: %s", toString(missing_vars)),
      call. = FALSE
    )
  }

  .validate_input_columns(
    data = data,
    y_name = y_name,
    d_name = d_name,
    covariates = covariates,
    group = group,
    time = time,
    group_forward = group_forward,
    cluster = cluster
  )

  cov_types <- .infer_covariate_types(data, covariates)

  df <- data
  if (is.null(group_forward)) {
    complete_mask <- stats::complete.cases(df[, role_vars, drop = FALSE])
  } else {
    other_vars <- unique(c(y_name, d_name, time, covariates, cluster))
    other_vars <- other_vars[!is.na(other_vars) & nzchar(other_vars)]
    complete_other <- stats::complete.cases(df[, other_vars, drop = FALSE])
    valid_group <- !(is.na(df[[group]]) & is.na(df[[group_forward]]))
    complete_mask <- complete_other & valid_group
  }

  df_used <- df[complete_mask, , drop = FALSE]

  if (nrow(df_used) == 0L) {
    stop("No observations remain after dropping missing values.", call. = FALSE)
  }

  list(
    data = df,
    used = df_used,
    mask = complete_mask,
    y_name = y_name,
    d_name = d_name,
    covariates = covariates,
    covariate_types = cov_types,
    group = group,
    time = time,
    group_forward = group_forward,
    cluster = cluster
  )
}

.to_column_matrix <- function(values, rownames_vec) {
  if (length(values) == 0L) {
    return(NULL)
  }
  mat <- matrix(values, ncol = 1)
  rownames(mat) <- rownames_vec
  mat
}

.to_ci_matrix <- function(low, high, rownames_vec) {
  if (length(low) == 0L) {
    return(NULL)
  }
  mat <- cbind(low, high)
  colnames(mat) <- c("lower", "upper")
  rownames(mat) <- rownames_vec
  mat
}

.safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

.safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  stats::quantile(x, probs = p, names = FALSE, type = 1)
}

.recategorize_treatment <- function(d, newcateg) {
  d <- as.numeric(d)
  if (is.null(newcateg)) return(d)

  min_d <- min(d, na.rm = TRUE)
  max_d <- max(d, na.rm = TRUE)

  if (newcateg[[1L]] < min_d || newcateg[[1L]] >= max_d) {
    stop("First `newcateg` value must be >= min(D) and < max(D).", call. = FALSE)
  }
  if (newcateg[[length(newcateg)]] < max_d) {
    stop("Last `newcateg` value must be >= max(D).", call. = FALSE)
  }

  out <- rep.int(NA_real_, length(d))
  for (i in seq_along(d)) {
    out[[i]] <- which(d[[i]] <= newcateg)[[1L]]
  }
  as.numeric(out)
}

.counterfactual_quantile_map <- function(y_target, y_d00, y_d01) {
  if (length(y_d00) == 0L || length(y_d01) == 0L || length(y_target) == 0L) {
    return(rep.int(NA_real_, length(y_target)))
  }
  Fy <- stats::ecdf(y_d00)
  probs <- Fy(y_target)
  probs[probs < 0] <- 0
  probs[probs > 1] <- 1
  stats::quantile(y_d01, probs = probs, names = FALSE, type = 1)
}

.cdf_at <- function(sample, x) {
  if (length(sample) == 0L) {
    return(rep(NA_real_, length(x)))
  }
  vapply(x, function(xx) mean(sample <= xx), numeric(1))
}

.counterfactual_inverse_map <- function(y_target, y_from, y_to) {
  if (length(y_from) == 0L || length(y_to) == 0L || length(y_target) == 0L) {
    return(rep.int(NA_real_, length(y_target)))
  }

  probs <- .cdf_at(y_from, y_target)
  out <- stats::quantile(y_to, probs = probs, names = FALSE, type = 1)
  lower <- min(y_to, na.rm = TRUE) - max(1, diff(range(y_to, na.rm = TRUE)))
  out[probs <= 0] <- lower
  as.numeric(out)
}

.rearrange_cdf <- function(x) {
  sort(pmin(pmax(x, 0), 1), na.last = TRUE)
}

.is_binary <- function(x) {
  ux <- sort(unique(x[!is.na(x)]))
  length(ux) == 2L
}

.validate_group_values <- function(df, group, group_forward) {
  if (is.null(group_forward)) {
    bad <- !(df[[group]] %in% c(0, 1) | is.na(df[[group]]))
    if (any(bad)) {
      stop(
        "When only one group variable is provided, it must take values in {0, 1, NA}.",
        call. = FALSE
      )
    }
    if (!any(df[[group]] == 0, na.rm = TRUE)) {
      stop("The group variable must take value 0 for some observations.", call. = FALSE)
    }
    return(invisible(NULL))
  }

  for (nm in c(group, group_forward)) {
    bad <- !(df[[nm]] %in% c(-3, -1, 0, 1) | is.na(df[[nm]]))
    if (any(bad)) {
      stop(sprintf("The group variable `%s` takes values outside {-3,-1,0,1,NA}.", nm), call. = FALSE)
    }
    if (!any(df[[nm]] == 0, na.rm = TRUE)) {
      stop("The group variables must take value 0 for some observations.", call. = FALSE)
    }
  }

  invisible(NULL)
}

.validate_native_request <- function(df, prepared, opts) {
  time_vals <- sort(unique(df[[prepared$time]]))
  if (length(time_vals) < 2L) {
    stop("The fuzzydid package requires at least two time periods.", call. = FALSE)
  }

  if (is.null(prepared$group_forward) && length(time_vals) > 2L) {
    stop(
      "With more than two time periods, both backward and forward group identifiers must be used.",
      call. = FALSE
    )
  }

  .validate_group_values(df, prepared$group, prepared$group_forward)

  d_vals <- as.numeric(df[[prepared$d_name]])
  ncateg <- length(unique(d_vals[!is.na(d_vals)]))
  if (ncateg < 2L) {
    stop("The treatment variable must take at least two different values.", call. = FALSE)
  }

  has_cov <- length(prepared$covariates) > 0L
  has_cont <- length(prepared$covariate_types$continuous) > 0L

  if (has_cov && (opts$cic || opts$lqte)) {
    stop("continuous and qualitative cannot be used in combination with cic or lqte.", call. = FALSE)
  }

  if (!has_cov && !is.null(opts$modelx)) {
    stop("`modelx` requires covariates in the formula.", call. = FALSE)
  }

  if (opts$sieves && has_cov && !has_cont) {
    stop("sieves requires continuous covariates.", call. = FALSE)
  }

  if (opts$lqte) {
    if (ncateg != 2L) {
      stop("Local quantile treatment effects can only be estimated with a binary treatment variable.", call. = FALSE)
    }
    if (length(time_vals) > 2L) {
      stop("lqte cannot be used when there are more than two periods.", call. = FALSE)
    }
  }

  if (opts$partial) {
    if (has_cov) {
      stop("partial can only be used without covariates.", call. = FALSE)
    }
    if (length(time_vals) > 2L) {
      stop("partial can only be used with two time periods.", call. = FALSE)
    }
  }

  if (!is.null(opts$modelx)) {
    if (length(opts$modelx) == 2L && ncateg != 2L) {
      stop(
        "As the treatment variable takes more than two values, `modelx` must have three methods.",
        call. = FALSE
      )
    }
  }

  invisible(NULL)
}

.build_pair_frames <- function(df, group, time, group_forward) {
  t_vals <- sort(unique(df[[time]]))
  out <- list()

  if (length(t_vals) == 2L && is.null(group_forward)) {
    sub <- df
    sub$.g_star <- as.numeric(sub[[group]])
    sub$.pair_id <- sprintf("%s_%s", t_vals[[1L]], t_vals[[2L]])
    sub$.t_prev <- t_vals[[1L]]
    sub$.t_curr <- t_vals[[2L]]
    out[[1L]] <- sub[!is.na(sub$.g_star), , drop = FALSE]
    return(out)
  }

  for (i in seq.int(2L, length(t_vals))) {
    t_prev <- t_vals[[i - 1L]]
    t_curr <- t_vals[[i]]
    sub <- df[df[[time]] %in% c(t_prev, t_curr), , drop = FALSE]

    if (is.null(group_forward)) {
      g_star <- sub[[group]]
    } else {
      g_star <- ifelse(sub[[time]] == t_curr, sub[[group]], sub[[group_forward]])
    }

    sub$.g_star <- as.numeric(g_star)
    sub$.pair_id <- sprintf("%s_%s", t_prev, t_curr)
    sub$.t_prev <- t_prev
    sub$.t_curr <- t_curr
    sub <- sub[!is.na(sub$.g_star), , drop = FALSE]

    out[[length(out) + 1L]] <- sub
  }

  out
}

.extract_pair_subdesigns <- function(pair_df, time_var) {
  signs <- sort(unique(pair_df$.g_star[pair_df$.g_star != 0 & pair_df$.g_star != -3]))
  if (length(signs) == 0L) {
    return(list())
  }

  t_vals <- sort(unique(pair_df[[time_var]]))
  if (length(t_vals) != 2L) {
    return(list())
  }
  t0 <- t_vals[[1L]]
  t1 <- t_vals[[2L]]

  out <- list()
  for (sg in signs) {
    sub <- pair_df[pair_df$.g_star %in% c(0, sg), , drop = FALSE]
    sub$.g_binary <- as.integer(sub$.g_star == sg)
    sub$.t_binary <- as.integer(sub[[time_var]] == t1)

    n00 <- sum(sub$.g_binary == 0 & sub$.t_binary == 0)
    n01 <- sum(sub$.g_binary == 0 & sub$.t_binary == 1)
    n10 <- sum(sub$.g_binary == 1 & sub$.t_binary == 0)
    n11 <- sum(sub$.g_binary == 1 & sub$.t_binary == 1)

    if (any(c(n00, n01, n10, n11) == 0L)) {
      next
    }

    out[[length(out) + 1L]] <- list(
      data = sub,
      sign = sg,
      pair_id = unique(sub$.pair_id)[[1L]],
      counts = c(n11 = n11, n10 = n10, n01 = n01, n00 = n00)
    )
  }

  out
}

.build_feature_frame <- function(df, continuous, qualitative, order = 1L) {
  out <- data.frame(.row_id = seq_len(nrow(df)))
  order <- as.integer(order)
  if (is.na(order) || order < 1L) {
    order <- 1L
  }

  if (length(continuous) > 0L) {
    for (v in continuous) {
      x <- as.numeric(df[[v]])
      for (p in seq_len(order)) {
        out[[sprintf("%s_p%s", v, p)]] <- x^p
      }
    }
  }

  if (length(qualitative) > 0L) {
    for (v in qualitative) {
      x <- as.character(df[[v]])
      lvls <- sort(unique(x[!is.na(x)]))
      if (length(lvls) <= 1L) {
        next
      }
      for (lv in lvls[-1L]) {
        out[[sprintf("%s__%s", v, make.names(lv))]] <- as.numeric(x == lv)
      }
    }
  }

  out$.row_id <- NULL
  out
}

.sieve_basis_cap <- function(n_obs) {
  as.integer(min(4800L, floor(n_obs / 5)))
}

.max_sieve_order <- function(df, cov_types) {
  n_obs <- nrow(df)
  n_cont <- length(cov_types$continuous)
  if (n_cont == 0L) {
    return(NA_integer_)
  }

  cap <- .sieve_basis_cap(n_obs)
  n_qual_terms <- ncol(.build_feature_frame(df, character(0), cov_types$qualitative, order = 1L))
  room <- cap - n_qual_terms
  if (room < n_cont) {
    return(NA_integer_)
  }
  as.integer(floor(room / n_cont))
}

.ensure_sieve_basis_within_cap <- function(df, cov_types, y_order, d_order) {
  cap <- .sieve_basis_cap(nrow(df))
  if (cap < 1L) {
    stop("Sieve estimation requires at least 5 observations in each estimation sub-design.", call. = FALSE)
  }

  y_terms <- ncol(.build_feature_frame(df, cov_types$continuous, cov_types$qualitative, order = y_order))
  d_terms <- ncol(.build_feature_frame(df, cov_types$continuous, cov_types$qualitative, order = d_order))
  if (y_terms > cap || d_terms > cap) {
    stop(
      sprintf(
        "Sieve basis has too many terms (%s for outcome, %s for treatment); cap is %s = min(4800, floor(n/5)).",
        y_terms,
        d_terms,
        cap
      ),
      call. = FALSE
    )
  }
}

.cv_score_sieve_order <- function(df, y, d, cov_types, order, k_folds = 5L) {
  x <- .build_feature_frame(df, cov_types$continuous, cov_types$qualitative, order = order)
  n <- nrow(df)
  k <- min(as.integer(k_folds), n)
  if (k < 2L) {
    return(Inf)
  }

  fold_id <- ((seq_len(n) - 1L) %% k) + 1L
  fold_losses <- numeric(k)

  for (fold in seq_len(k)) {
    test_idx <- fold_id == fold
    train_idx <- !test_idx

    y_pred <- .fit_predict_vector(
      train_y = y[train_idx],
      train_x = x[train_idx, , drop = FALSE],
      pred_x = x[test_idx, , drop = FALSE],
      method = "ols"
    )
    d_pred <- .fit_predict_vector(
      train_y = d[train_idx],
      train_x = x[train_idx, , drop = FALSE],
      pred_x = x[test_idx, , drop = FALSE],
      method = "ols"
    )

    y_true <- y[test_idx]
    d_true <- d[test_idx]
    y_ok <- is.finite(y_true) & is.finite(y_pred)
    d_ok <- is.finite(d_true) & is.finite(d_pred)
    if (!any(y_ok) || !any(d_ok)) {
      return(Inf)
    }

    mse_y <- mean((y_true[y_ok] - y_pred[y_ok])^2)
    mse_d <- mean((d_true[d_ok] - d_pred[d_ok])^2)
    if (!is.finite(mse_y) || !is.finite(mse_d)) {
      return(Inf)
    }
    fold_losses[[fold]] <- mse_y + mse_d
  }

  mean(fold_losses)
}

.resolve_sieve_orders <- function(df, y, d, cov_types, opts) {
  if (!opts$sieves) {
    return(list(y_order = 1L, d_order = 1L, selected_by_cv = FALSE))
  }

  if (!is.null(opts$sieveorder)) {
    so <- as.integer(opts$sieveorder)
    if (length(so) == 1L) {
      y_order <- so[[1L]]
      d_order <- so[[1L]]
    } else {
      y_order <- so[[1L]]
      d_order <- so[[2L]]
    }
    .ensure_sieve_basis_within_cap(df, cov_types, y_order = y_order, d_order = d_order)
    return(list(y_order = y_order, d_order = d_order, selected_by_cv = FALSE))
  }

  max_order <- .max_sieve_order(df, cov_types)
  if (!is.finite(max_order) || max_order < 2L) {
    stop(
      "Unable to choose sieve order by CV: basis cap min(4800, floor(n/5)) leaves no feasible order >= 2.",
      call. = FALSE
    )
  }

  candidates <- seq.int(2L, max_order)
  scores <- vapply(
    candidates,
    function(ord) .cv_score_sieve_order(df = df, y = y, d = d, cov_types = cov_types, order = ord),
    numeric(1)
  )
  if (!any(is.finite(scores))) {
    stop("Unable to choose sieve order by CV due to non-finite fold losses.", call. = FALSE)
  }

  best <- candidates[[which.min(scores)]]
  .ensure_sieve_basis_within_cap(df, cov_types, y_order = best, d_order = best)
  list(y_order = best, d_order = best, selected_by_cv = TRUE)
}

.resolve_model_methods <- function(modelx, y_vals, d_vals) {
  if (is.null(modelx)) {
    return(list(y = "ols", d1 = "ols", d2 = "ols"))
  }

  if (length(modelx) == 2L) {
    methods <- list(y = modelx[[1L]], d1 = modelx[[2L]], d2 = modelx[[2L]])
  } else if (length(modelx) == 3L) {
    methods <- list(y = modelx[[1L]], d1 = modelx[[2L]], d2 = modelx[[3L]])
  } else {
    stop("`modelx` must contain two or three methods.", call. = FALSE)
  }

  if (!.is_binary(y_vals) && methods$y != "ols") {
    stop("The outcome variable is not binary, so E(Y|X) cannot use probit/logit.", call. = FALSE)
  }

  if (!.is_binary(d_vals) && methods$d1 != "ols") {
    stop("The treatment variable is not binary, so E(D|X) cannot use probit/logit.", call. = FALSE)
  }

  methods
}

.fit_predict_vector <- function(train_y, train_x, pred_x, method) {
  n_pred <- nrow(pred_x)
  if (n_pred == 0L) {
    return(numeric(0))
  }

  out <- rep(NA_real_, n_pred)

  if (length(train_y) == 0L) {
    return(out)
  }

  train_ok <- is.finite(train_y)
  if (ncol(train_x) > 0L) {
    train_ok <- train_ok & stats::complete.cases(train_x)
  }

  train_y <- train_y[train_ok]
  train_x <- train_x[train_ok, , drop = FALSE]

  if (length(train_y) == 0L) {
    return(out)
  }

  pred_ok <- if (ncol(pred_x) > 0L) stats::complete.cases(pred_x) else rep(TRUE, n_pred)
  pred_x_use <- pred_x[pred_ok, , drop = FALSE]

  if (length(unique(train_y)) == 1L) {
    out[pred_ok] <- unique(train_y)[[1L]]
    return(out)
  }

  train_df <- data.frame(y = train_y, train_x, check.names = FALSE)
  pred_df <- if (ncol(pred_x_use) == 0L) {
    data.frame(row.names = seq_len(nrow(pred_x_use)))
  } else {
    data.frame(pred_x_use, check.names = FALSE)
  }

  fit <- suppressWarnings(tryCatch(
    {
      if (method == "ols") {
        stats::lm(y ~ ., data = train_df)
      } else if (method == "logit") {
        stats::glm(y ~ ., data = train_df, family = stats::binomial(link = "logit"))
      } else if (method == "probit") {
        stats::glm(y ~ ., data = train_df, family = stats::binomial(link = "probit"))
      } else {
        NULL
      }
    },
    error = function(...) NULL
  ))

  if (is.null(fit)) {
    return(out)
  }

  pred_vals <- suppressWarnings(tryCatch(
    stats::predict(fit, newdata = pred_df, type = "response"),
    error = function(...) rep(NA_real_, nrow(pred_df))
  ))

  pred_vals <- as.numeric(pred_vals)
  if (method != "ols") {
    pred_vals <- pmin(pmax(pred_vals, 0), 1)
  }

  out[pred_ok] <- pred_vals
  out
}

.fit_predict_mean <- function(train_y, train_x, pred_x, method) {
  .safe_mean(.fit_predict_vector(train_y, train_x, pred_x, method))
}

.estimate_pair_no_cov <- function(sub_df, y_name, d_true_name, d_tc_name, opts) {
  y <- as.numeric(sub_df[[y_name]])
  d_true <- as.numeric(sub_df[[d_true_name]])
  d_tc <- as.numeric(sub_df[[d_tc_name]])
  g <- as.integer(sub_df$.g_binary)
  t <- as.integer(sub_df$.t_binary)

  m_y11 <- .safe_mean(y[g == 1 & t == 1])
  m_y10 <- .safe_mean(y[g == 1 & t == 0])
  m_y01 <- .safe_mean(y[g == 0 & t == 1])
  m_y00 <- .safe_mean(y[g == 0 & t == 0])

  m_d11 <- .safe_mean(d_true[g == 1 & t == 1])
  m_d10 <- .safe_mean(d_true[g == 1 & t == 0])
  m_d01 <- .safe_mean(d_true[g == 0 & t == 1])
  m_d00 <- .safe_mean(d_true[g == 0 & t == 0])

  did_num <- m_y11 - m_y10 - (m_y01 - m_y00)
  did_den <- m_d11 - m_d10 - (m_d01 - m_d00)

  out <- list(
    did_num = did_num,
    did_den = did_den,
    tc_num = NA_real_,
    tc_den = NA_real_,
    tc_inf = NA_real_,
    tc_sup = NA_real_,
    cic_num = NA_real_,
    cic_den = NA_real_,
    special_case = FALSE
  )

  special_case <- length(unique(d_tc[g == 0])) == 1L
  out$special_case <- special_case

  if (opts$tc || opts$partial) {
    if (special_case) {
      out$tc_num <- did_num
      out$tc_den <- did_den
    } else {
      levels_10 <- sort(unique(d_tc[g == 1 & t == 0]))
      tc_den <- m_d11 - m_d10

      if (opts$partial) {
        min_y <- min(y, na.rm = TRUE)
        max_y <- max(y, na.rm = TRUE)
        levels_all <- sort(unique(d_tc))
        build_inf <- 0
        build_sup <- 0

        .partial_pctile <- function(vals, pct) {
          n <- length(vals)
          if (n == 0L || !is.finite(pct)) return(NA_real_)
          if (pct <= (100 / n)) return(min(vals))
          if (pct >= 100) return(max(vals))
          as.numeric(stats::quantile(vals, probs = pct / 100, type = 2, names = FALSE))
        }

        if (length(levels_all) == 2L) {
          lv1 <- levels_all[[1L]]
          lv2 <- levels_all[[2L]]

          p_low_01 <- mean(d_tc[g == 0 & t == 1] == lv1)
          p_low_00 <- mean(d_tc[g == 0 & t == 0] == lv1)
          p_high_01 <- 1 - p_low_01
          p_high_00 <- 1 - p_low_00

          for (lv in levels_10) {
            p_d10 <- mean(d_tc[g == 1 & t == 0] == lv)
            y01_d <- y[g == 0 & t == 1 & d_tc == lv]
            y00_d <- y[g == 0 & t == 0 & d_tc == lv]
            mean_y00_d <- .safe_mean(y00_d)

            delta_inf <- NA_real_
            delta_sup <- NA_real_

            if (lv == lv1) {
              ratio <- p_low_01 / p_low_00
              y01_mean <- .safe_mean(y01_d)
              delta_inf <- (1 - ratio) * min_y + ratio * y01_mean
              delta_sup <- (1 - ratio) * max_y + ratio * y01_mean
            } else if (lv == lv2) {
              ratio <- p_high_01 / p_high_00
              qsup <- .partial_pctile(y01_d, 100 * (1 - 1 / ratio))
              qinf <- .partial_pctile(y01_d, 100 / ratio)
              if (is.finite(qinf)) {
                delta_inf <- .safe_mean(y01_d[y01_d <= qinf])
              }
              if (is.finite(qsup)) {
                delta_sup <- .safe_mean(y01_d[y01_d >= qsup])
              }
            }

            delta_inf <- delta_inf - mean_y00_d
            delta_sup <- delta_sup - mean_y00_d

            build_inf <- build_inf + p_d10 * delta_sup
            build_sup <- build_sup + p_d10 * delta_inf
          }
        } else {
          for (lv in levels_10) {
            p_d10 <- mean(d_tc[g == 1 & t == 0] == lv)
            y01_d <- y[g == 0 & t == 1 & d_tc == lv]
            y00_d <- y[g == 0 & t == 0 & d_tc == lv]

            if (length(y01_d) == 0L || length(y00_d) == 0L) {
              delta_low <- min_y - max_y
              delta_high <- max_y - min_y
            } else {
              delta <- mean(y01_d) - mean(y00_d)
              delta_low <- delta
              delta_high <- delta
            }

            build_inf <- build_inf + p_d10 * delta_high
            build_sup <- build_sup + p_d10 * delta_low
          }
        }

        num_base <- m_y11 - m_y10
        tc_inf <- (num_base - build_inf) / tc_den
        tc_sup <- (num_base - build_sup) / tc_den
        bounds <- sort(c(tc_inf, tc_sup))
        out$tc_inf <- bounds[[1L]]
        out$tc_sup <- bounds[[2L]]
        out$tc_den <- tc_den
      } else {
        build <- 0
        valid <- TRUE

        for (lv in levels_10) {
          p_d10 <- mean(d_tc[g == 1 & t == 0] == lv)
          y01_d <- .safe_mean(y[g == 0 & t == 1 & d_tc == lv])
          y00_d <- .safe_mean(y[g == 0 & t == 0 & d_tc == lv])
          if (is.na(y01_d) || is.na(y00_d)) {
            valid <- FALSE
            break
          }
          build <- build + p_d10 * (y01_d - y00_d)
        }

        out$tc_num <- if (valid) m_y11 - m_y10 - build else NA_real_
        out$tc_den <- tc_den
      }
    }
  }

  if (opts$cic) {
    levels_10 <- sort(unique(d_tc[g == 1 & t == 0]))
    q_vals <- numeric(0)
    ctrl_level <- if (special_case) unique(d_tc[g == 0])[1L] else NA_real_

    for (lv in levels_10) {
      y10 <- y[g == 1 & t == 0 & d_tc == lv]
      if (special_case) {
        y00 <- y[g == 0 & t == 0 & d_tc == ctrl_level]
        y01 <- y[g == 0 & t == 1 & d_tc == ctrl_level]
      } else {
        y00 <- y[g == 0 & t == 0 & d_tc == lv]
        y01 <- y[g == 0 & t == 1 & d_tc == lv]
      }

      mapped <- .counterfactual_quantile_map(y_target = y10, y_d00 = y00, y_d01 = y01)
      q_vals <- c(q_vals, mapped)
    }

    out$cic_num <- m_y11 - .safe_mean(q_vals)
    out$cic_den <- m_d11 - m_d10
  }

  out
}

.estimate_pair_with_cov <- function(sub_df, y_name, d_true_name, d_tc_name, cov_types, opts) {
  y <- as.numeric(sub_df[[y_name]])
  d_true <- as.numeric(sub_df[[d_true_name]])
  d_tc <- as.numeric(sub_df[[d_tc_name]])
  g <- as.integer(sub_df$.g_binary)
  t <- as.integer(sub_df$.t_binary)

  methods <- .resolve_model_methods(opts$modelx, y_vals = y, d_vals = d_true)

  resolved_sieve <- .resolve_sieve_orders(
    df = sub_df,
    y = y,
    d = d_true,
    cov_types = cov_types,
    opts = opts
  )

  x_y <- .build_feature_frame(
    sub_df,
    cov_types$continuous,
    cov_types$qualitative,
    order = resolved_sieve$y_order
  )
  x_d <- .build_feature_frame(
    sub_df,
    cov_types$continuous,
    cov_types$qualitative,
    order = resolved_sieve$d_order
  )

  idx11 <- g == 1 & t == 1
  idx10 <- g == 1 & t == 0
  idx01 <- g == 0 & t == 1
  idx00 <- g == 0 & t == 0

  m_y11 <- .safe_mean(y[idx11])
  m_d11 <- .safe_mean(d_true[idx11])

  m_y10 <- .fit_predict_mean(y[idx10], x_y[idx10, , drop = FALSE], x_y[idx11, , drop = FALSE], methods$y)
  m_d10 <- .fit_predict_mean(d_true[idx10], x_d[idx10, , drop = FALSE], x_d[idx11, , drop = FALSE], methods$d1)

  out <- list(
    did_num = NA_real_,
    did_den = NA_real_,
    tc_num = NA_real_,
    tc_den = NA_real_,
    tc_inf = NA_real_,
    tc_sup = NA_real_,
    cic_num = NA_real_,
    cic_den = NA_real_,
    special_case = FALSE,
    sieveorder_selected = c(
      outcome = resolved_sieve$y_order,
      treatment = resolved_sieve$d_order
    ),
    sieveorder_selected_by_cv = isTRUE(resolved_sieve$selected_by_cv)
  )

  if (opts$did) {
    m_y01 <- .fit_predict_mean(y[idx01], x_y[idx01, , drop = FALSE], x_y[idx11, , drop = FALSE], methods$y)
    m_y00 <- .fit_predict_mean(y[idx00], x_y[idx00, , drop = FALSE], x_y[idx11, , drop = FALSE], methods$y)

    m_d01 <- .fit_predict_mean(d_true[idx01], x_d[idx01, , drop = FALSE], x_d[idx11, , drop = FALSE], methods$d1)
    m_d00 <- .fit_predict_mean(d_true[idx00], x_d[idx00, , drop = FALSE], x_d[idx11, , drop = FALSE], methods$d1)

    out$did_num <- m_y11 - m_y10 - (m_y01 - m_y00)
    out$did_den <- m_d11 - m_d10 - (m_d01 - m_d00)
  }

  if (opts$tc) {
    special_case <- length(unique(d_tc[g == 0])) == 1L
    out$special_case <- special_case

    if (special_case) {
      if (is.finite(out$did_num) && is.finite(out$did_den)) {
        out$tc_num <- out$did_num
        out$tc_den <- out$did_den
      } else {
        out$tc_num <- m_y11 - m_y10
        out$tc_den <- m_d11 - m_d10
      }
    } else {
      levels_10 <- sort(unique(d_tc[idx10]))
      build_vec <- rep(0, sum(idx11))
      valid <- TRUE

      for (lv in levels_10) {
        p_pred <- .fit_predict_vector(
          train_y = as.numeric(d_tc[idx10] == lv),
          train_x = x_d[idx10, , drop = FALSE],
          pred_x = x_d[idx11, , drop = FALSE],
          method = methods$d2
        )
        p_pred <- pmin(pmax(p_pred, 0), 1)

        y01_mask <- idx01 & d_tc == lv
        y00_mask <- idx00 & d_tc == lv

        y01_pred <- .fit_predict_vector(
          train_y = y[y01_mask],
          train_x = x_y[y01_mask, , drop = FALSE],
          pred_x = x_y[idx11, , drop = FALSE],
          method = methods$y
        )
        y00_pred <- .fit_predict_vector(
          train_y = y[y00_mask],
          train_x = x_y[y00_mask, , drop = FALSE],
          pred_x = x_y[idx11, , drop = FALSE],
          method = methods$y
        )

        if (
          !all(is.finite(p_pred)) ||
          !all(is.finite(y01_pred)) ||
          !all(is.finite(y00_pred))
        ) {
          valid <- FALSE
          break
        }

        build_vec <- build_vec + p_pred * (y01_pred - y00_pred)
      }

      out$tc_num <- if (valid) m_y11 - m_y10 - .safe_mean(build_vec) else NA_real_
      out$tc_den <- m_d11 - m_d10
    }
  }

  out
}

.estimate_lqte_no_cov <- function(sub_df, y_name, d_true_name) {
  y <- as.numeric(sub_df[[y_name]])
  d_true <- as.numeric(sub_df[[d_true_name]])
  d <- as.integer(factor(d_true, levels = sort(unique(d_true))))
  g <- as.integer(sub_df$.g_binary)
  t <- as.integer(sub_df$.t_binary)

  idx11 <- g == 1 & t == 1
  idx10 <- g == 1 & t == 0
  idx01 <- g == 0 & t == 1
  idx00 <- g == 0 & t == 0

  q_grid <- seq(0.05, 0.95, by = 0.05)

  mean_d11 <- .safe_mean(d_true[idx11])
  mean_d10 <- .safe_mean(d_true[idx10])
  den <- mean_d11 - mean_d10

  if (!is.finite(den) || den == 0 || length(unique(d)) != 2L) {
    out <- rep(NA_real_, length(q_grid))
    names(out) <- sprintf("%.2f", q_grid)
    return(out)
  }

  y_range <- range(y, na.rm = TRUE)
  grid <- sort(unique(c(y, y_range[[1L]] - max(1, diff(y_range)))))
  max_y11 <- max(y[idx11], na.rm = TRUE) - 0.001

  q_mat <- matrix(NA_real_, nrow = length(q_grid), ncol = 2L)

  for (status in 1:2) {
    f_status <- rep(0, length(grid))

    y11 <- y[d == status & idx11]
    if (length(y11) > 0L) {
      f_11 <- .cdf_at(y11, grid)
      if (status == 1L) {
        f_status <- (1 - mean_d11) * f_11 / (1 - mean_d11 - (1 - mean_d10))
      } else {
        f_status <- mean_d11 * f_11 / den
      }
    }

    y10 <- y[d == status & idx10]
    y00 <- y[d == status & idx00]
    y01 <- y[d == status & idx01]
    if (length(y10) == 0L || length(y00) == 0L || length(y01) == 0L) {
      out <- rep(NA_real_, length(q_grid))
      names(out) <- sprintf("%.2f", q_grid)
      return(out)
    }

    inv_q <- .counterfactual_inverse_map(
      y_target = grid,
      y_from = y01,
      y_to = y00
    )
    f_10 <- .cdf_at(y10, inv_q)
    if (status == 1L) {
      f_status <- f_status - (1 - mean_d10) * f_10 / (1 - mean_d11 - (1 - mean_d10))
    } else {
      f_status <- f_status - mean_d10 * f_10 / den
    }

    f_status <- .rearrange_cdf(f_status)
    f_status[grid >= max_y11] <- 1

    for (i in seq_along(q_grid)) {
      idx <- which(f_status >= q_grid[[i]])
      q_mat[[i, status]] <- if (length(idx) == 0L) NA_real_ else min(grid[idx])
    }
  }

  out <- q_mat[, 2L] - q_mat[, 1L]

  names(out) <- sprintf("%.2f", q_grid)
  out
}

.extract_scalar <- function(x, name) {
  val <- x[[name]]
  if (is.null(val) || length(val) == 0L) return(NA_real_)
  as.numeric(val[[1L]])
}

.aggregate_late <- function(pair_results, opts) {
  late <- numeric(0)

  if (opts$did) {
    did_num <- vapply(pair_results, .extract_scalar, numeric(1), name = "did_num")
    did_den <- vapply(pair_results, .extract_scalar, numeric(1), name = "did_den")
    keep <- is.finite(did_num) & is.finite(did_den) & did_den != 0
    if (!any(keep)) {
      stop("Given the data structure, impossible to estimate DID.", call. = FALSE)
    }

    if (opts$numerator) {
      late[["DID_num"]] <- sum(did_num[keep])
    } else {
      late[["W_DID"]] <- sum(did_num[keep]) / sum(did_den[keep])
    }
  }

  if (opts$tc || opts$partial) {
    if (opts$partial) {
      tc_inf <- vapply(pair_results, .extract_scalar, numeric(1), name = "tc_inf")
      tc_sup <- vapply(pair_results, .extract_scalar, numeric(1), name = "tc_sup")
      tc_den <- vapply(pair_results, .extract_scalar, numeric(1), name = "tc_den")
      keep <- is.finite(tc_inf) & is.finite(tc_sup) & is.finite(tc_den) & tc_den != 0
      if (!any(keep)) {
        stop("Given the data structure, partial cannot be estimated.", call. = FALSE)
      }
      inf_est <- sum(tc_inf[keep] * tc_den[keep]) / sum(tc_den[keep])
      sup_est <- sum(tc_sup[keep] * tc_den[keep]) / sum(tc_den[keep])
      bounds <- sort(c(inf_est, sup_est))
      late[["TC_inf"]] <- bounds[[1L]]
      late[["TC_sup"]] <- bounds[[2L]]
    } else {
      tc_num <- vapply(pair_results, .extract_scalar, numeric(1), name = "tc_num")
      tc_den <- vapply(pair_results, .extract_scalar, numeric(1), name = "tc_den")
      keep <- is.finite(tc_num) & is.finite(tc_den) & tc_den != 0
      if (!any(keep)) {
        stop("Given the data structure, impossible to estimate TC.", call. = FALSE)
      }

      if (opts$numerator) {
        late[["TC_num"]] <- sum(tc_num[keep])
      } else {
        late[["W_TC"]] <- sum(tc_num[keep]) / sum(tc_den[keep])
      }
    }
  }

  if (opts$cic) {
    cic_num <- vapply(pair_results, .extract_scalar, numeric(1), name = "cic_num")
    cic_den <- vapply(pair_results, .extract_scalar, numeric(1), name = "cic_den")
    keep <- is.finite(cic_num) & is.finite(cic_den) & cic_den != 0
    if (!any(keep)) {
      stop("Given the data structure, impossible to estimate CIC.", call. = FALSE)
    }

    if (opts$numerator) {
      late[["CIC_num"]] <- sum(cic_num[keep])
    } else {
      late[["W_CIC"]] <- sum(cic_num[keep]) / sum(cic_den[keep])
    }
  }

  late
}

.estimate_all_native <- function(df, prepared, opts) {
  .validate_native_request(df, prepared, opts)

  df$.d_true_native <- as.numeric(df[[prepared$d_name]])
  df$.d_tc_native <- .recategorize_treatment(df$.d_true_native, opts$newcateg)

  pair_frames <- .build_pair_frames(
    df = df,
    group = prepared$group,
    time = prepared$time,
    group_forward = prepared$group_forward
  )

  if (opts$lqte) {
    for (pair_df in pair_frames) {
      signs <- sort(unique(pair_df$.g_star[pair_df$.g_star != 0 & pair_df$.g_star != -3]))
      if (length(signs) != 1L) {
        stop("lqte cannot be used when there are more than two groups.", call. = FALSE)
      }
    }
  }

  if (opts$partial) {
    for (pair_df in pair_frames) {
      signs <- sort(unique(pair_df$.g_star[pair_df$.g_star != 0 & pair_df$.g_star != -3]))
      if (length(signs) != 1L) {
        stop("partial can only be used with two treatment groups.", call. = FALSE)
      }
    }
  }

  pair_results <- list()
  lqte <- numeric(0)
  sieve_selection <- NULL

  for (pair_df in pair_frames) {
    subdesigns <- .extract_pair_subdesigns(pair_df, prepared$time)
    if (length(subdesigns) == 0L) {
      next
    }

    for (sub in subdesigns) {
      if (length(prepared$covariates) == 0L) {
        est <- .estimate_pair_no_cov(
          sub_df = sub$data,
          y_name = prepared$y_name,
          d_true_name = ".d_true_native",
          d_tc_name = ".d_tc_native",
          opts = opts
        )
      } else {
        est <- .estimate_pair_with_cov(
          sub_df = sub$data,
          y_name = prepared$y_name,
          d_true_name = ".d_true_native",
          d_tc_name = ".d_tc_native",
          cov_types = prepared$covariate_types,
          opts = opts
        )
      }

      if (opts$partial && isTRUE(est$special_case)) {
        stop(
          "Given the data structure, partial cannot be used in partially-sharp treatment designs.",
          call. = FALSE
        )
      }

      est$pair_id <- sub$pair_id
      est$sign <- sub$sign
      est$counts <- sub$counts
      pair_results[[length(pair_results) + 1L]] <- est
      if (!is.null(est$sieveorder_selected)) {
        sieve_selection <- est$sieveorder_selected
      }

      if (opts$lqte && length(lqte) == 0L) {
        lqte <- .estimate_lqte_no_cov(
          sub_df = sub$data,
          y_name = prepared$y_name,
          d_true_name = ".d_true_native"
        )
      }
    }
  }

  if (length(pair_results) == 0L) {
    stop(
      "Given the data structure, impossible to estimate tc, cic or lqte. Try recategorizing treatment with `newcateg`.",
      call. = FALSE
    )
  }

  late <- .aggregate_late(pair_results, opts)

  if (opts$lqte && length(lqte) == 0L) {
    stop("lqte could not be estimated from the supplied design.", call. = FALSE)
  }

  list(
    late = late,
    lqte = lqte,
    pair_results = pair_results,
    sieveorder_selected = sieve_selection
  )
}

.bootstrap_native <- function(df, prepared, opts, point) {
  boot_sentinel <- 1000000000000000
  if (opts$nose) return(NULL)
  # Lock bootstrap RNG semantics to avoid cross-version drift.
  old_rng <- RNGkind()
  restore_seed <- !is.null(opts$seed)
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (restore_seed && old_seed_exists) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    do.call(RNGkind, as.list(old_rng))
    if (restore_seed) {
      if (old_seed_exists) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }
  }, add = TRUE)
  RNGversion("4.3.0")
  if (restore_seed) {
    set.seed(as.integer(opts$seed))
  }

  late_names <- names(point$late)
  lqte_names <- names(point$lqte)
  n_reps <- as.integer(opts$breps)
  n_misreps <- 0L

  reps_late <- if (length(late_names) > 0L) {
    matrix(NA_real_, nrow = n_reps, ncol = length(late_names), dimnames = list(NULL, late_names))
  } else {
    NULL
  }

  reps_lqte <- if (length(lqte_names) > 0L) {
    matrix(NA_real_, nrow = n_reps, ncol = length(lqte_names), dimnames = list(NULL, lqte_names))
  } else {
    NULL
  }

  n <- nrow(df)
  for (b in seq_len(n_reps)) {
    degenerate_correction <- if (stats::runif(1) < 0.5) -boot_sentinel else boot_sentinel

    if (is.null(opts$cluster)) {
      idx <- sample.int(n, size = n, replace = TRUE)
    } else {
      cl <- unique(df[[opts$cluster]])
      sampled_cl <- sample(cl, size = length(cl), replace = TRUE)
      idx <- unlist(lapply(sampled_cl, function(cc) which(df[[opts$cluster]] == cc)), use.names = FALSE)
    }

    est_b <- tryCatch(
      .estimate_all_native(df[idx, , drop = FALSE], prepared = prepared, opts = opts),
      error = function(...) NULL
    )

    late_vals <- NULL
    lqte_vals <- NULL
    invalid_rep <- is.null(est_b)

    if (!invalid_rep && !is.null(reps_late)) {
      late_vals <- as.numeric(est_b$late[late_names])
      invalid_rep <- length(late_vals) != length(late_names) || !all(is.finite(late_vals))
    }

    if (!invalid_rep && !is.null(reps_lqte)) {
      lqte_vals <- as.numeric(est_b$lqte[lqte_names])
      invalid_rep <- length(lqte_vals) != length(lqte_names) || !all(is.finite(lqte_vals))
    }

    if (invalid_rep) {
      n_misreps <- n_misreps + 1L
      if (!is.null(reps_late)) reps_late[b, ] <- degenerate_correction
      if (!is.null(reps_lqte)) reps_lqte[b, ] <- degenerate_correction
      next
    }

    if (!is.null(reps_late)) {
      reps_late[b, ] <- late_vals
    }

    if (!is.null(reps_lqte)) {
      reps_lqte[b, ] <- lqte_vals
    }
  }

  list(
    reps_late = reps_late,
    reps_lqte = reps_lqte,
    n_reps = n_reps,
    n_misreps = n_misreps,
    share_failures = if (n_reps > 0L) n_misreps / n_reps else NA_real_
  )
}

.pairwise_eqtest <- function(estimates_named) {
  nm <- names(estimates_named)
  if (length(nm) < 2L) return(NULL)
  out <- list()
  for (i in seq_len(length(nm) - 1L)) {
    for (j in seq.int(i + 1L, length(nm))) {
      key <- paste0(nm[[i]], "_", nm[[j]])
      out[[key]] <- estimates_named[[i]] - estimates_named[[j]]
    }
  }
  unlist(out)
}

.calc_boot_summary <- function(reps) {
  if (is.null(reps) || ncol(reps) == 0L) {
    return(list(se = numeric(0), ci = matrix(numeric(0), ncol = 2)))
  }

  reps <- .recode_bootstrap_sentinel(reps)

  se <- apply(reps, 2, function(x) {
    if (all(is.na(x))) return(NA_real_)
    stats::sd(x, na.rm = TRUE)
  })

  ci <- t(apply(reps, 2, function(x) {
    if (all(is.na(x))) return(c(NA_real_, NA_real_))
    stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  }))

  list(se = as.numeric(se), ci = ci)
}

.calc_boot_vcov <- function(reps_late, reps_lqte) {
  reps <- NULL
  if (!is.null(reps_late) && ncol(reps_late) > 0L) {
    reps <- reps_late
  }
  if (!is.null(reps_lqte) && ncol(reps_lqte) > 0L) {
    lqte_reps <- reps_lqte
    colnames(lqte_reps) <- paste0("LQTE_", colnames(lqte_reps))
    reps <- if (is.null(reps)) lqte_reps else cbind(reps, lqte_reps)
  }
  if (is.null(reps) || ncol(reps) == 0L) {
    return(NULL)
  }

  reps <- .recode_bootstrap_sentinel(reps)
  stats::cov(reps, use = "pairwise.complete.obs")
}

.recode_bootstrap_sentinel <- function(x) {
  if (is.null(x)) return(x)
  boot_sentinel <- 1000000000000000
  x[x == boot_sentinel | x == -boot_sentinel] <- NA_real_
  x
}

.compute_counts <- function(df, prepared) {
  n11 <- NA_integer_
  n10 <- NA_integer_
  n01 <- NA_integer_
  n00 <- NA_integer_

  if (is.null(prepared$group_forward)) {
    g_vals <- sort(unique(df[[prepared$group]]))
    t_vals <- sort(unique(df[[prepared$time]]))

    if (length(g_vals) == 2L && length(t_vals) == 2L && all(g_vals %in% c(0, 1))) {
      g0 <- g_vals[[1L]]
      g1 <- g_vals[[2L]]
      t0 <- t_vals[[1L]]
      t1 <- t_vals[[2L]]
      n11 <- sum(df[[prepared$group]] == g1 & df[[prepared$time]] == t1)
      n10 <- sum(df[[prepared$group]] == g1 & df[[prepared$time]] == t0)
      n01 <- sum(df[[prepared$group]] == g0 & df[[prepared$time]] == t1)
      n00 <- sum(df[[prepared$group]] == g0 & df[[prepared$time]] == t0)
    }
  }

  list(n11 = n11, n10 = n10, n01 = n01, n00 = n00)
}

.run_native_backend <- function(prepared, opts) {
  df <- prepared$used

  point <- .estimate_all_native(df = df, prepared = prepared, opts = opts)
  bt <- .bootstrap_native(df = df, prepared = prepared, opts = opts, point = point)

  late <- data.frame(
    estimator = names(point$late),
    estimate = as.numeric(point$late),
    std.error = rep(NA_real_, length(point$late)),
    conf.low = rep(NA_real_, length(point$late)),
    conf.high = rep(NA_real_, length(point$late)),
    stringsAsFactors = FALSE
  )

  lqte <- NULL
  if (length(point$lqte) > 0L) {
    lqte <- data.frame(
      quantile = as.numeric(names(point$lqte)),
      estimate = as.numeric(point$lqte),
      std.error = rep(NA_real_, length(point$lqte)),
      conf.low = rep(NA_real_, length(point$lqte)),
      conf.high = rep(NA_real_, length(point$lqte)),
      stringsAsFactors = FALSE
    )
  }

  matrices <- list(
    b_LATE = .to_column_matrix(late$estimate, late$estimator)
  )
  n_reps <- NA_integer_
  n_misreps <- NA_integer_
  share_failures <- NA_real_
  vcov_mat <- NULL

  if (!opts$nose) {
    n_reps <- bt$n_reps
    n_misreps <- bt$n_misreps
    share_failures <- bt$share_failures
    vcov_mat <- .calc_boot_vcov(bt$reps_late, bt$reps_lqte)
    late_boot <- .calc_boot_summary(bt$reps_late)
    late$std.error <- late_boot$se
    if (length(late_boot$se) > 0L) {
      late$conf.low <- late_boot$ci[, 1]
      late$conf.high <- late_boot$ci[, 2]
    }

    matrices$se_LATE <- .to_column_matrix(late$std.error, late$estimator)
    matrices$ci_LATE <- .to_ci_matrix(late$conf.low, late$conf.high, late$estimator)

    if (!is.null(lqte)) {
      lqte_boot <- .calc_boot_summary(bt$reps_lqte)
      lqte$std.error <- lqte_boot$se
      if (length(lqte_boot$se) > 0L) {
        lqte$conf.low <- lqte_boot$ci[, 1]
        lqte$conf.high <- lqte_boot$ci[, 2]
      }

      rn <- sprintf("%.2f", lqte$quantile)
      matrices$b_LQTE <- .to_column_matrix(lqte$estimate, rn)
      matrices$se_LQTE <- .to_column_matrix(lqte$std.error, rn)
      matrices$ci_LQTE <- .to_ci_matrix(lqte$conf.low, lqte$conf.high, rn)
    }
  } else if (!is.null(lqte)) {
    rn <- sprintf("%.2f", lqte$quantile)
    matrices$b_LQTE <- .to_column_matrix(lqte$estimate, rn)
  }

  eqtest <- NULL
  if (opts$eqtest && !opts$numerator && nrow(late) > 1L) {
    eq_point <- .pairwise_eqtest(stats::setNames(late$estimate, late$estimator))
    if (!is.null(eq_point)) {
      eqtest <- data.frame(
        contrast = names(eq_point),
        estimate = as.numeric(eq_point),
        std.error = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        stringsAsFactors = FALSE
      )

      if (!opts$nose && !is.null(bt$reps_late)) {
        reps_late_clean <- .recode_bootstrap_sentinel(bt$reps_late)
        eq_reps <- t(apply(reps_late_clean, 1, function(row) {
          .pairwise_eqtest(stats::setNames(as.numeric(row), colnames(reps_late_clean)))
        }))

        if (!is.null(eq_reps) && ncol(eq_reps) > 0L) {
          eq_boot <- .calc_boot_summary(eq_reps)
          eqtest$std.error <- eq_boot$se
          eqtest$conf.low <- eq_boot$ci[, 1]
          eqtest$conf.high <- eq_boot$ci[, 2]
        }
      }

      matrices$b_LATE_eqtest <- .to_column_matrix(eqtest$estimate, eqtest$contrast)
      if (!opts$nose) {
        matrices$se_LATE_eqtest <- .to_column_matrix(eqtest$std.error, eqtest$contrast)
        matrices$ci_LATE_eqtest <- .to_ci_matrix(eqtest$conf.low, eqtest$conf.high, eqtest$contrast)
      }
    }
  }

  counts <- .compute_counts(df, prepared)
  tag_mask <- if (opts$tagobs) prepared$mask else NULL

  list(
    late = late,
    eqtest = eqtest,
    lqte = lqte,
    matrices = matrices,
    vcov = vcov_mat,
    tagobs = tag_mask,
    n = nrow(df),
    n11 = counts$n11,
    n10 = counts$n10,
    n01 = counts$n01,
    n00 = counts$n00,
    sieveorder_selected = point$sieveorder_selected,
    n_reps = n_reps,
    n_misreps = n_misreps,
    share_failures = share_failures
  )
}

#' @title summary.fuzzydid
#' @description Print a compact summary table for fuzzydid results.
#' @param object A fuzzydid object.
#' @param ... Unused.
#' @return The input \code{object}, returned invisibly with class
#'   \code{"fuzzydid"}, after printing the available estimator tables. This
#'   method is called for its side effect of displaying the \code{late},
#'   \code{eqtest}, and \code{lqte} components in a compact tabular form.
#' @examples
#' make_example_cell <- function(g, t, ones, n_cell = 20L) {
#'   data.frame(
#'     g = rep.int(g, n_cell),
#'     t = rep.int(t, n_cell),
#'     d = c(rep.int(1L, ones), rep.int(0L, n_cell - ones))
#'   )
#' }
#'
#' df <- rbind(
#'   make_example_cell(0L, 0L, 4L),
#'   make_example_cell(0L, 1L, 8L),
#'   make_example_cell(1L, 0L, 6L),
#'   make_example_cell(1L, 1L, 16L)
#' )
#' df$id <- seq_len(nrow(df))
#' df$y <- 1 + 0.5 * df$g + 0.4 * df$t + 2 * df$d + sin(df$id / 7)
#'
#' fit <- fuzzydid(
#'   data = df[, c("y", "g", "t", "d")],
#'   formula = y ~ d,
#'   group = "g",
#'   time = "t",
#'   did = TRUE,
#'   nose = TRUE
#' )
#'
#' summary(fit)
#' @export
summary.fuzzydid <- function(object, ...) {
  if (!inherits(object, "fuzzydid")) {
    stop("`object` must be a fuzzydid object.", call. = FALSE)
  }

  if (!is.null(object$late) && nrow(object$late) > 0L) {
    cat("LATE estimators\n")
    print(knitr::kable(object$late, row.names = FALSE))
  } else {
    cat("No LATE estimators returned.\n")
  }

  if (!is.null(object$eqtest) && nrow(object$eqtest) > 0L) {
    cat("\nEquality tests\n")
    print(knitr::kable(object$eqtest, row.names = FALSE))
  }

  if (!is.null(object$lqte) && nrow(object$lqte) > 0L) {
    cat("\nLQTE estimators\n")
    print(knitr::kable(object$lqte, row.names = FALSE))
  }

  invisible(object)
}

.coef_fuzzydid <- function(x) {
  vals <- numeric(0)
  if (!is.null(x$late) && nrow(x$late) > 0L) {
    vals <- stats::setNames(x$late$estimate, x$late$estimator)
  }
  if (!is.null(x$lqte) && nrow(x$lqte) > 0L) {
    lqte_vals <- stats::setNames(
      x$lqte$estimate,
      paste0("LQTE_", sprintf("%.2f", x$lqte$quantile))
    )
    vals <- c(vals, lqte_vals)
  }
  vals
}

.confint_fuzzydid <- function(x) {
  out <- matrix(
    numeric(0),
    nrow = 0L,
    ncol = 2L,
    dimnames = list(character(0), c("2.5 %", "97.5 %"))
  )
  if (!is.null(x$late) && nrow(x$late) > 0L) {
    out <- rbind(
      out,
      matrix(
        c(x$late$conf.low, x$late$conf.high),
        ncol = 2L,
        dimnames = list(x$late$estimator, c("2.5 %", "97.5 %"))
      )
    )
  }
  if (!is.null(x$lqte) && nrow(x$lqte) > 0L) {
    rn <- paste0("LQTE_", sprintf("%.2f", x$lqte$quantile))
    out <- rbind(
      out,
      matrix(
        c(x$lqte$conf.low, x$lqte$conf.high),
        ncol = 2L,
        dimnames = list(rn, c("2.5 %", "97.5 %"))
      )
    )
  }
  out
}

#' @title print.fuzzydid
#' @description Print a compact fuzzydid object header and estimator table.
#' @param x A fuzzydid object.
#' @param ... Unused.
#' @return The input \code{x}, returned invisibly.
#' @examples
#' df <- expand.grid(i = seq_len(20), g = 0:1, t = 0:1)
#' df$d <- as.integer(df$i <= c(4, 8, 6, 16)[1 + df$t + 2 * df$g])
#' df$y <- 1 + 0.5 * df$g + 0.4 * df$t + 2 * df$d + sin(df$i / 7)
#'
#' fit <- fuzzydid(df, y ~ d, group = "g", time = "t", did = TRUE, nose = TRUE)
#' print(fit)
#' @export
print.fuzzydid <- function(x, ...) {
  if (!inherits(x, "fuzzydid")) {
    stop("`x` must be a fuzzydid object.", call. = FALSE)
  }
  cat(sprintf("fuzzydid fit: %s observations", x$n), "\n")
  if (!is.null(x$late) && nrow(x$late) > 0L) {
    print(x$late, row.names = FALSE)
  }
  if (!is.null(x$lqte) && nrow(x$lqte) > 0L) {
    cat("\nLQTE estimators: ", nrow(x$lqte), " quantiles\n", sep = "")
  }
  invisible(x)
}

#' @title coef.fuzzydid
#' @description Extract fuzzydid point estimates.
#' @param object A fuzzydid object.
#' @param ... Unused.
#' @return A named numeric vector of LATE and LQTE point estimates.
#' @examples
#' df <- expand.grid(i = seq_len(20), g = 0:1, t = 0:1)
#' df$d <- as.integer(df$i <= c(4, 8, 6, 16)[1 + df$t + 2 * df$g])
#' df$y <- 1 + 0.5 * df$g + 0.4 * df$t + 2 * df$d + sin(df$i / 7)
#'
#' fit <- fuzzydid(df, y ~ d, group = "g", time = "t", did = TRUE, nose = TRUE)
#' coef(fit)
#' @export
coef.fuzzydid <- function(object, ...) {
  if (!inherits(object, "fuzzydid")) {
    stop("`object` must be a fuzzydid object.", call. = FALSE)
  }
  .coef_fuzzydid(object)
}

#' @title confint.fuzzydid
#' @description Extract stored percentile bootstrap confidence intervals.
#' @param object A fuzzydid object.
#' @param parm Optional parameter subset.
#' @param level Confidence level. Only \code{0.95} is currently stored.
#' @param ... Unused.
#' @return A matrix with lower and upper confidence limits.
#' @examples
#' df <- expand.grid(i = seq_len(20), g = 0:1, t = 0:1)
#' df$d <- as.integer(df$i <= c(4, 8, 6, 16)[1 + df$t + 2 * df$g])
#' df$y <- 1 + 0.5 * df$g + 0.4 * df$t + 2 * df$d + sin(df$i / 7)
#'
#' fit <- fuzzydid(df, y ~ d, group = "g", time = "t", did = TRUE, nose = TRUE)
#' confint(fit)
#' @export
confint.fuzzydid <- function(object, parm, level = 0.95, ...) {
  if (!inherits(object, "fuzzydid")) {
    stop("`object` must be a fuzzydid object.", call. = FALSE)
  }
  if (!identical(as.numeric(level), 0.95)) {
    stop("Only the stored 95% confidence intervals are available.", call. = FALSE)
  }
  out <- .confint_fuzzydid(object)
  if (!missing(parm)) {
    out <- out[parm, , drop = FALSE]
  }
  out
}

#' @title nobs.fuzzydid
#' @description Extract the estimation sample size.
#' @param object A fuzzydid object.
#' @param ... Unused.
#' @return Integer number of observations used for estimation.
#' @examples
#' df <- expand.grid(i = seq_len(20), g = 0:1, t = 0:1)
#' df$d <- as.integer(df$i <= c(4, 8, 6, 16)[1 + df$t + 2 * df$g])
#' df$y <- 1 + 0.5 * df$g + 0.4 * df$t + 2 * df$d + sin(df$i / 7)
#'
#' fit <- fuzzydid(df, y ~ d, group = "g", time = "t", did = TRUE, nose = TRUE)
#' nobs(fit)
#' @export
nobs.fuzzydid <- function(object, ...) {
  if (!inherits(object, "fuzzydid")) {
    stop("`object` must be a fuzzydid object.", call. = FALSE)
  }
  as.integer(object$n)
}

#' @title formula.fuzzydid
#' @description Extract the formula used to fit a fuzzydid object.
#' @param x A fuzzydid object.
#' @param ... Unused.
#' @return The original formula.
#' @examples
#' df <- expand.grid(i = seq_len(20), g = 0:1, t = 0:1)
#' df$d <- as.integer(df$i <= c(4, 8, 6, 16)[1 + df$t + 2 * df$g])
#' df$y <- 1 + 0.5 * df$g + 0.4 * df$t + 2 * df$d + sin(df$i / 7)
#'
#' fit <- fuzzydid(df, y ~ d, group = "g", time = "t", did = TRUE, nose = TRUE)
#' formula(fit)
#' @export
formula.fuzzydid <- function(x, ...) {
  if (!inherits(x, "fuzzydid")) {
    stop("`x` must be a fuzzydid object.", call. = FALSE)
  }
  x$formula
}

#' @title vcov.fuzzydid
#' @description Extract the bootstrap covariance matrix for stored estimates.
#' @param object A fuzzydid object.
#' @param ... Unused.
#' @return A covariance matrix for \code{coef(object)}.
#' @examples
#' df <- expand.grid(i = seq_len(20), g = 0:1, t = 0:1)
#' df$d <- as.integer(df$i <= c(4, 8, 6, 16)[1 + df$t + 2 * df$g])
#' df$y <- 1 + 0.5 * df$g + 0.4 * df$t + 2 * df$d + sin(df$i / 7)
#'
#' fit <- fuzzydid(df, y ~ d, group = "g", time = "t", did = TRUE,
#'                 breps = 5, seed = 1)
#' vcov(fit)
#' @export
vcov.fuzzydid <- function(object, ...) {
  if (!inherits(object, "fuzzydid")) {
    stop("`object` must be a fuzzydid object.", call. = FALSE)
  }
  if (is.null(object$vcov)) {
    stop("Bootstrap covariance is unavailable; refit with `nose = FALSE`.", call. = FALSE)
  }
  object$vcov
}

#' @title plot.fuzzydid
#' @description Plot fuzzydid point estimates and stored confidence intervals.
#' @param x A fuzzydid object.
#' @param ... Unused.
#' @return The input \code{x}, returned invisibly.
#' @examples
#' df <- expand.grid(i = seq_len(20), g = 0:1, t = 0:1)
#' df$d <- as.integer(df$i <= c(4, 8, 6, 16)[1 + df$t + 2 * df$g])
#' df$y <- 1 + 0.5 * df$g + 0.4 * df$t + 2 * df$d + sin(df$i / 7)
#'
#' fit <- fuzzydid(df, y ~ d, group = "g", time = "t", did = TRUE, nose = TRUE)
#' plot(fit)
#' @export
plot.fuzzydid <- function(x, ...) {
  if (!inherits(x, "fuzzydid")) {
    stop("`x` must be a fuzzydid object.", call. = FALSE)
  }

  est <- .coef_fuzzydid(x)
  ci <- .confint_fuzzydid(x)
  if (length(est) == 0L) {
    stop("No estimates are available to plot.", call. = FALSE)
  }

  y_pos <- seq_along(est)
  xlim <- range(c(est, ci), na.rm = TRUE)
  graphics::plot(
    est,
    y_pos,
    xlim = xlim,
    yaxt = "n",
    ylab = "",
    xlab = "Estimate",
    pch = 19
  )
  graphics::axis(2, at = y_pos, labels = names(est), las = 1)
  ok <- stats::complete.cases(ci)
  graphics::segments(ci[ok, 1L], y_pos[ok], ci[ok, 2L], y_pos[ok])
  graphics::abline(v = 0, lty = 3, col = "gray50")
  invisible(x)
}

#' @title tidy.fuzzydid
#' @description Tidy extractor for `fuzzydid` objects.
#' @param x A fuzzydid object.
#' @param ... Unused.
#' @return A data frame with class \code{"data.frame"} and one row per
#'   available estimate or contrast. \code{component} identifies whether the
#'   row comes from the LATE table, equality-test table, or LQTE table.
#'   \code{model} and \code{term} label the estimator or contrast.
#'   \code{estimate} is the point estimate, while \code{std.error},
#'   \code{conf.low}, and \code{conf.high} contain bootstrap uncertainty
#'   summaries when available and \code{NA} otherwise. If no estimates are
#'   available, an empty data frame with the same columns is returned.
#' @examples
#' make_example_cell <- function(g, t, ones, n_cell = 20L) {
#'   data.frame(
#'     g = rep.int(g, n_cell),
#'     t = rep.int(t, n_cell),
#'     d = c(rep.int(1L, ones), rep.int(0L, n_cell - ones))
#'   )
#' }
#'
#' df <- rbind(
#'   make_example_cell(0L, 0L, 4L),
#'   make_example_cell(0L, 1L, 8L),
#'   make_example_cell(1L, 0L, 6L),
#'   make_example_cell(1L, 1L, 16L)
#' )
#' df$id <- seq_len(nrow(df))
#' df$y <- 1 + 0.5 * df$g + 0.4 * df$t + 2 * df$d + sin(df$id / 7)
#'
#' fit <- fuzzydid(
#'   data = df[, c("y", "g", "t", "d")],
#'   formula = y ~ d,
#'   group = "g",
#'   time = "t",
#'   did = TRUE,
#'   nose = TRUE
#' )
#'
#' generics::tidy(fit)
#' @export
tidy.fuzzydid <- function(x, ...) {
  if (!inherits(x, "fuzzydid")) {
    stop("`x` must be a fuzzydid object.", call. = FALSE)
  }

  rows <- list()

  if (!is.null(x$late) && nrow(x$late) > 0L) {
    rows[[length(rows) + 1L]] <- data.frame(
      component = "late",
      model = x$late$estimator,
      term = x$late$estimator,
      estimate = x$late$estimate,
      std.error = x$late$std.error,
      conf.low = x$late$conf.low,
      conf.high = x$late$conf.high,
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(x$eqtest) && nrow(x$eqtest) > 0L) {
    rows[[length(rows) + 1L]] <- data.frame(
      component = "eqtest",
      model = "eqtest",
      term = x$eqtest$contrast,
      estimate = x$eqtest$estimate,
      std.error = x$eqtest$std.error,
      conf.low = x$eqtest$conf.low,
      conf.high = x$eqtest$conf.high,
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(x$lqte) && nrow(x$lqte) > 0L) {
    rows[[length(rows) + 1L]] <- data.frame(
      component = "lqte",
      model = "lqte",
      term = x$lqte$quantile,
      estimate = x$lqte$estimate,
      std.error = x$lqte$std.error,
      conf.low = x$lqte$conf.low,
      conf.high = x$lqte$conf.high,
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0L) {
    return(data.frame(
      component = character(0),
      model = character(0),
      term = character(0),
      estimate = numeric(0),
      std.error = numeric(0),
      conf.low = numeric(0),
      conf.high = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, rows)
}

#' @title glance.fuzzydid
#' @description One-row summary for `fuzzydid` objects.
#' @param x A fuzzydid object.
#' @param ... Unused.
#' @return A one-row data frame with class \code{"data.frame"} summarizing the
#'   fitted object. \code{backend} reports the computation path and
#'   \code{Num.Obs.} reports the estimation sample size. \code{N.11},
#'   \code{N.10}, \code{N.01}, and \code{N.00} give the four design cell
#'   counts. \code{N.reps}, \code{N.misreps}, and \code{Share.failures}
#'   describe bootstrap replication totals and failure rates, or are
#'   \code{NA} when \code{nose = TRUE}.
#' @examples
#' make_example_cell <- function(g, t, ones, n_cell = 20L) {
#'   data.frame(
#'     g = rep.int(g, n_cell),
#'     t = rep.int(t, n_cell),
#'     d = c(rep.int(1L, ones), rep.int(0L, n_cell - ones))
#'   )
#' }
#'
#' df <- rbind(
#'   make_example_cell(0L, 0L, 4L),
#'   make_example_cell(0L, 1L, 8L),
#'   make_example_cell(1L, 0L, 6L),
#'   make_example_cell(1L, 1L, 16L)
#' )
#' df$id <- seq_len(nrow(df))
#' df$y <- 1 + 0.5 * df$g + 0.4 * df$t + 2 * df$d + sin(df$id / 7)
#'
#' fit <- fuzzydid(
#'   data = df[, c("y", "g", "t", "d")],
#'   formula = y ~ d,
#'   group = "g",
#'   time = "t",
#'   did = TRUE,
#'   nose = TRUE
#' )
#'
#' generics::glance(fit)
#' @export
glance.fuzzydid <- function(x, ...) {
  if (!inherits(x, "fuzzydid")) {
    stop("`x` must be a fuzzydid object.", call. = FALSE)
  }

  data.frame(
    backend = x$backend,
    Num.Obs. = x$n,
    N.11 = x$n11,
    N.10 = x$n10,
    N.01 = x$n01,
    N.00 = x$n00,
    N.reps = x$n_reps,
    N.misreps = x$n_misreps,
    Share.failures = x$share_failures,
    stringsAsFactors = FALSE
  )
}
