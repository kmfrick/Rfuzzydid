#' @title fuzzydid
#' @description
#' Formula-first interface for fuzzy difference-in-differences estimators.
#' Estimation is fully native in R.
#' @param data Data frame.
#' @param formula Formula of the form `y ~ d + covariates`.
#' @param group Name of the group variable (backward group for multi-period).
#' @param time Name of the time variable.
#' @param group_forward Optional name of the forward group variable for
#'   multi-period designs.
#' @param did,tc,cic,lqte Logical flags indicating estimators to compute.
#' @param newcateg Optional numeric vector of upper bounds used to recategorize
#'   treatment values for TC/CIC.
#' @param numerator Logical; return estimator numerators for DID/TC/CIC.
#' @param partial Logical; request TC partial-identification bounds.
#' @param nose Logical; skip bootstrap standard errors and confidence intervals.
#' @param cluster Optional name of cluster variable for clustered bootstrap.
#' @param breps Number of bootstrap replications.
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
#' @param backend One of `"auto"` or `"native"`. `"stata"` is no longer
#'   supported and errors.
#' @param seed Preserved for API compatibility. Bootstrap inference uses a
#'   fixed Stata-parity seed (`1`) when `nose = FALSE`.
#' @return An object of class `fuzzydid`.
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
  backend = c("auto", "native", "stata"),
  seed = 1
) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula like y ~ d + x1 + x2.", call. = FALSE)
  }

  backend <- match.arg(backend)
  breps_missing <- missing(breps)

  opts <- list(
    did = isTRUE(did),
    tc = isTRUE(tc),
    cic = isTRUE(cic),
    lqte = isTRUE(lqte),
    newcateg = newcateg,
    numerator = isTRUE(numerator),
    partial = isTRUE(partial),
    nose = isTRUE(nose),
    cluster = cluster,
    breps = breps,
    eqtest = isTRUE(eqtest),
    modelx = modelx,
    sieves = isTRUE(sieves),
    sieveorder = sieveorder,
    tagobs = isTRUE(tagobs),
    seed = seed
  )

  .validate_common_options(opts, breps_missing = breps_missing)
  parsed <- .parse_fuzzy_formula(formula, data)
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

  selected_backend <- .choose_backend(backend = backend)
  out <- .run_native_backend(
    prepared = prepared,
    opts = opts
  )

  out$backend <- selected_backend
  out$formula <- formula
  out$group <- group
  out$time <- time
  out$group_forward <- group_forward
  out$call <- match.call()
  out$options <- opts
  class(out) <- "fuzzydid"
  out
}

.validate_common_options <- function(opts, breps_missing) {
  if (!opts$did && !opts$tc && !opts$cic && !opts$lqte) {
    stop(
      "At least one estimator must be requested: did, tc, cic, or lqte.",
      call. = FALSE
    )
  }

  if (opts$nose && !breps_missing) {
    stop("`nose = TRUE` cannot be used with explicit `breps`.", call. = FALSE)
  }

  if (opts$nose && !is.null(opts$cluster)) {
    stop("`nose = TRUE` cannot be used with `cluster`.", call. = FALSE)
  }

  if (opts$eqtest && (opts$numerator || opts$nose)) {
    stop("`eqtest` cannot be used with `numerator` or `nose`.", call. = FALSE)
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
    if (any(!is.finite(opts$sieveorder))) {
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

  if (!opts$nose) {
    if (!is.numeric(opts$breps) || length(opts$breps) != 1L || is.na(opts$breps) || opts$breps < 2) {
      stop("`breps` must be a numeric scalar >= 2.", call. = FALSE)
    }
  }
}

.parse_fuzzy_formula <- function(formula, data) {
  tr <- stats::terms(formula, data = data)
  rhs_terms <- attr(tr, "term.labels")
  if (length(rhs_terms) < 1L) {
    stop("Formula RHS must include treatment variable `d`.", call. = FALSE)
  }

  lhs_vars <- all.vars(formula[[2L]])
  if (length(lhs_vars) != 1L) {
    stop("Formula LHS must contain exactly one outcome variable.", call. = FALSE)
  }

  y_name <- lhs_vars[[1L]]
  d_name <- rhs_terms[[1L]]
  covariates <- rhs_terms[-1L]

  needed <- c(y_name, d_name, covariates)
  missing_vars <- setdiff(needed, names(data))
  if (length(missing_vars) > 0L) {
    stop(
      sprintf("Variables not found in `data`: %s", paste(missing_vars, collapse = ", ")),
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

.prepare_input_data <- function(data, y_name, d_name, covariates, group, time, group_forward, cluster) {
  role_vars <- unique(c(y_name, d_name, group, time, group_forward, covariates, cluster))
  role_vars <- role_vars[!is.na(role_vars) & nzchar(role_vars)]

  missing_vars <- setdiff(role_vars, names(data))
  if (length(missing_vars) > 0L) {
    stop(
      sprintf("Required variable(s) not found in `data`: %s", paste(missing_vars, collapse = ", ")),
      call. = FALSE
    )
  }

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

.choose_backend <- function(backend) {
  if (backend == "stata") {
    stop(
      "`backend = \"stata\"` is no longer supported. Use `backend = \"native\"`.",
      call. = FALSE
    )
  }
  "native"
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
    if (sum(df[[group]] == 0, na.rm = TRUE) == 0L) {
      stop("The group variable must take value 0 for some observations.", call. = FALSE)
    }
    return(invisible(NULL))
  }

  for (nm in c(group, group_forward)) {
    bad <- !(df[[nm]] %in% c(-3, -1, 0, 1) | is.na(df[[nm]]))
    if (any(bad)) {
      stop(sprintf("The group variable `%s` takes values outside {-1,0,1,NA}.", nm), call. = FALSE)
    }
    if (sum(df[[nm]] == 0, na.rm = TRUE) == 0L) {
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

  if (opts$sieves && !has_cov) {
    stop("`sieves` requires covariates in the formula.", call. = FALSE)
  }

  if (opts$sieves && !has_cont) {
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
  if (all(!is.finite(scores))) {
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

        if (any(!is.finite(p_pred)) || any(!is.finite(y01_pred)) || any(!is.finite(y00_pred))) {
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
  d <- as.numeric(sub_df[[d_true_name]])
  g <- as.integer(sub_df$.g_binary)
  t <- as.integer(sub_df$.t_binary)

  idx11 <- g == 1 & t == 1
  idx10 <- g == 1 & t == 0
  idx01 <- g == 0 & t == 1
  idx00 <- g == 0 & t == 0

  den <- .safe_mean(d[idx11]) - .safe_mean(d[idx10])
  q_grid <- seq(0.05, 0.95, by = 0.05)

  if (!is.finite(den) || den == 0) {
    out <- rep(NA_real_, length(q_grid))
    names(out) <- sprintf("%.2f", q_grid)
    return(out)
  }

  mapped <- .counterfactual_quantile_map(
    y_target = y[idx10],
    y_d00 = y[idx00],
    y_d01 = y[idx01]
  )

  out <- vapply(
    q_grid,
    function(q) {
      q11 <- .safe_quantile(y[idx11], q)
      qcf <- .safe_quantile(mapped, q)
      (q11 - qcf) / den
    },
    numeric(1)
  )

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
  # Stata parity: bootstrap always runs with seed(1), independent of user seed.
  set.seed(1L)

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
      invalid_rep <- length(late_vals) != length(late_names) || any(!is.finite(late_vals))
    }

    if (!invalid_rep && !is.null(reps_lqte)) {
      lqte_vals <- as.numeric(est_b$lqte[lqte_names])
      invalid_rep <- length(lqte_vals) != length(lqte_names) || any(!is.finite(lqte_vals))
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

  if (!opts$nose) {
    n_reps <- bt$n_reps
    n_misreps <- bt$n_misreps
    share_failures <- bt$share_failures
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
#' @exportS3Method Rfuzzydid::summary
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

#' @title tidy.fuzzydid
#' @description Tidy extractor for `fuzzydid` objects.
#' @param x A fuzzydid object.
#' @param ... Unused.
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
