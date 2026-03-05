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
#' @param modelx Reserved for native covariate-adjusted methods (`ols`,
#'   `logit`, `probit`); currently unsupported.
#' @param sieves Reserved for native sieve estimation; currently unsupported.
#' @param sieveorder Optional length-2 numeric vector for sieve orders.
#' @param tagobs Logical; return logical mask of observations used.
#' @param backend One of `"auto"` or `"native"`. `"stata"` is no longer
#'   supported and errors.
#' @param seed Bootstrap/random seed.
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
    allowed <- c("ols", "logit", "probit")
    if (!all(opts$modelx %in% allowed)) {
      stop("`modelx` entries must be one of: ols, logit, probit.", call. = FALSE)
    }
  }

  if (!is.null(opts$sieveorder)) {
    if (!is.numeric(opts$sieveorder) || length(opts$sieveorder) != 2L) {
      stop("`sieveorder` must be a numeric vector of length 2.", call. = FALSE)
    }
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
  role_vars <- c(y_name, d_name, group, time, group_forward, covariates, cluster)
  role_vars <- role_vars[nzchar(role_vars)]
  role_vars <- unique(role_vars)

  missing_vars <- setdiff(role_vars, names(data))
  if (length(missing_vars) > 0L) {
    stop(
      sprintf("Required variable(s) not found in `data`: %s", paste(missing_vars, collapse = ", ")),
      call. = FALSE
    )
  }

  df <- data
  complete_mask <- stats::complete.cases(df[, role_vars, drop = FALSE])
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
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

.validate_native_design <- function(df, group, time) {
  g_vals <- sort(unique(df[[group]]))
  t_vals <- sort(unique(df[[time]]))

  if (length(g_vals) != 2L || length(t_vals) != 2L) {
    stop(
      "Native backend currently supports exactly two group values and two time periods.",
      call. = FALSE
    )
  }

  list(g0 = g_vals[[1L]], g1 = g_vals[[2L]], t0 = t_vals[[1L]], t1 = t_vals[[2L]])
}

.recategorize_treatment <- function(d, newcateg) {
  if (is.null(newcateg)) return(as.numeric(d))
  d <- as.numeric(d)
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

.native_point_estimates <- function(df, names, opts) {
  y <- df[[names$y]]
  g <- df[[names$g]]
  t <- df[[names$t]]
  d_true <- df[[names$d_true]]
  d_tc <- df[[names$d_tc]]

  design <- .validate_native_design(df, names$g, names$t)
  g0 <- design$g0
  g1 <- design$g1
  t0 <- design$t0
  t1 <- design$t1

  m_y11 <- .safe_mean(y[g == g1 & t == t1])
  m_y10 <- .safe_mean(y[g == g1 & t == t0])
  m_y01 <- .safe_mean(y[g == g0 & t == t1])
  m_y00 <- .safe_mean(y[g == g0 & t == t0])

  m_d11 <- .safe_mean(d_true[g == g1 & t == t1])
  m_d10 <- .safe_mean(d_true[g == g1 & t == t0])
  m_d01 <- .safe_mean(d_true[g == g0 & t == t1])
  m_d00 <- .safe_mean(d_true[g == g0 & t == t0])

  did_num <- m_y11 - m_y10 - (m_y01 - m_y00)
  did_den <- m_d11 - m_d10 - (m_d01 - m_d00)
  did_est <- did_num / did_den

  out <- list()
  if (opts$did) {
    if (opts$numerator) {
      out$DID_num <- did_num
    } else {
      out$W_DID <- did_est
    }
  }

  if (opts$tc) {
    special_case <- length(unique(d_tc[g == g0])) == 1L
    if (special_case) {
      if (opts$numerator) {
        out$TC_num <- did_num
      } else {
        out$W_TC <- did_est
      }
    } else {
      levels_d <- sort(unique(d_tc))
      build <- 0
      valid <- TRUE
      for (lv in levels_d) {
        p_d10 <- mean(d_tc[g == g1 & t == t0] == lv)
        y01_d <- .safe_mean(y[g == g0 & t == t1 & d_tc == lv])
        y00_d <- .safe_mean(y[g == g0 & t == t0 & d_tc == lv])
        if (is.na(y01_d) || is.na(y00_d)) {
          valid <- FALSE
          break
        }
        build <- build + p_d10 * (y01_d - y00_d)
      }
      tc_num <- if (valid) m_y11 - m_y10 - build else NA_real_
      tc_den <- m_d11 - m_d10
      tc_est <- tc_num / tc_den
      if (opts$numerator) {
        out$TC_num <- tc_num
      } else {
        out$W_TC <- tc_est
      }
    }
  }

  if (opts$cic) {
    special_case <- length(unique(d_tc[g == g0])) == 1L
    levels_10 <- sort(unique(d_tc[g == g1 & t == t0]))
    q_vals <- numeric(0)
    ctrl_level <- if (special_case) unique(d_tc[g == g0])[1L] else NA_real_

    for (lv in levels_10) {
      y10 <- y[g == g1 & t == t0 & d_tc == lv]
      if (special_case) {
        y00 <- y[g == g0 & t == t0 & d_tc == ctrl_level]
        y01 <- y[g == g0 & t == t1 & d_tc == ctrl_level]
      } else {
        y00 <- y[g == g0 & t == t0 & d_tc == lv]
        y01 <- y[g == g0 & t == t1 & d_tc == lv]
      }

      mapped <- .counterfactual_quantile_map(y_target = y10, y_d00 = y00, y_d01 = y01)
      q_vals <- c(q_vals, mapped)
    }

    cic_num <- m_y11 - .safe_mean(q_vals)
    cic_den <- m_d11 - m_d10
    cic_est <- cic_num / cic_den
    if (opts$numerator) {
      out$CIC_num <- cic_num
    } else {
      out$W_CIC <- cic_est
    }
  }

  unlist(out)
}

.bootstrap_native <- function(df, names, opts) {
  if (opts$nose) return(NULL)
  set.seed(as.integer(opts$seed))

  point <- .native_point_estimates(df, names, opts)
  stat_names <- names(point)
  reps <- matrix(NA_real_, nrow = as.integer(opts$breps), ncol = length(stat_names))
  colnames(reps) <- stat_names

  n <- nrow(df)
  for (b in seq_len(as.integer(opts$breps))) {
    if (is.null(opts$cluster)) {
      idx <- sample.int(n, size = n, replace = TRUE)
    } else {
      cl <- unique(df[[opts$cluster]])
      sampled_cl <- sample(cl, size = length(cl), replace = TRUE)
      idx <- unlist(lapply(sampled_cl, function(cc) which(df[[opts$cluster]] == cc)), use.names = FALSE)
    }
    est_b <- .native_point_estimates(df[idx, , drop = FALSE], names, opts)
    reps[b, ] <- as.numeric(est_b[stat_names])
  }

  list(point = point, reps = reps)
}

.pairwise_eqtest <- function(estimates_named) {
  nm <- names(estimates_named)
  if (length(nm) < 2L) return(NULL)
  out <- list()
  for (i in seq_len(length(nm) - 1L)) {
    for (j in seq.int(i + 1L, length(nm))) {
      key <- paste0(sub("^W_", "", nm[[i]]), "_", sub("^W_", "", nm[[j]]))
      out[[key]] <- estimates_named[[i]] - estimates_named[[j]]
    }
  }
  unlist(out)
}

.run_native_backend <- function(prepared, opts) {
  if (length(prepared$covariates) > 0L) {
    stop("Native backend does not currently support covariates.", call. = FALSE)
  }
  if (!is.null(prepared$group_forward)) {
    stop("Native backend does not support `group_forward`.", call. = FALSE)
  }
  if (!is.null(opts$modelx)) {
    stop("Requested options are not available in native backend (`modelx`).", call. = FALSE)
  }
  if (opts$sieves) {
    stop("Requested options are not available in native backend (`sieves`).", call. = FALSE)
  }
  if (!is.null(opts$sieveorder)) {
    stop("Requested options are not available in native backend (`sieveorder`).", call. = FALSE)
  }
  if (opts$lqte) {
    stop("Requested options are not available in native backend (`lqte`).", call. = FALSE)
  }
  if (opts$partial) {
    stop("Requested options are not available in native backend (`partial`).", call. = FALSE)
  }

  df <- prepared$used
  d_tc <- .recategorize_treatment(df[[prepared$d_name]], opts$newcateg)
  df$d_true_native <- as.numeric(df[[prepared$d_name]])
  df$d_tc_native <- d_tc
  df$g_native <- as.numeric(df[[prepared$group]])
  df$t_native <- as.numeric(df[[prepared$time]])
  df$y_native <- as.numeric(df[[prepared$y_name]])

  names_map <- list(
    y = "y_native",
    g = "g_native",
    t = "t_native",
    d_true = "d_true_native",
    d_tc = "d_tc_native"
  )

  point <- .native_point_estimates(df, names_map, opts)
  bt <- .bootstrap_native(df, names_map, opts)
  if (!is.null(bt)) {
    point <- bt$point
    reps <- bt$reps
  } else {
    reps <- NULL
  }

  late_terms <- names(point)
  late <- data.frame(
    estimator = late_terms,
    estimate = as.numeric(point),
    std.error = NA_real_,
    conf.low = NA_real_,
    conf.high = NA_real_,
    stringsAsFactors = FALSE
  )

  matrices <- list(
    b_LATE = .to_column_matrix(late$estimate, late$estimator)
  )

  if (!opts$nose) {
    se <- apply(reps, 2, stats::sd, na.rm = TRUE)
    ci <- t(apply(reps, 2, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE))
    late$std.error <- as.numeric(se)
    late$conf.low <- as.numeric(ci[, 1])
    late$conf.high <- as.numeric(ci[, 2])
    matrices$se_LATE <- .to_column_matrix(late$std.error, late$estimator)
    matrices$ci_LATE <- .to_ci_matrix(late$conf.low, late$conf.high, late$estimator)
  }

  eqtest <- NULL
  if (opts$eqtest && !opts$numerator) {
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

      if (!opts$nose) {
        eq_reps <- t(apply(reps, 1, function(row) .pairwise_eqtest(stats::setNames(as.numeric(row), colnames(reps)))))
        if (!is.null(eq_reps)) {
          se_eq <- apply(eq_reps, 2, stats::sd, na.rm = TRUE)
          ci_eq <- t(apply(eq_reps, 2, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE))
          eqtest$std.error <- as.numeric(se_eq)
          eqtest$conf.low <- as.numeric(ci_eq[, 1])
          eqtest$conf.high <- as.numeric(ci_eq[, 2])
        }
      }

      matrices$b_LATE_eqtest <- .to_column_matrix(eqtest$estimate, eqtest$contrast)
      if (!opts$nose) {
        matrices$se_LATE_eqtest <- .to_column_matrix(eqtest$std.error, eqtest$contrast)
        matrices$ci_LATE_eqtest <- .to_ci_matrix(eqtest$conf.low, eqtest$conf.high, eqtest$contrast)
      }
    }
  }

  # Binary-cell counts where meaningful.
  g_vals <- sort(unique(df$g_native))
  t_vals <- sort(unique(df$t_native))
  n11 <- NA_integer_
  n10 <- NA_integer_
  n01 <- NA_integer_
  n00 <- NA_integer_
  if (length(g_vals) == 2L && length(t_vals) == 2L) {
    g0 <- g_vals[[1L]]
    g1 <- g_vals[[2L]]
    t0 <- t_vals[[1L]]
    t1 <- t_vals[[2L]]
    n11 <- sum(df$g_native == g1 & df$t_native == t1)
    n10 <- sum(df$g_native == g1 & df$t_native == t0)
    n01 <- sum(df$g_native == g0 & df$t_native == t1)
    n00 <- sum(df$g_native == g0 & df$t_native == t0)
  }

  tag_mask <- if (opts$tagobs) prepared$mask else NULL

  list(
    late = late,
    eqtest = eqtest,
    lqte = NULL,
    matrices = matrices,
    tagobs = tag_mask,
    n = nrow(df),
    n11 = n11,
    n10 = n10,
    n01 = n01,
    n00 = n00
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
    stringsAsFactors = FALSE
  )
}

#' @title tidy
#' @description Generic tidy method.
#' @param x Object to tidy.
#' @param ... Additional arguments.
#' @export
tidy <- function(x, ...) {
  UseMethod("tidy")
}

#' @title glance
#' @description Generic glance method.
#' @param x Object to summarize.
#' @param ... Additional arguments.
#' @export
glance <- function(x, ...) {
  UseMethod("glance")
}
