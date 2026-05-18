#' srr_stats_native_tests
#'
#' @srrstats {G5.2} This file tests representative errors and warnings for the
#' native API.
#' @srrstats {RE7.3} This file tests the model-like extractors and plot method
#' required by the regression standards.
#' @noRd
NULL

make_native_fixture <- function() {
  n_cell <- 80L

  mk <- function(g, t, ones) {
    d <- c(rep.int(1L, ones), rep.int(0L, n_cell - ones))
    id <- seq_len(n_cell)
    y <- 1 + 0.7 * g + 0.5 * t + 1.8 * d + 0.2 * sin(id / 5)
    data.frame(y = y, g = g, t = t, d = d)
  }

  rbind(
    mk(0L, 0L, 16L),
    mk(0L, 1L, 28L),
    mk(1L, 0L, 24L),
    mk(1L, 1L, 56L)
  )
}

make_covariate_fixture <- function() {
  df <- make_native_fixture()
  n <- nrow(df)
  idx <- seq_len(n)
  df$x1 <- sin(idx / 11) + cos(idx / 17)
  df$x2 <- factor(ifelse((idx %% 2L) == 0L, "a", "b"))
  df
}

make_group_forward_fixture <- function() {
  n_cell <- 70L

  mk <- function(gb, gf, t, ones, shift = 0) {
    d <- c(rep.int(1L, ones), rep.int(0L, n_cell - ones))
    id <- seq_len(n_cell)
    y <- 1 + 0.35 * t + 1.6 * d + shift + 0.25 * cos(id / 6)
    data.frame(y = y, d = d, gb = gb, gf = gf, t = t)
  }

  rbind(
    mk(0L, 0L, 0L, 14L, 0),
    mk(0L, 0L, 1L, 21L, 0),
    mk(0L, 0L, 2L, 25L, 0),
    mk(0L, 1L, 0L, 21L, 0.2),
    mk(1L, 1L, 1L, 42L, 0.8),
    mk(1L, 0L, 2L, 53L, 1.1)
  )
}

test_that("formula API works on native backend", {
  df <- make_native_fixture()

  fit <- fuzzydid(
    data = df,
    formula = y ~ d,
    group = "g",
    time = "t",
    did = TRUE,
    tc = TRUE,
    cic = TRUE,
    nose = TRUE,
    backend = "native",
    seed = 1
  )

  expect_s3_class(fit, "fuzzydid")
  expect_identical(sort(fit$late$estimator), sort(c("W_DID", "W_TC", "W_CIC")))
  expect_true("estimate" %in% names(fit$late))
  expect_false(is.null(fit$matrices$b_LATE))
})

test_that("summary, tidy, and glance methods expose stable outputs", {
  df <- make_native_fixture()

  fit <- fuzzydid(
    data = df,
    formula = y ~ d,
    group = "g",
    time = "t",
    did = TRUE,
    tc = TRUE,
    nose = TRUE,
    backend = "native"
  )

  invisible(capture.output(out <- summary(fit)))
  expect_identical(out, fit)

  tidy_out <- generics::tidy(fit)
  expect_s3_class(tidy_out, "data.frame")
  expect_true(all(c("component", "term", "estimate") %in% names(tidy_out)))
  expect_true(all(tidy_out$component %in% c("late", "eqtest", "lqte")))

  glance_out <- generics::glance(fit)
  expect_s3_class(glance_out, "data.frame")
  expect_identical(nrow(glance_out), 1L)
  expect_true(all(c("backend", "Num.Obs.", "N.reps") %in% names(glance_out)))

  expect_named(coef(fit), c("W_DID", "W_TC"))
  expect_identical(nobs(fit), nrow(df))
  expect_identical(formula(fit), y ~ d)
  expect_error(vcov(fit), "Bootstrap covariance is unavailable")
  expect_output(print(fit), "fuzzydid fit")
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(fit))
})

test_that("bootstrap fit exposes coefficient intervals and covariance", {
  df <- make_native_fixture()

  fit <- fuzzydid(
    data = df,
    formula = y ~ d,
    group = "g",
    time = "t",
    did = TRUE,
    tc = TRUE,
    breps = 10,
    seed = 1,
    backend = "native"
  )

  ci <- confint(fit)
  expect_true(is.matrix(ci))
  expect_identical(rownames(ci), names(coef(fit)))
  expect_error(confint(fit, level = 0.9), "95%")

  cov <- vcov(fit)
  expect_true(is.matrix(cov))
  expect_identical(rownames(cov), names(coef(fit)))
  expect_identical(colnames(cov), names(coef(fit)))

  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot_out <- plot(fit))
  expect_identical(plot_out, fit)
})

test_that("treatment is explicit and robust to RHS ordering", {
  df <- make_covariate_fixture()

  fit_explicit_1 <- fuzzydid(
    data = df,
    formula = y ~ d + x1,
    treatment = "d",
    group = "g",
    time = "t",
    did = TRUE,
    nose = TRUE,
    backend = "native"
  )

  fit_explicit_2 <- fuzzydid(
    data = df,
    formula = y ~ x1 + d,
    treatment = "d",
    group = "g",
    time = "t",
    did = TRUE,
    nose = TRUE,
    backend = "native"
  )

  expect_equal(fit_explicit_1$late$estimate, fit_explicit_2$late$estimate, tolerance = 1e-12)

  expect_error(
    fuzzydid(
      data = df,
      formula = y ~ x1 + x2,
      group = "g",
      time = "t",
      did = TRUE,
      nose = TRUE,
      backend = "native"
    ),
    "Unable to infer treatment"
  )
})


test_that("bootstrap summary drops Stata sentinel reps", {
  reps <- matrix(
    c(
      1, 2,
      1000000000000000, -1000000000000000,
      3, 4
    ),
    ncol = 2,
    byrow = TRUE
  )

  calc_boot_summary <- getFromNamespace(".calc_boot_summary", "Rfuzzydid")
  out <- calc_boot_summary(reps)
  expect_equal(out$se, c(stats::sd(c(1, 3)), stats::sd(c(2, 4))), tolerance = 1e-12)
})

test_that("group_forward supports multi-period DID/TC", {
  df <- make_group_forward_fixture()

  fit <- fuzzydid(
    data = df,
    formula = y ~ d,
    group = "gb",
    group_forward = "gf",
    time = "t",
    did = TRUE,
    tc = TRUE,
    nose = TRUE,
    backend = "native"
  )

  expect_true(all(c("W_DID", "W_TC") %in% fit$late$estimator))
})

test_that("native backend supports lqte under binary two-period design", {
  df <- make_native_fixture()

  fit <- fuzzydid(
    data = df,
    formula = y ~ d,
    group = "g",
    time = "t",
    lqte = TRUE,
    nose = TRUE,
    backend = "native"
  )

  expect_false(is.null(fit$lqte))
  expect_identical(nrow(fit$lqte), 19L)
})

test_that("partial returns TC bounds under valid design", {
  df <- make_native_fixture()

  fit <- fuzzydid(
    data = df,
    formula = y ~ d,
    group = "g",
    time = "t",
    tc = TRUE,
    partial = TRUE,
    nose = TRUE,
    backend = "native"
  )

  expect_identical(fit$late$estimator, c("TC_inf", "TC_sup"))
})

test_that("eqtest supports point contrasts when standard errors are skipped", {
  df <- make_native_fixture()

  fit <- fuzzydid(
    data = df,
    formula = y ~ d,
    group = "g",
    time = "t",
    did = TRUE,
    tc = TRUE,
    cic = TRUE,
    eqtest = TRUE,
    nose = TRUE,
    backend = "native"
  )

  expect_s3_class(fit$eqtest, "data.frame")
  expect_identical(
    fit$eqtest$contrast,
    c("W_DID_W_TC", "W_DID_W_CIC", "W_TC_W_CIC")
  )
  expect_true(all(is.na(fit$eqtest$std.error)))
})

test_that("sieves is a no-op without formula covariates", {
  df <- make_native_fixture()

  fit_plain <- fuzzydid(
    data = df,
    formula = y ~ d,
    group = "g",
    time = "t",
    did = TRUE,
    tc = TRUE,
    cic = TRUE,
    lqte = TRUE,
    nose = TRUE,
    backend = "native"
  )

  fit_sieves <- fuzzydid(
    data = df,
    formula = y ~ d,
    group = "g",
    time = "t",
    did = TRUE,
    tc = TRUE,
    cic = TRUE,
    lqte = TRUE,
    sieves = TRUE,
    nose = TRUE,
    backend = "native"
  )

  expect_equal(fit_sieves$late$estimate, fit_plain$late$estimate)
  expect_equal(fit_sieves$lqte$estimate, fit_plain$lqte$estimate)
  expect_null(fit_sieves$sieveorder_selected)
})

test_that("covariates with modelx and sieves are supported for DID/TC", {
  df <- make_covariate_fixture()

  fit_modelx <- fuzzydid(
    data = df,
    formula = y ~ d + x1 + x2,
    group = "g",
    time = "t",
    did = TRUE,
    tc = TRUE,
    modelx = c("ols", "ols"),
    nose = TRUE,
    backend = "native"
  )

  expect_true(all(c("W_DID", "W_TC") %in% fit_modelx$late$estimator))

  fit_sieve <- fuzzydid(
    data = df,
    formula = y ~ d + x1,
    group = "g",
    time = "t",
    did = TRUE,
    tc = TRUE,
    sieves = TRUE,
    nose = TRUE,
    backend = "native"
  )

  expect_true(all(c("W_DID", "W_TC") %in% fit_sieve$late$estimator))
})

test_that("sieveorder defaults to deterministic CV selection and supports legacy length-2", {
  df <- make_covariate_fixture()

  fit_cv_1 <- fuzzydid(
    data = df,
    formula = y ~ d + x1,
    group = "g",
    time = "t",
    did = TRUE,
    tc = TRUE,
    sieves = TRUE,
    nose = TRUE,
    backend = "native"
  )

  fit_cv_2 <- fuzzydid(
    data = df,
    formula = y ~ d + x1,
    group = "g",
    time = "t",
    did = TRUE,
    tc = TRUE,
    sieves = TRUE,
    nose = TRUE,
    backend = "native"
  )

  expect_false(is.null(fit_cv_1$sieveorder_selected))
  expect_identical(as.integer(fit_cv_1$sieveorder_selected), as.integer(fit_cv_2$sieveorder_selected))
  expect_true(all(as.integer(fit_cv_1$sieveorder_selected) >= 2L))

  expect_warning(
    fuzzydid(
      data = df,
      formula = y ~ d + x1,
      group = "g",
      time = "t",
      did = TRUE,
      tc = TRUE,
      sieves = TRUE,
      sieveorder = c(2, 2),
      nose = TRUE,
      backend = "native"
    ),
    "deprecated"
  )
})

test_that("strict restrictions mirror Stata-style constraints", {
  df <- make_covariate_fixture()

  expect_error(
    fuzzydid(
      data = df,
      formula = y ~ d + x1,
      group = "g",
      time = "t",
      tc = TRUE,
      partial = TRUE,
      nose = TRUE,
      backend = "native"
    ),
    "without covariates"
  )

  expect_error(
    fuzzydid(
      data = df,
      formula = y ~ d + x1,
      group = "g",
      time = "t",
      did = TRUE,
      tc = TRUE,
      sieves = TRUE,
      sieveorder = 1,
      nose = TRUE,
      backend = "native"
    ),
    ">= 2"
  )

  expect_error(
    fuzzydid(
      data = df,
      formula = y ~ d + x1,
      group = "g",
      time = "t",
      did = TRUE,
      tc = TRUE,
      sieves = TRUE,
      sieveorder = 200,
      nose = TRUE,
      backend = "native"
    ),
    "min\\(4800, floor\\(n/5\\)\\)"
  )

  expect_error(
    fuzzydid(
      data = make_group_forward_fixture(),
      formula = y ~ d,
      group = "gb",
      group_forward = "gf",
      time = "t",
      lqte = TRUE,
      nose = TRUE,
      backend = "native"
    ),
    "more than two periods"
  )

  expect_error(
    fuzzydid(
      data = make_native_fixture(),
      formula = y ~ d,
      group = "g",
      time = "t",
      did = TRUE,
      breps = 2.5,
      backend = "native"
    ),
    "integer scalar >= 2"
  )

  expect_error(
    fuzzydid(
      data = make_native_fixture(),
      formula = y ~ d,
      group = "g",
      time = "t",
      did = 1,
      nose = TRUE,
      backend = "native"
    ),
    "`did` must be TRUE or FALSE",
    fixed = TRUE
  )
})

test_that("explicit breps is ignored when standard errors are skipped", {
  df <- make_native_fixture()

  fit <- fuzzydid(
    data = df,
    formula = y ~ d,
    group = "g",
    time = "t",
    did = TRUE,
    breps = 2L,
    nose = TRUE,
    backend = "native"
  )

  expect_s3_class(fit, "fuzzydid")
  expect_true(is.na(fit$late$std.error))
})

test_that("group_forward value check message includes -3", {
  df <- make_group_forward_fixture()
  df$gb[1] <- 2

  expect_error(
    fuzzydid(
      data = df,
      formula = y ~ d,
      group = "gb",
      group_forward = "gf",
      time = "t",
      did = TRUE,
      nose = TRUE,
      backend = "native"
    ),
    "\\{-3,-1,0,1,NA\\}"
  )
})

test_that("tagobs mask is returned without mutating input data", {
  df <- make_native_fixture()
  df$z <- NA_real_
  df$z[1:10] <- 1

  fit <- fuzzydid(
    data = df,
    formula = y ~ d,
    group = "g",
    time = "t",
    did = TRUE,
    nose = TRUE,
    tagobs = TRUE,
    backend = "native"
  )

  expect_type(fit$tagobs, "logical")
  expect_length(fit$tagobs, nrow(df))
  expect_false("tagobs" %in% names(df))
})

test_that("analysis roles reject unsupported and non-finite inputs", {
  df <- make_native_fixture()

  bad_outcome <- df
  bad_outcome$y <- as.character(bad_outcome$y)
  expect_error(
    fuzzydid(bad_outcome, y ~ d, group = "g", time = "t", did = TRUE, nose = TRUE),
    "`y` must be a numeric vector",
    fixed = TRUE
  )

  bad_treatment <- df
  bad_treatment$d <- factor(bad_treatment$d)
  expect_error(
    fuzzydid(bad_treatment, y ~ d, group = "g", time = "t", did = TRUE, nose = TRUE),
    "`d` must be a numeric vector",
    fixed = TRUE
  )

  bad_inf <- df
  bad_inf$y[1] <- Inf
  expect_error(
    fuzzydid(bad_inf, y ~ d, group = "g", time = "t", did = TRUE, nose = TRUE),
    "`y` must not contain Inf or -Inf",
    fixed = TRUE
  )

  bad_covariate <- make_covariate_fixture()
  bad_covariate$x_bad <- I(as.list(seq_len(nrow(bad_covariate))))
  expect_error(
    fuzzydid(
      bad_covariate,
      y ~ d + x_bad,
      group = "g",
      time = "t",
      did = TRUE,
      nose = TRUE
    ),
    "Covariate `x_bad` must be a numeric, factor, character, or logical vector",
    fixed = TRUE
  )
})

test_that("missing values are complete-case filtered and exposed by tagobs", {
  df <- make_native_fixture()
  df$y[1] <- NA_real_
  df$d[2] <- NaN

  fit <- fuzzydid(
    data = df,
    formula = y ~ d,
    group = "g",
    time = "t",
    did = TRUE,
    nose = TRUE,
    tagobs = TRUE,
    backend = "native"
  )

  expect_identical(fit$n, nrow(df) - 2L)
  expect_false(fit$tagobs[1])
  expect_false(fit$tagobs[2])
})

test_that("exact covariate relationships remain computable", {
  df <- make_covariate_fixture()
  df$x_exact <- df$y
  df$x_duplicate <- df$x1

  fit <- fuzzydid(
    data = df,
    formula = y ~ d + x1 + x_duplicate + x_exact,
    group = "g",
    time = "t",
    did = TRUE,
    tc = TRUE,
    nose = TRUE,
    backend = "native"
  )

  expect_true(all(is.finite(fit$late$estimate)))
})
