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
  expect_equal(sort(fit$late$estimator), sort(c("W_DID", "W_TC", "W_CIC")))
  expect_true("estimate" %in% names(fit$late))
  expect_true(!is.null(fit$matrices$b_LATE))
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

  out <- Rfuzzydid:::.calc_boot_summary(reps)
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

  expect_true(!is.null(fit$lqte))
  expect_equal(nrow(fit$lqte), 19L)
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

  expect_equal(fit$late$estimator, c("TC_inf", "TC_sup"))
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

  expect_true(!is.null(fit_cv_1$sieveorder_selected))
  expect_equal(as.integer(fit_cv_1$sieveorder_selected), as.integer(fit_cv_2$sieveorder_selected))
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

  expect_true(is.logical(fit$tagobs))
  expect_equal(length(fit$tagobs), nrow(df))
  expect_false("tagobs" %in% names(df))
})

test_that("stata backend option is rejected explicitly", {
  df <- make_native_fixture()

  expect_error(
    fuzzydid(
      data = df,
      formula = y ~ d,
      group = "g",
      time = "t",
      did = TRUE,
      nose = TRUE,
      backend = "stata"
    ),
    "no longer supported"
  )
})
