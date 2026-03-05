make_native_fixture <- function() {
  set.seed(42)
  n_cell <- 80L

  mk <- function(g, t, p) {
    d <- rbinom(n_cell, size = 1, prob = p)
    y <- 1 + 0.7 * g + 0.5 * t + 1.8 * d + rnorm(n_cell, sd = 0.2)
    data.frame(y = y, g = g, t = t, d = d)
  }

  rbind(
    mk(0, 0, 0.20),
    mk(0, 1, 0.35),
    mk(1, 0, 0.30),
    mk(1, 1, 0.70)
  )
}

make_covariate_fixture <- function() {
  df <- make_native_fixture()
  set.seed(123)
  df$x1 <- rnorm(nrow(df))
  df$x2 <- factor(ifelse(runif(nrow(df)) > 0.5, "a", "b"))
  df
}

make_group_forward_fixture <- function() {
  set.seed(7)
  n_cell <- 70L

  mk <- function(gb, gf, t, p, shift = 0) {
    d <- rbinom(n_cell, size = 1, prob = p)
    y <- 1 + 0.35 * t + 1.6 * d + shift + rnorm(n_cell, sd = 0.25)
    data.frame(y = y, d = d, gb = gb, gf = gf, t = t)
  }

  rbind(
    mk(0, 0, 0, 0.20, 0),
    mk(0, 0, 1, 0.30, 0),
    mk(0, 0, 2, 0.35, 0),
    mk(0, 1, 0, 0.30, 0.2),
    mk(1, 1, 1, 0.60, 0.8),
    mk(1, 0, 2, 0.75, 1.1)
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
    breps = 10,
    eqtest = TRUE,
    backend = "native",
    seed = 1
  )

  expect_s3_class(fit, "fuzzydid")
  expect_equal(sort(fit$late$estimator), sort(c("W_DID", "W_TC", "W_CIC")))
  expect_true(all(c("estimate", "std.error", "conf.low", "conf.high") %in% names(fit$late)))
  expect_true(!is.null(fit$eqtest))
  expect_true(!is.null(fit$matrices$b_LATE))
  expect_true(!is.null(fit$matrices$se_LATE))
  expect_true(!is.null(fit$matrices$ci_LATE))
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
    sieveorder = c(2, 2),
    nose = TRUE,
    backend = "native"
  )

  expect_true(all(c("W_DID", "W_TC") %in% fit_sieve$late$estimator))
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
      backend = "native"
    ),
    "without covariates"
  )

  expect_error(
    fuzzydid(
      data = make_group_forward_fixture(),
      formula = y ~ d,
      group = "gb",
      group_forward = "gf",
      time = "t",
      lqte = TRUE,
      backend = "native"
    ),
    "more than two periods"
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
      backend = "stata"
    ),
    "no longer supported"
  )
})
