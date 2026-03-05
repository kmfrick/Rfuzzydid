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

test_that("native backend enforces unsupported features clearly", {
  df <- make_native_fixture()

  expect_error(
    fuzzydid(
      data = df,
      formula = y ~ d,
      group = "g",
      time = "t",
      lqte = TRUE,
      backend = "native"
    ),
    "not available in native backend"
  )

  expect_error(
    fuzzydid(
      data = df,
      formula = y ~ d,
      group = "g",
      time = "t",
      did = TRUE,
      modelx = c("ols", "ols"),
      backend = "native"
    ),
    "not available in native backend"
  )

  expect_error(
    fuzzydid(
      data = df,
      formula = y ~ d,
      group = "g",
      time = "t",
      did = TRUE,
      sieves = TRUE,
      backend = "native"
    ),
    "not available in native backend"
  )
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
