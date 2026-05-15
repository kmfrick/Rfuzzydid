make_bootstrap_seed_fixture <- function() {
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

test_that("bootstrap seed is optional and user controlled", {
  df <- make_bootstrap_seed_fixture()

  fit_default_seed <- fuzzydid(
    data = df,
    formula = y ~ d,
    group = "g",
    time = "t",
    did = TRUE,
    breps = 10,
    backend = "native"
  )

  set.seed(123)
  fit_seed_1_a <- fuzzydid(
    data = df,
    formula = y ~ d,
    group = "g",
    time = "t",
    did = TRUE,
    breps = 10,
    backend = "native",
    seed = 1
  )

  set.seed(999)
  fit_seed_1_b <- fuzzydid(
    data = df,
    formula = y ~ d,
    group = "g",
    time = "t",
    did = TRUE,
    breps = 10,
    backend = "native",
    seed = 1
  )

  expect_s3_class(fit_default_seed, "fuzzydid")
  expect_equal(fit_seed_1_a$late$std.error, fit_seed_1_b$late$std.error, tolerance = 1e-12)
  expect_equal(fit_seed_1_a$late$conf.low, fit_seed_1_b$late$conf.low, tolerance = 1e-12)
  expect_equal(fit_seed_1_a$late$conf.high, fit_seed_1_b$late$conf.high, tolerance = 1e-12)
})
