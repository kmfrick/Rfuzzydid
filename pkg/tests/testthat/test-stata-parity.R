make_parity_fixture <- function() {
  n_cell <- 120L

  make_cell <- function(g, t, ones) {
    data.frame(
      g = rep.int(g, n_cell),
      t = rep.int(t, n_cell),
      d = c(rep.int(1L, ones), rep.int(0L, n_cell - ones))
    )
  }

  df <- rbind(
    make_cell(0L, 0L, 24L),
    make_cell(0L, 1L, 48L),
    make_cell(1L, 0L, 36L),
    make_cell(1L, 1L, 96L)
  )

  df$id <- seq_len(nrow(df))
  df$y <- 1 + 0.5 * df$g + 0.4 * df$t + 2 * df$d + sin(df$id / 7)
  df$cl <- rep(seq_len(24L), each = 20L)

  df[, c("y", "g", "t", "d", "cl")]
}

test_that("Rfuzzydid matches frozen Stata parity goldens on core estimates and counts", {
  df <- make_parity_fixture()
  stata_out <- read_stata_parity_golden("stata-parity-core-golden.csv")

  invisible(capture.output({
    r_fit <- fuzzydid(
      data = df,
      formula = y ~ d,
      group = "g",
      time = "t",
      did = TRUE,
      tc = TRUE,
      cic = TRUE,
      breps = 5,
      backend = "native",
      seed = 1
    )
  }))

  r_estimates <- c(
    did = unname(r_fit$late$estimate[r_fit$late$estimator == "W_DID"]),
    tc = unname(r_fit$late$estimate[r_fit$late$estimator == "W_TC"]),
    cic = unname(r_fit$late$estimate[r_fit$late$estimator == "W_CIC"])
  )

  stata_estimates <- setNames(stata_out$estimate, stata_out$estimator)[c("did", "tc", "cic")]
  expect_equal(r_estimates, stata_estimates, tolerance = 1e-6)

  expected_counts <- unique(stata_out[c("n11", "n10", "n01", "n00")])
  expect_equal(nrow(expected_counts), 1)
  expect_equal(as.integer(r_fit$n11), as.integer(expected_counts$n11[[1]]))
  expect_equal(as.integer(r_fit$n10), as.integer(expected_counts$n10[[1]]))
  expect_equal(as.integer(r_fit$n01), as.integer(expected_counts$n01[[1]]))
  expect_equal(as.integer(r_fit$n00), as.integer(expected_counts$n00[[1]]))
})

test_that("Rfuzzydid matches frozen parity goldens on eqtest outputs", {
  df <- make_parity_fixture()
  eq_gold <- read_stata_parity_golden("stata-parity-eqtest-golden.csv")

  fit <- fuzzydid(
    data = df,
    formula = y ~ d,
    group = "g",
    time = "t",
    did = TRUE,
    tc = TRUE,
    cic = TRUE,
    eqtest = TRUE,
    breps = 5,
    backend = "native",
    seed = 1
  )

  expect_equal(fit$eqtest$contrast, eq_gold$contrast)
  expect_equal(fit$eqtest$estimate, eq_gold$estimate, tolerance = 1e-6)
  expect_equal(fit$eqtest$std.error, eq_gold$std.error, tolerance = 1e-6)
  expect_equal(fit$eqtest$conf.low, eq_gold$conf.low, tolerance = 1e-6)
  expect_equal(fit$eqtest$conf.high, eq_gold$conf.high, tolerance = 1e-6)
})

test_that("Rfuzzydid matches frozen parity goldens on lqte point estimates", {
  df <- make_parity_fixture()
  lqte_gold <- read_stata_parity_golden("stata-parity-lqte-golden.csv")

  fit <- fuzzydid(
    data = df,
    formula = y ~ d,
    group = "g",
    time = "t",
    lqte = TRUE,
    nose = TRUE,
    backend = "native",
    seed = 1
  )

  expect_equal(fit$lqte$quantile, lqte_gold$quantile, tolerance = 1e-12)
  expect_equal(fit$lqte$estimate, lqte_gold$estimate, tolerance = 1e-6)
})

test_that("Rfuzzydid matches frozen parity goldens on partial bounds", {
  set.seed(50321)
  n_cell <- 120
  make_cell <- function(g, t, prob_d) {
    d <- rbinom(n_cell, size = 1, prob = prob_d)
    y <- 1 + 0.5 * g + 0.4 * t + 2.0 * d + sin(seq_len(n_cell) / 7)
    data.frame(y = y, g = g, t = t, d = d)
  }

  df <- rbind(
    make_cell(g = 0, t = 0, prob_d = 0.20),
    make_cell(g = 0, t = 1, prob_d = 0.40),
    make_cell(g = 1, t = 0, prob_d = 0.30),
    make_cell(g = 1, t = 1, prob_d = 0.80)
  )

  partial_gold <- read_stata_parity_golden("stata-parity-partial-golden.csv")
  fit <- fuzzydid(
    data = df,
    formula = y ~ d,
    group = "g",
    time = "t",
    tc = TRUE,
    partial = TRUE,
    breps = 5,
    backend = "native",
    seed = 1
  )

  r_partial <- setNames(fit$late$estimate, fit$late$estimator)
  gold_partial <- setNames(partial_gold$estimate, partial_gold$estimator)
  expect_equal(
    unname(r_partial[names(gold_partial)]),
    unname(gold_partial),
    tolerance = 1e-6
  )
})

test_that("Rfuzzydid matches frozen parity goldens on clustered bootstrap diagnostics", {
  df <- make_parity_fixture()
  cluster_gold <- read_stata_parity_golden("stata-parity-cluster-golden.csv")

  fit <- fuzzydid(
    data = df,
    formula = y ~ d,
    group = "g",
    time = "t",
    did = TRUE,
    tc = TRUE,
    cic = TRUE,
    cluster = "cl",
    breps = 20,
    backend = "native",
    seed = 1
  )

  expect_equal(as.integer(fit$n_reps), as.integer(cluster_gold$n_reps[[1]]))
  expect_equal(as.integer(fit$n_misreps), as.integer(cluster_gold$n_misreps[[1]]))
  expect_equal(as.numeric(fit$share_failures), as.numeric(cluster_gold$share_failures[[1]]), tolerance = 1e-12)
})
