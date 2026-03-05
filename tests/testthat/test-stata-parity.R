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

stata_core_golden <- c(
  did = 2.17430884241085,
  tc = 2.18553890987232,
  cic = 2.31982316267324
)

stata_cell_counts <- c(n11 = 120L, n10 = 120L, n01 = 120L, n00 = 120L)

stata_eqtest_golden <- data.frame(
  contrast = c("W_DID_W_TC", "W_DID_W_CIC", "W_TC_W_CIC"),
  estimate = c(-0.0112300674614629, -0.145514320262384, -0.134284252800921),
  std.error = c(0.161742439552579, 0.202172232202299, 0.0718583842239164),
  conf.low = c(-0.0435594151111732, -0.217510914771174, -0.179774021038072),
  conf.high = c(0.294748443505267, 0.240994490317917, -0.0170124474154521),
  stringsAsFactors = FALSE
)

stata_lqte_golden <- data.frame(
  quantile = seq(0.05, 0.95, by = 0.05),
  estimate = c(
    -1.42410977726763, -1.16020396591842, -0.155955670531315,
    1.16671849288358, 1.69883708320675, 1.46646296054033,
    1.53079877262529, 1.67443332467136, 1.88539747538913,
    0.697452557886513, 1.06766117780084, 1.63197830474743,
    2.08945772994696, 2.42192910221658, 2.77891151846807,
    2.92337346861517, 2.50163057321635, 1.60868067363183,
    0.973322738624446
  )
)

stata_partial_golden <- c(TC_inf = 0.01220596000000, TC_sup = 3.08177130000000)
stata_cluster_golden <- c(n_reps = 20L, n_misreps = 3L, share_failures = 0.15)

test_that("Rfuzzydid matches frozen Stata parity goldens on core estimates and counts", {
  df <- make_parity_fixture()

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

  expect_equal(r_estimates, stata_core_golden, tolerance = 1e-6)
  expect_equal(as.integer(r_fit$n11), stata_cell_counts[["n11"]])
  expect_equal(as.integer(r_fit$n10), stata_cell_counts[["n10"]])
  expect_equal(as.integer(r_fit$n01), stata_cell_counts[["n01"]])
  expect_equal(as.integer(r_fit$n00), stata_cell_counts[["n00"]])
})

test_that("Rfuzzydid matches frozen parity goldens on eqtest outputs", {
  df <- make_parity_fixture()

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

  expect_equal(fit$eqtest$contrast, stata_eqtest_golden$contrast)
  expect_equal(fit$eqtest$estimate, stata_eqtest_golden$estimate, tolerance = 1e-6)
  expect_equal(fit$eqtest$std.error, stata_eqtest_golden$std.error, tolerance = 1e-6)
  expect_equal(fit$eqtest$conf.low, stata_eqtest_golden$conf.low, tolerance = 1e-6)
  expect_equal(fit$eqtest$conf.high, stata_eqtest_golden$conf.high, tolerance = 1e-6)
})

test_that("Rfuzzydid matches frozen parity goldens on lqte point estimates", {
  df <- make_parity_fixture()

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

  expect_equal(fit$lqte$quantile, stata_lqte_golden$quantile, tolerance = 1e-12)
  expect_equal(fit$lqte$estimate, stata_lqte_golden$estimate, tolerance = 1e-6)
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
  expect_equal(
    unname(r_partial[names(stata_partial_golden)]),
    unname(stata_partial_golden),
    tolerance = 1e-6
  )
})

test_that("Rfuzzydid matches frozen parity goldens on clustered bootstrap diagnostics", {
  df <- make_parity_fixture()

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

  expect_equal(as.integer(fit$n_reps), as.integer(stata_cluster_golden[["n_reps"]]))
  expect_equal(as.integer(fit$n_misreps), as.integer(stata_cluster_golden[["n_misreps"]]))
  expect_equal(as.numeric(fit$share_failures), as.numeric(stata_cluster_golden[["share_failures"]]), tolerance = 1e-12)
})
