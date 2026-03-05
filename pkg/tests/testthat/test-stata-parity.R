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

  df[, c("y", "g", "t", "d")]
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
