read_stata_parity_golden <- function(filename) {
  candidates <- c(
    file.path("pkg", "tests", "testthat", filename),
    file.path("tests", "testthat", filename),
    file.path("..", "tests", "testthat", filename),
    file.path("..", "..", "tests", "testthat", filename)
  )

  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(utils::read.csv(candidate, stringsAsFactors = FALSE))
    }
  }

  stop(sprintf("Parity golden file not found: %s", filename), call. = FALSE)
}
