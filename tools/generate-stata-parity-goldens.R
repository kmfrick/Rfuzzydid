#!/usr/bin/env Rscript

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

stata_bin <- Sys.which("stata")
if (!nzchar(stata_bin)) {
  stop("Could not find `stata` on PATH.", call. = FALSE)
}

tmp_dir <- tempfile("stata-parity-")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

fixture_path <- file.path(tmp_dir, "parity.csv")
utils::write.csv(make_parity_fixture(), fixture_path, row.names = FALSE)

do_path <- file.path(tmp_dir, "parity.do")
fixture_path_stata <- gsub("\\\\", "/", normalizePath(fixture_path, winslash = "/"))

do_lines <- c(
  "clear all",
  "set more off",
  "capture which fuzzydid",
  "if _rc {",
  "  di as error \"CODEX_ERROR|missing_fuzzydid\"",
  "  exit 111",
  "}",
  sprintf("import delimited using \"%s\", clear", fixture_path_stata),
  "quietly fuzzydid y g t d, did tc cic nose",
  "display \"CODEX_RESULT|core|did|\" %21.15f el(e(b_LATE),1,1)",
  "display \"CODEX_RESULT|core|tc|\" %21.15f el(e(b_LATE),2,1)",
  "display \"CODEX_RESULT|core|cic|\" %21.15f el(e(b_LATE),3,1)",
  "quietly fuzzydid y g t d, tc partial nose",
  "display \"CODEX_RESULT|partial|TC_inf|\" %21.15f el(e(b_LATE),1,1)",
  "display \"CODEX_RESULT|partial|TC_sup|\" %21.15f el(e(b_LATE),2,1)",
  "quietly fuzzydid y g t d, lqte nose",
  "forvalues i=1/19 {",
  "  local q = `i' * 5",
  "  display \"CODEX_RESULT|lqte|q_`q'|\" %21.15f el(e(b_LQTE),`i',1)",
  "}",
  "exit, clear"
)

writeLines(do_lines, do_path)

stata_out <- system2(
  stata_bin,
  c("-q", "do", do_path),
  stdout = TRUE,
  stderr = TRUE
)

status <- attr(stata_out, "status")
if (!is.null(status) && status != 0L) {
  stop(
    paste(c("Stata command failed.", stata_out), collapse = "\n"),
    call. = FALSE
  )
}

result_lines <- grep("^CODEX_RESULT\\|", stata_out, value = TRUE)
if (length(result_lines) == 0L) {
  stop(
    paste(c("No `CODEX_RESULT` markers were found in Stata output.", stata_out), collapse = "\n"),
    call. = FALSE
  )
}

parsed <- do.call(
  rbind,
  lapply(
    result_lines,
    function(line) {
      parts <- strsplit(line, "\\|", fixed = FALSE)[[1L]]
      data.frame(
        section = parts[[2L]],
        key = parts[[3L]],
        value = as.numeric(trimws(parts[[4L]])),
        stringsAsFactors = FALSE
      )
    }
  )
)

core <- parsed[parsed$section == "core", ]
partial <- parsed[parsed$section == "partial", ]
lqte <- parsed[parsed$section == "lqte", ]

core_vec <- setNames(core$value, core$key)
partial_vec <- setNames(partial$value, partial$key)
lqte_df <- data.frame(
  quantile = as.numeric(sub("^q_", "", lqte$key)) / 100,
  estimate = lqte$value
)

cat("# Paste into tests/testthat/test-stata-parity.R\n\n")
cat("stata_core_golden <- ")
dput(core_vec)
cat("\n")

cat("stata_partial_golden <- ")
dput(partial_vec)
cat("\n")

cat("stata_lqte_reference <- ")
dput(lqte_df)
cat("\n")

cat(
  "# Note: Stata's `fuzzydid` command does not expose bootstrap failure counts\n",
  "# directly, so `n_misreps`/`share_failures` remain native regression checks.\n",
  sep = ""
)
