#!/usr/bin/env Rscript

status_msg <- function(fmt, ...) {
  cat(sprintf(paste0("[generate-stata-parity-goldens] ", fmt, "\n"), ...), file = stderr())
}

time_step <- function(label, expr) {
  start <- Sys.time()
  status_msg("%s...", label)
  value <- force(expr)
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  status_msg("%s finished in %.2fs", label, elapsed)
  invisible(value)
}

tail_lines <- function(path, n = 40L) {
  if (!file.exists(path)) {
    return(character(0))
  }

  lines <- readLines(path, warn = FALSE)
  utils::tail(lines, n)
}

stata_failure <- function(status, log_path) {
  msg <- c(sprintf("Stata command failed with exit status %s.", status))

  if (file.exists(log_path)) {
    msg <- c(msg, sprintf("Stata log: %s", log_path), "Log tail:", tail_lines(log_path))
  }

  stop(paste(msg, collapse = "\n"), call. = FALSE)
}

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
do_path <- file.path(tmp_dir, "parity.do")
log_path <- file.path(tmp_dir, "parity.log")

time_step("Writing parity fixture", {
  utils::write.csv(make_parity_fixture(), fixture_path, row.names = FALSE)
})

fixture_path_stata <- gsub("\\\\", "/", normalizePath(fixture_path, winslash = "/"))
log_path_stata <- gsub("\\\\", "/", normalizePath(log_path, winslash = "/", mustWork = FALSE))

do_lines <- c(
  "clear all",
  "set more off",
  sprintf("log using \"%s\", replace text name(codexlog)", log_path_stata),
  "capture which fuzzydid",
  "if _rc {",
  "  di as error \"CODEX_ERROR|missing_fuzzydid\"",
  "  capture log close codexlog",
  "  exit 111",
  "}",
  sprintf("import delimited using \"%s\", clear", fixture_path_stata),
  "display as text \"[stata-parity] running core estimates\"",
  "timer clear",
  "timer on 1",
  "quietly fuzzydid y g t d, did tc cic nose",
  "timer off 1",
  "display \"CODEX_RESULT|core|did|\" %21.15f el(e(b_LATE),1,1)",
  "display \"CODEX_RESULT|core|tc|\" %21.15f el(e(b_LATE),2,1)",
  "display \"CODEX_RESULT|core|cic|\" %21.15f el(e(b_LATE),3,1)",
  "display as text \"[stata-parity] running partial bounds\"",
  "timer on 2",
  "quietly fuzzydid y g t d, tc partial nose",
  "timer off 2",
  "display \"CODEX_RESULT|partial|TC_inf|\" %21.15f el(e(b_LATE),1,1)",
  "display \"CODEX_RESULT|partial|TC_sup|\" %21.15f el(e(b_LATE),2,1)",
  "display as text \"[stata-parity] running lqte estimates\"",
  "timer on 3",
  "quietly fuzzydid y g t d, lqte nose",
  "timer off 3",
  "forvalues i=1/19 {",
  "  local q = `i' * 5",
  "  display \"CODEX_RESULT|lqte|q_`q'|\" %21.15f el(e(b_LQTE),`i',1)",
  "}",
  "timer list",
  "log close codexlog",
  "exit, clear"
)

time_step("Writing Stata driver", {
  writeLines(do_lines, do_path)
})

status_msg("Temporary Stata log: %s", log_path)

stata_status <- time_step("Running Stata parity job", {
  system2(
    stata_bin,
    c("-b", "do", do_path),
    stdout = "",
    stderr = ""
  )
})

if (!identical(as.integer(stata_status), 0L)) {
  stata_failure(stata_status, log_path = log_path)
}

result_lines <- time_step("Reading Stata result markers", {
  if (!file.exists(log_path)) {
    stata_failure("missing-log", log_path = log_path)
  }

  grep("^CODEX_RESULT\\|", readLines(log_path, warn = FALSE), value = TRUE)
})

if (length(result_lines) == 0L) {
  stata_failure("no-markers", log_path = log_path)
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

status_msg("Parsed %d Stata markers", length(result_lines))

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
