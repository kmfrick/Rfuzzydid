# Contributing

## Development Workflow

- Install package dependencies with `rv sync`, then run `R CMD check .` from the
  `rv` environment before opening a PR.
- Regenerate documentation with a temporary script that calls
  `roxygen2::roxygenise()` after updating roxygen comments.
- Run a temporary script that calls `testthat::test_dir("tests/testthat")` for a
  fast local test pass.
- For SRR readiness, run `srr::srr_stats_pre_submit(".")` and `pkgcheck`.

## Stata Parity Goldens

- Regular tests do not require Stata.
- Frozen Stata parity constants live directly in
  `tests/testthat/test-stata-parity.R`.
- If those constants need to be recovered, run Stata locally with the same
  fixture construction as the test and paste the resulting values into the
  hardcoded constants. Do not add committed Stata execution code to the test
  suite.
- The Stata commands used to recover the currently asserted values are:
  `fuzzydid y g t d, did tc cic nose` for `stata_core_golden`, then
  `matrix list e(b_LATE)`; and `fuzzydid y g t d, tc partial nose` for
  `stata_partial_golden`, then `matrix list e(b_LATE)`.

## Pull Requests

- Keep behavior changes covered by `testthat`.
- Update `NEWS.md` for user-visible changes.
- Prefer targeted patches over large formatting-only rewrites.
