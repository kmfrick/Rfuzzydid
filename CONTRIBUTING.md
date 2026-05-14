# Contributing

## Development Workflow

- Install package dependencies, then run `R CMD check .` before opening a PR.
- Regenerate documentation with `Rscript -e "roxygen2::roxygenise()"` after
  updating roxygen comments.
- Run `Rscript -e "testthat::test_dir('tests/testthat')"` for a fast local test
  pass.

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
