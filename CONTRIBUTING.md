# Contributing

## Development Workflow

- Install package dependencies, then run `R CMD check .` before opening a PR.
- Regenerate documentation with `Rscript -e "roxygen2::roxygenise()"` after
  updating roxygen comments.
- Run `Rscript -e "testthat::test_dir('tests/testthat')"` for a fast local test
  pass.

## Stata Parity Goldens

- Regular tests do not require Stata.
- To refresh the frozen Stata parity constants, run
  `Rscript tools/generate-stata-parity-goldens.R`.
- The helper assumes the `stata` CLI and the Stata `fuzzydid` command are
  already installed locally.

## Pull Requests

- Keep behavior changes covered by `testthat`.
- Update `NEWS.md` for user-visible changes.
- Prefer targeted patches over large formatting-only rewrites.
