# Rfuzzydid 1.1.1

- Fixed CRAN no-Suggests checks by making vignette rendering work when
  optional vignette-only packages are unavailable.
- Added rOpenSci-oriented package checks, coverage, pkgdown, and CI
  improvements.
- Added a package citation and clarified license metadata.
- Documented reproducible bootstrap seeds and added regression coverage.
- Simplified Stata parity maintenance by keeping frozen constants in the test
  suite and dropping committed generator tooling.

# Rfuzzydid 1.1.0

- Removed the vendored Stata package sources from the installed package.
- Kept Stata parity checks as frozen test goldens instead of a runtime test
  dependency.
- Added missing package metadata and user-facing contribution guidance.
- Expanded documentation for `fuzzydid()` and its S3 methods.
