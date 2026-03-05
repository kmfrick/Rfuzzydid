# Rfuzzydid
R package for corrected estimators in fuzzy difference-in-differences designs proposed by de Chaisemartin and D'Haultfoeuille ([RES 2018](https://academic.oup.com/restud/article-abstract/85/2/999/4096388)).

## Current status

- Formula-first API: `fuzzydid(data, y ~ d + x, group = ..., time = ...)`.
- Native R backend supports `did`, `tc`, `cic`, `lqte`, and `partial` (with design restrictions matching Stata-style checks).
- Multi-period designs are supported through `group` + `group_forward`.
- Covariate-adjusted DID/TC supports `modelx`, `sieves`, and `sieveorder`.
- Bootstrap inference follows Stata parity defaults: fixed bootstrap seed `1`, failed/degenerate rep tracking (`n_reps`, `n_misreps`, `share_failures`), and sentinel-drop handling before SE/CI computation.
- A vendored Stata snapshot is shipped under `inst/stata/fuzzydid` for parity reference.

See `?fuzzydid` for option-level constraints and behavior.
