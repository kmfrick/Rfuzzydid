# Rfuzzydid
R package for corrected estimators in fuzzy difference-in-differences designs proposed by de Chaisemartin and D'Haultfoeuille ([RES 2018](https://academic.oup.com/restud/article-abstract/85/2/999/4096388)).

## Current status

- Formula-first API: `fuzzydid(data, y ~ d, group = ..., time = ...)`.
- Native R backend for two-period `did`/`tc`/`cic` estimation.
- Stata delegation has been removed from runtime paths.
- Parity tests compare native estimates against frozen Stata golden fixtures.

See `?fuzzydid` for options and backend behavior.
