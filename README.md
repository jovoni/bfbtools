
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bfbtools

<!-- badges: start -->

[![R-CMD-check](https://github.com/jovoni/bfbtools/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jovoni/bfbtools/actions/workflows/R-CMD-check.yaml)
[![License:
GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

Linear-time algorithms deciding whether a chromosomal segment
copy-number count-vector is consistent with the breakage-fusion-bridge
(BFB) mechanism, and finding the closest BFB-admitting count-vector to a
noisy observed one. Ported from the reference implementation of Zakov,
Kinsella & Bafna (2013, *PNAS* 110(14):5546-5551), building on Kinsella
& Bafna (2012, *J. Comput. Biol.* 19(6):662-678).

`bfbtools` is deliberately sim-agnostic: it works on plain integer
count-vectors, with no dependency on any particular simulator or fitting
pipeline. It’s the detection engine underneath
[`bridges`](https://github.com/jovoni/bridges), which wraps it for use
directly on `cna_data` tibbles (per-cell copy-number profiles),
including chromosome-arm-aware handling.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("jovoni/bfbtools")
```

## Usage

``` r
library(bfbtools)

admits_bfb(c(6, 3, 5))               # TRUE
admits_bfb(c(14, 7, 18, 16, 9, 12))  # FALSE, a real near-miss example

nearest_bfb(c(14, 7, 18, 16, 9, 12))
#> $weight
#> [1] 0.9114639
#> $vector
#> [1] 14  7 19 17  9 13

bfb_distance_weighted(c(14, 7, 18, 16, 9, 12))
#> [1] 0.09270324
```

`admits_bfb()` is an exact, linear-time decision, not an approximation
or statistical test, and never needs a search budget. `nearest_bfb()`
finds the closest BFB-admitting vector under a Poisson noise model, for
scoring real (noisy) copy-number data.

See `vignette("bfbtools-algorithms")` for the full walkthrough,
including a worked-through explanation of what “admits BFB” actually
means and a caveat about trivial (non-amplified) count-vectors.

## Part of a two-package family

- **`bfbtools`** (this package), the core detection algorithms, GPL-3.
- [**`bridges`**](https://github.com/jovoni/bridges), simulation,
  phylogenetic fitting, and batch/tibble-aware detection built on top of
  this package. MIT licensed; depends on `bfbtools` via a normal package
  dependency, which does not make `bridges` GPL.

## Why a separate package?

`bfbtools`’s core algorithm (`src/decision.cpp`) is a line-by-line port
of the reference Java implementation at
[shay-zakov/BFBFinder](https://github.com/shay-zakov/BFBFinder), which
is GPL-3 licensed.

## License

GPL-3.
