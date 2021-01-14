<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- Add README.Rmd to .Rbuildignore -->

## echoseq: Faithful replication and simulation of molecular and clinical data

<!-- Run usethis::use_github_actions() to set up the pipeline and then -->

[![R build
status](https://github.com/hruffieux/echoseq/workflows/R-CMD-check/badge.svg)](https://github.com/hruffieux/echoseq/actions)
[![](https://travis-ci.org/hruffieux/echoseq.svg?branch=devel)](https://travis-ci.org/hruffieux/echoseq)
[![](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![License: GPL
v3](https://img.shields.io/badge/license-GPL%20v3-blue.svg)](https://cran.r-project.org/web/licenses/GPL%20v3)
[![](https://img.shields.io/badge/devel%20version-0.3.1-blue.svg)](https://github.com/hruffieux/echoseq)
[![](https://img.shields.io/github/languages/code-size/hruffieux/echoseq.svg)](https://github.com/hruffieux/echoseq)

## Overview

**echoseq** is an R package providing functions to emulate molecular
quantitative trait locus data, and clinical data associated with genetic
variants under user-specified association patterns. The data can be
obtained by pure simulation or can replicate real datasets supplied by
the user. Datasets simulated using **echoseq** may replace real datasets
when these cannot be shared for diverse privacy reasons. The data
generation schemes are based on generally accepted principles of
population genetics (Hardy–Weinberg equilibrium, linkage-disequilibrium,
natural selection, pleiotropic control, sparsity assumptions, epigenetic
control, etc).

## Warning

**This is a development branch**, it is not guaranteed to be stable at
any given time and features are subject to change. Please use the
[stable version](https://github.com/hruffieux/echoseq), unless you want
to test and report issues.

## Installation

To install, run the following commands in R:

``` r
if(!require(devtools)) install.packages("devtools")
devtools::install_github("hruffieux/echoseq", ref = "devel")
```

## License and authors

This software uses the GPL v3 license, see [LICENSE](LICENSE). Authors
and copyright are provided in [DESCRIPTION](DESCRIPTION).

## Issues

To report an issue, please use the [echoseq issue
tracker](https://github.com/hruffieux/echoseq/issues) at github.com.
