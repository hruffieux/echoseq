<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- First time: run usethis::use_readme_rmd() to create a pre-commit hook that 
prevents from committing if the README.Rmd has changed, but has not been 
re-knitted to generate an updated README.md -->

## echoseq: Faithful replication and simulation of molecular and clinical data <img src="man/figures/echoseq_logo.png" align="right" height="150"/>

<!-- Run for the R CMD checks, run usethis::use_github_actions() to set up the pipeline, possibly modify the .yaml file and then: -->

[![](https://travis-ci.org/hruffieux/echoseq.svg?branch=master)](https://travis-ci.org/hruffieux/echoseq)
[![R build
status](https://github.com/hruffieux/echoseq/workflows/R-CMD-check/badge.svg)](https://github.com/hruffieux/echoseq/actions)
[![License: GPL
v3](https://img.shields.io/badge/license-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
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

## Installation

To install, run the following commands in R:

``` r
if(!require(remotes)) install.packages("remotes")
remotes::install_github("hruffieux/echoseq")
```

## License and authors

This software uses the GPL v3 license, see [LICENSE](LICENSE). Authors
and copyright are provided in [DESCRIPTION](DESCRIPTION).

## Issues

To report an issue, please use the [echoseq issue
tracker](https://github.com/hruffieux/echoseq/issues) at github.com.
