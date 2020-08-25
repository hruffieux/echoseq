## echoseq: Faithful replication and simulation of molecular and clinical data

[![Travis-CI Build Status](https://travis-ci.org/hruffieux/echoseq.svg?branch=master)](https://travis-ci.org/hruffieux/echoseq)

## Overview

**echoseq** is an R package providing functions to emulate molecular quantitative
trait locus data, and clinical data associated with genetic variants under
user-specified association patterns. The data can be obtained by pure simulation 
or can replicate real datasets supplied by the user. Datasets simulated using 
**echoseq** may replace real datasets when these cannot be shared for diverse 
privacy reasons. The data generation schemes are based on generally accepted 
principles of population genetics (Hardy–Weinberg equilibrium, 
linkage-disequilibrium, natural selection, pleiotropic control, sparsity 
assumptions, epigenetic control, etc).


## Installation

To install, run the following commands in R:

``` r
if(!require(devtools)) install.packages("devtools")
devtools::install_github("hruffieux/echoseq")
```

## License and authors

This software uses the GPL v3 license, see [LICENSE](LICENSE).
Authors and copyright are provided in [DESCRIPTION](DESCRIPTION).

## Issues

To report an issue, please use the [echoseq issue tracker](https://github.com/hruffieux/echoseq/issues) at github.com.