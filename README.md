## mimeo - Faithful simulations and replications of molecular and clinical data

[![Travis-CI Build Status](https://travis-ci.org/hruffieux/mimeo.svg?branch=master)](https://travis-ci.org/hruffieux/mimeo)
 
## Overview

**mimeo** is an R package providing functions to emulate molecular quantitative
trait locus data and clinical data associated with genetic variants under
user-specified association patterns. The data can be the result of pure 
simulation or can replicate real datasets supplied by the user (which may
replace real data when these cannot be shared for diverse privacy reasons). The 
data generation schemes are based on generally accepted principles of population 
genetics (Hardy–Weinberg equilibrium, linkage-disequilibrium, natural selection,
pleiotropic control, sparsity assumptions).

Details in H. Ruffieux, A. C. Davison, J. Hager, I. Irincheeva, Efficient 
inference for genetic association studies with multiple outcomes, *Biostatistics*, 
2017. 

## Installation

To install, run the following commands in R:

``` r
install.packages("devtools")
devtools::install_github("hruffieux/mimeo")
```

## License and authors

This software uses the GPL v2 license, see [LICENSE](LICENSE).
Authors and copyright are provided in [DESCRIPTION](DESCRIPTION).

## Issues

To report an issue, please use the [locus issue tracker](https://github.com/hruffieux/mimeo/issues) at github.com.