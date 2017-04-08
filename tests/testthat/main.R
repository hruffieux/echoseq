rm(list = ls())

set.seed(123)

############################
## Simulate basic dataset ##
############################

n <- 100; p <- 300; d <- 100

## SNPs
##
sim_snps <- generate_snps(n, p)
repl_snps <- replicate_real_snps(n, sim_snps$snps, bl_lgth = p)

## Phenotypes
##
var_err <- runif(d, min = 0.1, max = 0.5)
sim_phenos <- generate_phenos(n, d, var_err, vec_rho = vec_rho)
repl_phenos <- replicate_real_phenos(n, sim_phenos$phenos)
