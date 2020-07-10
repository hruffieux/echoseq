# This file is part of the `echoseq` R package:
#     https://github.com/hruffieux/echoseq
#

#' Generate SNP data with prespecified spatial correlation structure.
#'
#' This function generates SNPs under Hardy-Weinberg equilibrium with specific
#' block correlation structure and minor allele frequencies.
#'
#' @param n Number of observations.
#' @param p Number of SNPs.
#' @param cor_type String describing the type of dependence structure. The SNPs
#'   can \code{autocorrelated}, \code{equicorrelated}. Set to \code{NULL} for
#'   independent SNPs.
#' @param vec_rho Vector of correlation coefficients. Its length determines the
#'   number of blocks of correlated SNPs. Must be smaller than p. Set to
#'   \code{NULL} if independent SNPs.
#' @param vec_maf Vector of length p containing the reference minor allele
#'   frequencies used to generate the SNPs. If \code{NULL}, the minor allele
#'   frequencies drawn uniformly at random between 0.05 and 0.5.
#' @param n_cpus Number of CPUs used when simulating correlated SNP blocks.
#'   Ignored if independent SNPs. Set to 1 for serial execution.
#' @param user_seed Seed set for reproducibility. Default is \code{NULL}, no
#'   seed set.
#'
#' @return An object of class "\code{sim_snps}".
#'  \item{snps}{Matrix containing the generated SNP data.}
#'  \item{vec_maf}{Vector containing the SNP sample minor allele frequencies.}
#'
#' @examples
#' user_seed <- 123; set.seed(user_seed)
#' n <- 500; p <- 10000
#' cor_type <- "autocorrelated"; vec_rho <- runif(100, min = 0.25, max = 0.95)
#' list_snps <- generate_snps(n, p, cor_type, vec_rho, n_cpus = 2,
#'                            user_seed = user_seed)
#'
#' @seealso \code{\link{replicate_real_snps}}, \code{\link{generate_phenos}},
#'   \code{\link{replicate_real_phenos}}, \code{\link{generate_dependence}}
#'
#' @export
#'
generate_snps <- function(n, p, cor_type = NULL, vec_rho = NULL, vec_maf = NULL,
                          n_cpus = 1, user_seed = NULL) {

  check_structure_(user_seed, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(user_seed)){
    RNGkind("L'Ecuyer-CMRG") # ensure reproducibility when using mclapply
    set.seed(user_seed)
  }

  check_structure_(n, "vector", "numeric", 1)
  check_natural_(n)

  check_structure_(p, "vector", "numeric", 1)
  check_natural_(p)

  if (!is.null(cor_type))
    stopifnot(cor_type %in% c("autocorrelated", "equicorrelated"))

  if(is.null(vec_maf)) {
    vec_maf <- runif(p, 0.05, 0.5)
  } else {
    check_structure_(vec_maf, "vector", "numeric", p)
    check_zero_one_(vec_maf)
  }

  if (is.null(cor_type)) {

    if (n_cpus > 1)
      warning("n_cpus is ignored when the SNPs are generated independently of one another.")

    snps <- sapply(vec_maf, function(maf) rbinom(n, 2, maf)) # Hardy-Weinberg equilibrium

  } else {

    check_structure_(vec_rho, "vector", "numeric")
    if(cor_type == "equicorrelated") check_zero_one_(vec_rho)
    else check_zero_one_(abs(vec_rho))

    if(length(vec_rho) > p)
      stop(paste("Provided number of blocks of correlated SNPs, length(vec_rho), ",
                 "must be smaller than the number of SNPs, p: ",
                 length(vec_rho), " > ", p, sep = ""))

    check_structure_(n_cpus, "vector", "numeric", 1)
    check_natural_(n_cpus)

    if (n_cpus > 1) {
      n_cpus_avail <- parallel::detectCores()
      if (n_cpus > n_cpus_avail) {
        n_cpus <- n_cpus_avail
        warning(paste("The number of CPUs specified exceeds the number of CPUs ",
                      "available on the machine. The latter has been used instead.", sep=""))
      }
    }

    n_bl <- length(vec_rho)
    ind_bl <- make_chunks_(1:p, n_bl)

    snps <- parallel::mclapply(1:n_bl, function(bl) {

      p_bl <- length(ind_bl[[bl]])
      rho_bl <- vec_rho[[bl]]
      R <- matrix(NA, nrow = p_bl, ncol = p_bl)

      if (cor_type == "autocorrelated") {
        for( i in 1:p_bl ){
          for( j in i:p_bl ){
            R[i,j] <- rho_bl^abs(i - j)
            R[j,i] <- R[i,j]
          }
        }
      } else { # equicorrelated
        R[] <- rho_bl
      }

      diag(R) <- 1

      L <- t(chol(R))
      tZ <- matrix(rnorm(n * p_bl), ncol = n)
      X <- t(L %*% tZ)

      snps <- matrix(1, nrow = n, ncol = p_bl)

      for(j in 1:p_bl) {
        maf <- vec_maf[ind_bl[[bl]]][j]
        snps[X[,j] < qnorm((1 - maf)^2), j] <- 0
        snps[X[,j] > qnorm(1 - maf^2), j] <- 2
      }
      snps
    }, mc.cores = n_cpus)

    snps <- cbind_fill_matrix(snps)
  }

  snps <- as.matrix(snps)
  if (n == 1) snps <- t(snps)

  rownames(snps) <- paste("ind_", 1:n, sep = "")
  colnames(snps) <- paste("snp_", 1:p, sep = "")

  vec_maf <- apply(snps, 2, mean) / 2 # empirical maf
  names(vec_maf) <- colnames(snps)

  list_snps <- create_named_list_(snps, vec_maf)
  class(list_snps) <- "sim_snps"
  list_snps
}



#' Generate phenotypic data with prespecified block-wise correlation structure.
#'
#' This function generates Gaussian phenotypes with specific block correlation
#' structure. If binary phenotypes are wanted, must be used in conjunction with
#' the \code{\link{generate_dependence}} function: the present function
#' simulates latent Gaussian phenotypes which will be used in
#' \code{\link{generate_dependence}} to generate binary phenotypes associated
#' with SNPs from a probit model.
#'
#' @param n Number of observations.
#' @param d Number of phenos.
#' @param var_err Vector of length 1 or d containing the variances of the
#'   (latent) Gaussian distributions used to generate the phenotypes. If of
#'   length 1, the value is repeated d times. If binary data are targeted
#'   (using the probit modelling within the \code{\link{generate_dependence}}
#'   function), var_err is usually set to 1 (no impact on inference).
#' @param cor_type String describing the type of dependence structure. The
#'   phenotypes can \code{autocorrelated}, \code{equicorrelated}. Set to
#'   \code{NULL} for independent phenotypes.
#' @param vec_rho Vector of correlation coefficients. Its length determines the
#'   number of blocks of correlated phenotypes. Must be smaller than d. Set to
#'   \code{NULL} if independent phenotypes.
#' @param n_cpus Number of CPUs used when simulating correlated phenotype blocks.
#'   Ignored if independent phenotypes. Set to 1 for serial execution.
#' @param user_seed Seed set for reproducibility. Default is \code{NULL}, no
#'   seed set.
#'
#' @return An object of class "\code{sim_phenos}".
#'  \item{phenos}{Matrix containing the generated phenotypic data.}
#'  \item{var_err}{Vector containing the sample phenotypic variances.}
#'  \item{ind_bl}{List of length given by the number of blocks, containing the
#'                indices of the phenotypes in each block. Is \code{NULL} if
#'                \code{cor_type} is \code{NULL}.}
#'
#' @examples
#' user_seed <- 123; set.seed(user_seed)
#' n <- 500; d <- 10000
#' cor_type <- "equicorrelated"; vec_rho <- runif(100, min = 0.25, max = 0.95)
#'
#' list_phenos <- generate_phenos(n, d, cor_type = cor_type, vec_rho = vec_rho,
#'                                n_cpus = 2, user_seed = user_seed)
#'
#' @seealso \code{\link{generate_snps}}, \code{\link{replicate_real_snps}},
#'   \code{\link{replicate_real_phenos}}, \code{\link{generate_dependence}}
#'
#' @export
#'
generate_phenos <- function(n, d, var_err = 1, cor_type = NULL, vec_rho = NULL,
                            n_cpus = 1, user_seed = NULL) {

  check_structure_(user_seed, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(user_seed)){
    RNGkind("L'Ecuyer-CMRG") # ensure reproducibility when using mclapply
    set.seed(user_seed)
  }

  check_structure_(n, "vector", "numeric", 1)
  check_natural_(n)

  check_structure_(d, "vector", "numeric", 1)
  check_natural_(d)

  check_structure_(var_err, "vector", "numeric", c(1, d))
  check_positive_(var_err)
  if (length(var_err) == 1) var_err <- rep(var_err, d)

  if (!is.null(cor_type))
    stopifnot(cor_type %in% c("autocorrelated", "equicorrelated"))

  check_structure_(n_cpus, "vector", "numeric", 1)
  check_natural_(n_cpus)

  if (is.null(cor_type)) {

    if (n_cpus > 1)
      warning("n_cpus is ignored when the phenotypes are generated independently of one another.")

    phenos <- sapply(var_err, function(var_err) rnorm(n, 0, sqrt(var_err)))

    ind_bl <- NULL

  } else {

    check_structure_(vec_rho, "vector", "numeric")
    if(cor_type == "equicorrelated") check_zero_one_(vec_rho)
    else check_zero_one_(abs(vec_rho))

    if(length(vec_rho) > d)
      stop(paste("Provided number of blocks of correlated phenotypes, length(vec_rho), ",
                 "must be smaller than the number of phenotypes, d: ",
                 length(vec_rho), " > ", d, sep = ""))


    if (n_cpus > 1) {
      n_cpus_avail <- parallel::detectCores()
      if (n_cpus > n_cpus_avail) {
        n_cpus <- n_cpus_avail
        warning(paste("The number of CPUs specified exceeds the number of CPUs ",
                      "available on the machine. The latter has been used instead.", sep=""))
      }
    }

    n_bl <- length(vec_rho)
    ind_bl <- make_chunks_(1:d, n_bl)

    phenos <- parallel::mclapply(1:n_bl, function(bl) {

      d_bl <- length(ind_bl[[bl]])
      rho_bl <- vec_rho[[bl]]
      R <- matrix(NA, nrow = d_bl, ncol = d_bl)

      if (cor_type == "autocorrelated") {
        for( i in 1:d_bl ){
          for( j in i:d_bl ){
            R[i,j] <- rho_bl^abs(i - j)
            R[j,i] <- R[i,j]
          }
        }
      } else { # equicorrelated
        R[] <- rho_bl
      }
      diag(R) <- 1
      L <- t(chol(R))
      tZ <- matrix(sapply(var_err[ind_bl[[bl]]], function(ve) rnorm(n, 0, sqrt(ve))),
                   ncol = n, byrow = TRUE)
      as.matrix(t(L %*% tZ))

    }, mc.cores = n_cpus)

    phenos <- cbind_fill_matrix(phenos)
  }

  phenos <- as.matrix(phenos)
  if (n == 1) phenos <- t(phenos)

  rownames(phenos) <- paste("ind_", 1:n, sep = "")
  colnames(phenos) <- paste("pheno_", 1:d, sep = "")

  var_err <- apply(phenos, 2, var) # empirical error variance
  names(var_err) <- colnames(phenos)

  list_phenos <- create_named_list_(phenos, var_err, ind_bl)

  class(list_phenos) <- "sim_phenos"
  list_phenos
}
