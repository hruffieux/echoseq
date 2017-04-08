# This file is part of the `mimeo` R package:
#     https://github.com/hruffieux/mimeo
#

#' Generate SNPs emulating real SNP data at hand.
#'
#' This function simulates SNPs from real SNP data based on their sample minor
#' allele frequencies and correlation structure.
#'
#' @param n Number of observations.
#' @param real_snps Matrix of real SNPs (rows observations, columns SNP
#'   variables), without missing values. The entries must be 0, 1 or 2.
#' @param bl_lgth Number of variables per block for reproducing the dependence
#'   structure of real SNPs. Must be between 2 and p. Must be small enough
#'   (e.g. 1000) for tractability reasons.
#' @param p Number of SNPs. If \code{NULL}, the total number of SNPs
#'   available in the real_snps matrix is used (default).
#' @param maf_thres Lower bound for sample minor allele frequencies. Simulated
#'   SNPs with lower sample minor allele frequencies are excluded (as a result
#'   the number of simulated SNPs can be smaller than p). Default is \code{NULL}
#'   for no exclusion.
#' @param n_cpus Number of CPUs used when simulating SNP by blocks. Set to 1 for
#'   serial execution.
#' @param user_seed Seed set for reproducibility. Default is \code{NULL}, no
#'   seed set.
#'
#' @return An object of class "\code{sim_snps}".
#'  \item{snps}{Matrix containing the generated SNP data.}
#'  \item{vec_maf}{Vector containing the SNP sample minor allele frequencies.}
#'
#' @examples
#' user_seed <- 123; set.seed(user_seed)
#' n <- 500; p <- 7500
#' cor_type <- "autocorrelated"; vec_rho <- runif(100, min = 0.25, max = 0.95)
#' list_fake_real_snps <- generate_snps(n, p, cor_type, vec_rho, n_cpus = 2,
#'                                      user_seed = user_seed)
#' list_snps <- replicate_real_snps(n, list_fake_real_snps$snps, bl_lgth = 100,
#'                                  n_cpus = 2, user_seed = user_seed)
#'
#' @seealso \code{\link{generate_snps}}, \code{\link{generate_phenos}},
#'   \code{\link{replicate_real_phenos}}, \code{\link{generate_dependence}}
#'
#' @export
#'
replicate_real_snps <- function(n, real_snps, bl_lgth, p = NULL, maf_thres = NULL,
                                n_cpus = 1, user_seed = NULL) {


  check_structure_(user_seed, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(user_seed)){
    RNGkind("L'Ecuyer-CMRG") # ensure reproducibility when using mclapply
    set.seed(user_seed)
  }

  check_structure_(real_snps, "matrix", "numeric")

  check_structure_(maf_thres, "vector", "numeric", 1, null_ok = TRUE)
  if(!is.null(maf_thres)) check_zero_one_(maf_thres)

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

  # n can be larger than the number of available observations for the real snps
  check_structure_(n, "vector", "numeric", 1)
  check_natural_(n)

  check_structure_(p, "vector", "numeric", 1, null_ok = TRUE)
  if(!is.null(p)) check_natural_(p)

  if (!is.null(p)) {
    p_av <- ncol(real_snps)
    if(p > p_av)
      stop(paste("Provided n greater than number of snps available: ",
                 p, " > ", p_av, sep = ""))
    real_snps <- real_snps[, 1:p]
  } else {
    p <- ncol(real_snps)
  }

  check_structure_(bl_lgth, "vector", "numeric", 1)
  check_natural_(bl_lgth)

  if (bl_lgth == 1) stop("Provided block length must be larger than 1")
  if (bl_lgth > p)
    stop(paste("Provided block length must be smaller than the number of SNPs available: ",
               bl_lgth, " > ", p, sep = ""))

  vec_real_maf <- apply(real_snps, 2, mean) / 2

  n_real <- nrow(real_snps)
  n_bl <- floor(p / bl_lgth)
  ind_bl <- make_chunks_(1:p, n_bl)

  snps <- parallel::mclapply(1:n_bl, function(bl) {

    p_bl <- length(ind_bl[[bl]])

    # we add some noise to avoid undefined correlation in case of constant phenotypes.
    R <- cor(real_snps[, ind_bl[[bl]]] + matrix(rnorm(n_real*p_bl), nrow = n_real))
    R <- Matrix::nearPD(R, corr = TRUE, do2eigen = TRUE)$mat

    L <- t(chol(R))
    tZ <- matrix(rnorm(n * p_bl), nrow = p_bl, ncol = n)
    X <- t(L %*% tZ) # Gaussian variables

    snps <- matrix(1, nrow = n, ncol = p_bl)
    for(j in 1:p_bl) {
      maf <- vec_real_maf[ind_bl[[bl]]][j]
      snps[X[,j] < qnorm((1 - maf)^2), j] <- 0
      snps[X[,j] > qnorm(1 - maf^2), j] <- 2
    }
    snps
  }, mc.cores = n_cpus)

  snps <- do.call(cbind, snps)
  rownames(snps) <- paste("ind_", 1:n, sep = "")
  colnames(snps) <- paste("snp_", 1:p, sep = "")

  vec_maf <- apply(snps, 2, mean) / 2
  names(vec_maf) <- colnames(snps)

  if (!is.null(maf_thres)) {
    ind_rare <- vec_maf < maf_thres
    snps <- snps[, !ind_rare]
    vec_maf <- vec_maf[!ind_rare]
  }

  list_snps <- create_named_list_(snps, vec_maf)
  class(list_snps) <- "sim_snps"
  list_snps
}



#' Generate phenotypes emulating real phenotypic data at hand.
#'
#' This function simulates phenotypes from real phenotypic data based on their
#' sample correlation structure. If binary data are provided, will simulate
#' latent Gaussian data from them (binary data can then be obtained using the
#' \code{\link{generate_dependence}} function).
#'
#' @param n Number of observations.
#' @param real_phenos Matrix of real phenotypes (rows observations, columns
#'   phenotypic variables), without missing values.
#' @param input_family Phenotype distribution assumption for the phenotypes
#'   provided in real_phenos. Must be either "\code{gaussian}" or
#'   "\code{binomial}" for binary phenotypes.
#' @param bl_lgth Number of variables per block for reproducing the dependence
#'   structure of real phenotypes. Must be between 2 and d. Must be small enough
#'   (e.g. 1000) for tractability reasons. Default is \code{NULL} for a single
#'   block.
#' @param d Number of phenotypes. If \code{NULL}, the total number of phenotypes
#'   available in the real_phenos matrix is used (default).
#' @param n_cpus Number of CPUs used when simulating phenotypes by blocks. Set
#'   to 1 for serial execution.
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
#' n <- 500; d <- 1000
#' cor_type <- "equicorrelated"; vec_rho <- runif(100, min = 0.25, max = 0.95)
#'
#' # Provided phenotypes assumed to be normally distributed
#' var_err <- runif(d, min = 0.1, max = 0.4)
#' list_fake_real_phenos <- generate_phenos(n, d, var_err, cor_type = cor_type,
#'                                          vec_rho = vec_rho, n_cpus = 2,
#'                                          user_seed = user_seed)
#' list_phenos <- replicate_real_phenos(n, list_fake_real_phenos$phenos,
#'                                      input_family = "gaussian", bl_lgth = 100,
#'                                      n_cpus = 2, user_seed = user_seed)
#'
#' @seealso \code{\link{generate_snps}}, \code{\link{replicate_real_snps}},
#'   \code{\link{generate_phenos}}, \code{\link{generate_dependence}}
#'
#' @export
#'
replicate_real_phenos <- function(n, real_phenos, input_family = "gaussian",
                                  bl_lgth = NULL, d = NULL, n_cpus = 1,
                                  user_seed = NULL) {

  check_structure_(user_seed, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(user_seed)){
    RNGkind("L'Ecuyer-CMRG") # ensure reproducibility when using mclapply
    set.seed(user_seed)
  }

  stopifnot(input_family %in% c("gaussian", "binomial"))

  check_structure_(real_phenos, "matrix", "numeric")

  var_err <- apply(real_phenos, 2, var)

  if (input_family == "binomial") {
    if(!identical(as.vector(real_phenos), as.numeric(as.logical(real_phenos))))
      stop("real_phenos must be a binary matrix when family is set to binomial.")
  }

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

  check_structure_(n, "vector", "numeric", 1) # can be larger than the number of available observations in real_pheno!
  check_natural_(n)

  check_structure_(d, "vector", "numeric", 1, null_ok = TRUE)
  if(!is.null(d)) check_natural_(d)

  if (!is.null(d)) {
    d_av <- ncol(real_phenos)
    if(d > d_av)
      stop(paste("Provided n greater than number of phenotypes available: ",
                 d, " > ", d_av, sep = ""))
    real_phenos <- real_phenos[, 1:d]
  } else {
    d <- ncol(real_phenos)
  }

  if (is.null(bl_lgth)) {  # bl_lgth = NULL, means one block

    bl_lgth <- d

  } else {

    check_structure_(bl_lgth, "vector", "numeric", 1)
    check_natural_(bl_lgth)

    if (bl_lgth == 1) stop("Provided block length must be larger than 1")
    if (bl_lgth > d)
      stop(paste("Provided block length must be smaller or equal to the number ",
                 "of phenotypes available: ", bl_lgth, " > ", d, sep = ""))
  }

  n_real <- nrow(real_phenos)
  n_bl <- floor(d / bl_lgth)
  ind_bl <- make_chunks_(1:d, n_bl)

  phenos <- parallel::mclapply(1:n_bl, function(bl) {
    d_bl <- length(ind_bl[[bl]])

    # we add some noise to avoid undefined correlation in case of constant phenotypes.
    R <- cor(real_phenos[, ind_bl[[bl]]] + matrix(rnorm(n_real*d_bl), nrow = n_real))
    R <- Matrix::nearPD(R, corr = TRUE, do2eigen = TRUE)$mat
    L <- t(chol(R))
    tZ <- matrix(sapply(var_err[ind_bl[[bl]]], function(ve) rnorm(n, 0, sqrt(ve))),
                 ncol = n, byrow = TRUE)
    as.matrix(t(L %*% tZ))

  }, mc.cores = n_cpus)

  phenos <- do.call(cbind, phenos)

  rownames(phenos) <- paste("ind_", 1:n, sep = "")
  colnames(phenos) <- paste("pheno_", 1:d, sep = "")

  ind_bl <- NULL # block chuncks do not necessarily represent blocks of correlated phenotypes.

  var_err <- apply(phenos, 2, var) # empirical error variance
  names(var_err) <- colnames(phenos)

  list_phenos <- create_named_list_(phenos, var_err, ind_bl)
  class(list_phenos) <- "sim_phenos"
  list_phenos
}
