# This file is part of the `echoseq` R package:
#     https://github.com/hruffieux/echoseq
#

# Internal function setting the pattern of associations used by the
# `generate_dependence` function.
#
set_pattern_ <- function(d, p, ind_d0, ind_p0, vec_prob_sh, block_phenos, chunks_ph) {


  if (block_phenos) {

    if (is.null(chunks_ph)) { # no imposed correlation block structure (either indep
      # or correlation from real phenotypes). creates artificial chunks.
      n_chunks_ph <- length(vec_prob_sh)
      chunks_ph <- make_chunks_(1:d, n_chunks_ph)
    } else {
      n_chunks_ph <- length(chunks_ph)
    }

    pat <- matrix(FALSE, nrow = p, ncol = d)

    for(ind_j in ind_p0) {

      # random permutation, so that two different SNPs can be associated with a given
      # block with different probabilities
      vec_prob_sh_perm <- sample(vec_prob_sh, n_chunks_ph, replace = TRUE)

      for(ch in 1:n_chunks_ph) {
        ind_ch_d0 <- intersect(ind_d0, chunks_ph[[ch]])
        pat[ind_j, ind_ch_d0] <- sample(c(TRUE, FALSE),
                                        length(ind_ch_d0), replace = TRUE,
                                        prob=c(vec_prob_sh_perm[ch], 1 - vec_prob_sh_perm[ch]))
      }

      if (all(!pat[ind_j, ind_d0]))
        pat[ind_j, ind_d0][sample(1:length(ind_d0), 1)] <- TRUE # each active covariate must be
      # associated with at least one response
    }

  } else {

    d0 <- length(ind_d0)
    p0 <- length(ind_p0)

    if (length(vec_prob_sh) == 1)
      vec_prob_sh <- rep(vec_prob_sh, p0)

    pat <- matrix(FALSE, nrow = p, ncol = d)

    for(j in 1:p0) {

      ind_j <- ind_p0[j]

      pat[ind_j, ind_d0] <- sample(c(TRUE, FALSE), d0, replace = TRUE,
                                      prob=c(vec_prob_sh[j], 1 - vec_prob_sh[j]))

      if (all(!pat[ind_j, ind_d0]))
        pat[ind_j, ind_d0][sample(1:length(ind_d0), 1)] <- TRUE # each active covariate must be
                                                                # associated with at least one response
    }

  }

  for(ind_k in ind_d0) {
    if (all(!pat[ind_p0, ind_k]))
      # each active covariate must be associated with at least one response
      pat[ind_p0, ind_k][sample(1:length(ind_p0), 1)] <- TRUE
  }
  pat
}


# Internal function setting the effect sizes used by the `generate_dependence`
# function.
#
generate_eff_sizes_ <- function(d, phenos_act, snps_act, ind_d0, ind_p0, pat,
                                vec_prob_sh, vec_maf, pve_per_snp, max_tot_pve,
                                var_err, block_phenos, chunks_ph) {

  # pve_per_snp average variance explained per snp
  p <- length(vec_maf)

  var_snps_act <- apply(snps_act, 2, var)
  bool_cst <- var_snps_act == 0
  if (any(bool_cst)) { ### PROB VECTOR NEEDED? OR USE CUTOFF ZERO? ### centered logistic?

    if (sum(bool_cst) < 50) {
      message <- paste0("SNP(s) with id ", paste0(ind_p0, collapse = " "),
                       " constant. Effect(s) on the phenotypes removed.\n")

      if (all(bool_cst))
        stop(paste0(message, "No remaining ``active'' SNP, please change ind_p0 ",
                   "(now empty).\n"))
    } else {
      message <- paste0(sum(bool_cst), " SNPs constant. Effects on the phenotypes ",
                       "removed.\n")
      if (all(bool_cst))
        stop(paste0(message, "No remaining ``active'' SNP, please change ",
                   "ind_p0 (now empty).\n"))
    }

    warning(message)

    if (is.null(pat))
      ind_p0 <- ind_p0[!bool_cst]
  }

  var_phenos_act <- apply(phenos_act, 2, var)
  bool_cst <- var_phenos_act == 0
  if (any(bool_cst)) {

    if (sum(bool_cst) < 50) {
      message <- paste0("Phenotype(s) with id ", paste0(ind_d0, collapse = " "),
                       " constant. Association(s) with the SNPs removed.\n")

      if (all(bool_cst))
        stop(paste0(message, "No remaining ``active'' phenotype, please change ",
                   "ind_d0 (now empty).\n"))
    } else {
      message <- paste0(sum(bool_cst), " phenotype(s) constant. Association(s) ",
                       "with the SNPs removed.\n")

      if (all(bool_cst))
        stop(paste0(message, "No remaining ``active'' phenotype, please change ",
                   "ind_d0 (now empty).\n"))
    }

    warning(message)

    if (is.null(pat))
      ind_d0 <- ind_d0[!bool_cst]
  }

  if (is.null(pat))
    pat <- set_pattern_(d, p, ind_d0, ind_p0, vec_prob_sh, block_phenos, chunks_ph)

  check_structure_(vec_maf, "vector", "numeric", p)
  check_zero_one_(vec_maf)

  check_structure_(var_err, "vector", "numeric", d)
  check_positive_(var_err[ind_d0])

  max_per_resp <- max(colSums(pat))
  eps <- .Machine$double.eps^0.75
  if (is.null(pve_per_snp)) {
    # sets pve_per_snp to the max possible so that the tot_pve for all responses are below 1.
    if (is.null(max_tot_pve)) max_tot_pve <- 1 - eps

    pve_per_snp <- max_tot_pve / max_per_resp
  } else {
    check_structure_(pve_per_snp, "vector", "numeric", 1)
    check_zero_one_(pve_per_snp)

    if (max_per_resp * pve_per_snp > 1 - eps)
      stop(paste0("Provided average proportion of variance explained per SNP too ",
                 "high, would lead to a total genetic variance explained above ",
                 "100% for at least one response. \n Setting pve_per_snp < 1 / length(ind_p0) ",
                 "will work for any pattern. \n"))
  }

  beta <- matrix(0.0, nrow = p, ncol = d)

  # The proportion of phenotypic variance explained per SNP for phenotype k is
  # computed using the assumption of independent SNPs for simplicity (such an
  # assumption is commonly made for such computation). In this case we have, for
  # SNP X_j (assuming that the true effects \beta_{jk} are known):
  #
  # pve_j = var(X_j \beta_{jk} | \beta_{jk}) / (tot_var_expl_k + var_err_k),
  #
  # where tot_var_expl_k is the total variance explained by the SNPs for
  # phenotype k.
  #
  # We estimate the numerator using
  # var_hat(X_j \beta_{jk} | \beta_{jk}) = \beta_{jk}^2 * var_hat(X_j) =
  #                                        \beta_{jk}^2 * 2 * m_hat_j * (1 - m_hat_j),
  # where m_hat_j is the empricial minor allele frequency of X_j.
  #
  # We therefore set the regression effects as
  #
  # \beta_{jk} <- sqrt{ pve_j * (tot_var_expl_k + var_err_k) / [ 2 * m_hat_j * (1 - m_hat_j) ] }
  #
  # Natural selection: the relationship between \beta_{jk} and m_hat_j is
  # roughly inverse (for m_hat_j between 0 and 0.5).

  beta[, ind_d0] <- sapply(ind_d0, function(k) {

    p0_k <- sum(pat[,k])
    vec_pve_per_snp <- rbeta(p0_k, shape1 = 2, shape2 = 5) # positively skewed Beta distribution,
                                                           # to give more weight to smaller effect sizes
    vec_pve_per_snp <- vec_pve_per_snp / sum(vec_pve_per_snp) * pve_per_snp * p0_k

    tot_var_expl_k <- pve_per_snp * p0_k * var_err[k] / (1 - pve_per_snp * p0_k)

    vec_maf_act <- vec_maf[pat[,k]]
    vec_var_act <- 2 * vec_maf_act * (1 - vec_maf_act)

    beta_k <- rep(0.0, p)
    beta_k[pat[,k]] <- sqrt((tot_var_expl_k + var_err[k]) * vec_pve_per_snp / vec_var_act)

    # switches signs with probabilty 0.5
    beta_k[pat[,k]] <- sample(c(1, -1), p0_k, replace = TRUE) * beta_k[pat[,k]]

    beta_k
  })

  create_named_list_(beta, pat, pve_per_snp)

}

#' Generate pleiotropic associations between SNPs and phenotypes.
#'
#' This function sets the association pattern and the effect sizes between SNP
#' and phenotype objects previously obtained from the functions
#' \code{\link{generate_snps}} or \code{\link{replicate_real_snps}}, and
#' \code{\link{generate_phenos}} or \code{\link{replicate_real_phenos}}. It
#' therefore adds a genetic contribution to the phenotypic data.
#'
#' The user can provide using the argument \code{vec_prob_sh} a selection of
#' probabilities describing the propensity with which a given active SNP (i.e.,
#' associated with at least one phenotype) will be associated with active
#' phenotypes (i.e., associated with at least one SNP). If \code{block_phenos}
#' is \code{FALSE} (default), the association pattern is created independently
#' of any structure in the phenotype matrix. If \code{block_phenos} is
#' \code{TRUE}, then if the phenotypes have been generated with some block
#' correlation structure, this block structure will be used to specify the
#' correlation pattern, else, if the phenotypes were generated independently one
#' from another, blocks are defined articially and the number of blocks
#' corresponds to the length of the vector \code{vec_prob_sh}). More precisely,
#' for each active SNP and each phenotypic block, a value from this vector is
#' selected uniformly at random; for instance a large probability implies that
#' the SNPs is highly likely to be associated with each active phenotype in the
#' block. If a single value is provided, all active SNPs will have the same
#' probability to be associated with active phenotypes of all blocks.
#'
#' The user can provide either argument \code{pve_per_snp}, specifying the
#' average proportion of phenotypic variance explained per active SNP for a
#' given active phenotype, or \code{max_tot_pve}, specifying the maximum value
#' for an active phenotype of its proportion of variance explained by the
#' cummulated genetic effects. If both \code{pve_per_snp} and \code{max_tot_pve}
#' are \code{NULL}, the proportion of phenotypic variance explained per SNP is
#' set to its maximum value so that the total proportion of variance explained
#' for the phenotypes are all below 1. Individual proportions of variance
#' explained are drawn from a Beta distribution with shape parameters 2 and 5,
#' putting more weights on smaller effects.
#'
#' If family is "\code{binomial}", the phenotypes are generated from a probit
#' model, and the phenotypic variance explained by the SNPs is with respect to
#' the latent Gaussian variables involved in the probit model.
#'
#' @param list_snps An object of class "list_snps" or "sim_snps" containing
#'   SNPs and their corresponding sample minor allele frequencies. It must be
#'   obtained from the function \code{\link{convert_snps}},
#'   \code{\link{generate_snps}} or \code{\link{replicate_real_snps}}.
#' @param list_phenos An object of class "list_phenos" or "sim_phenos"
#'   containing phenotypic data variables, their sample variance and block
#'   structure information. It must be obtained from the function
#'   \code{\link{convert_phenos}}, \code{\link{generate_phenos}} or
#'   \code{\link{replicate_real_phenos}}.
#' @param ind_d0 A vector of indices specifying the position of the "active"
#'   phenotypes (i.e., which will be associated with at least one SNP). Must
#'   range between 1 and \code{ncol(list_phenos$phenos)}. Must be \code{NULL} if
#'   \code{pat} is supplied.
#' @param ind_p0 A vector of indices specifying the position of the "active"
#'   SNPs (i.e., which will be associated with at least one phenotype). Must
#'   range between 1 and \code{ncol(list_snps$snps)}. Must be \code{NULL} if
#'   \code{pat} is supplied.
#' @param vec_prob_sh If \code{block_phenos} is \code{FALSE} (default), vector
#'   of length 1 or \code{length(ind_p0)} providing the probabilities with which
#'   each active SNP will be associated with an additional active phenotype. If
#'   \code{block_phenos} is \code{TRUE}, the vector must have size between 1 and
#'   \code{ncol(list_phenos$phenos)} and gives the set of probabilities with
#'   which an active SNP is associated with an additional active phenotype is
#'   specific to each phenotypic block. Must be \code{NULL} if \code{pat} is
#'   supplied.
#' @param pat Boolean matrix of size \code{ncol(list_snps$snps)} x
#' \code{ncol(list_phenos$phenos)} which can be supplied to set the association
#'   pattern, instead of providing \code{ind_d0} and \code{ind_p0}. Must be
#'   \code{NULL} if \code{ind_d0} and \code{ind_p0} are provided.
#' @param family Distribution used to generate the phenotypes. Must be either
#' "\code{gaussian}" or "\code{binomial}" for binary phenotypes.
#' @param pve_per_snp Average proportion of phenotypic variance explained by
#'   each active SNP (for an active phenotype). Must be \code{NULL} if
#'   \code{max_tot_pve} is provided. See Details section.
#' @param max_tot_pve Maximum proportion of phenotypic variance explained by the
#'   active SNPs across all phenotypes. Must be \code{NULL} if
#'   \code{pve_per_snp} is provided. See Details section.
#' @param block_phenos Boolean for deciding whether the values in
#'   \code{vec_prob_sh} should be randomly selected and assigned differently to
#'   each block of phenotypes (if the phenotypes have no block-correlation
#'   structure, blocks are defined articially and the number of blocks
#'   corresponds to the length of the vector \code{vec_prob_sh}). Default is
#'   \code{FALSE}, no phenotypic block structure is used to create the
#'   association pattern. Not used if \code{pat} is supplied.
#' @param user_seed Seed set for reproducibility. Default is \code{NULL}, no
#'   seed set.
#'
#' @return An object of class "\code{sim_data}".
#'  \item{phenos}{Matrix containing the updated phenotypic data (whose variance
#'                is now partly explained by genetic effects).}
#'  \item{snps}{Matrix containing the original SNPs data.}
#'  \item{beta}{Matrix containing the generated effect sizes between the SNPs
#'             (rows) and phenotypes (columns).}
#'  \item{pat}{Matrix of booleans specifying the generated association pattern
#'             between the SNPs (rows) and phenotypes (columns).}
#'  \item{pve_per_snp}{Average proportion of phenotypic variance explained by
#'                     each active SNP (for an active phenotype).}
#'
#' @seealso \code{\link{convert_snps}}, \code{\link{generate_snps}},
#'   \code{\link{replicate_real_snps}}, \code{\link{convert_phenos}},
#'   \code{\link{generate_phenos}}, \code{\link{replicate_real_phenos}}
#'
#'
#' @examples
#' user_seed <- 123; set.seed(user_seed)
#' n <- 500; p <- 5000; p0 <- 200; d <- 500; d0 <- 400
#'
#' list_snps <- generate_snps(n = n, p = p)
#'
#' cor_type <- "equicorrelated"; vec_rho <- runif(100, min = 0.25, max = 0.95)
#'
#' list_phenos <- generate_phenos(n, d, cor_type = cor_type, vec_rho = vec_rho,
#'                                n_cpus = 1)
#'
#' # Gaussian phenotypes
#' dat_g <- generate_dependence(list_snps, list_phenos, ind_d0 = sample(1:d, d0),
#'                            ind_p0 = sample(1:p, p0), vec_prob_sh = 0.05,
#'                            family = "gaussian", max_tot_pve = 0.5)
#'
#' # Binary phenotypes
#' dat_b <- generate_dependence(list_snps, list_phenos, ind_d0 = sample(1:d, d0),
#'                            ind_p0 = sample(1:p, p0), vec_prob_sh = 0.05,
#'                            family = "binomial", max_tot_pve = 0.5)
#'
#' @export
#'
generate_dependence <- function(list_snps, list_phenos, ind_d0, ind_p0,
                                vec_prob_sh, pat = NULL, family = "gaussian",
                                pve_per_snp = NULL, max_tot_pve = 0.5,
                                block_phenos = FALSE, user_seed = NULL) {

  check_structure_(user_seed, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(user_seed)){
    RNGkind("L'Ecuyer-CMRG") # ensure reproducibility when using mclapply
    set.seed(user_seed)
  }

  stopifnot(family %in% c("gaussian", "binomial"))

  if (!is.null(pve_per_snp) & !is.null(max_tot_pve))
    stop("Either pve_per_snp or max_tot_pve must be NULL.")

  if (is.null(pve_per_snp) & is.null(max_tot_pve))
    warning(paste0("As both pve_per_snp or max_tot_pve were provided as NULL, the ",
                  "pve per SNP was set to its maximum value so that the total ",
                  "pve for the responses are all below 1."))

  if (!inherits(list_snps, "list_snps") & !inherits(list_snps, "sim_snps"))
    stop(paste0("The provided list_snps must be an object of class ``list_snps'' ",
               "or ``sim_snps''. \n",
               "*** You must either use the function convert_snps to obtain a ",
               "list_snps object, the function generate_snps to simulate snps ",
               "under Hardy-Weinberg equilibrium or the function replicate_real_snps ",
               "to simulate SNPs from real SNP data, by replicating their minor ",
               "allele frequencies and linkage disequilibrium structure. ***"))

  if (!inherits(list_phenos, "list_phenos") & !inherits(list_phenos, "sim_phenos"))
    stop(paste0("The provided list_phenos must be an object of class ```list_phenos'' ",
               "or `sim_phenos''. \n",
               "*** You must either use the function convert_phenos to obtain a ",
               "list_phenos object, the function generate_phenos to simulate ",
               "phenotypes from (possibly correlated) Gaussian variables or the ",
               "function replicate_real_phenos to simulate phenotypes from real ",
               "phenotypic data, by replicating their correlation structure. ***"))

  if (!is.null(pve_per_snp) & !is.null(max_tot_pve))
    stop("Either pve_per_snp or max_tot_pve must be NULL.")


  with(c(list_snps, list_phenos), {

    d <- ncol(phenos)
    n <- nrow(snps)
    p <- ncol(snps)


    check_structure_(pat, "matrix", "logical", c(p, d), null_ok = TRUE)

    if (is.null(pat)) {

      check_structure_(ind_d0, "vector", "numeric")
      check_structure_(ind_p0, "vector", "numeric")

      check_structure_(block_phenos, "vector", "logical", 1)

      if (block_phenos) {
        check_structure_(vec_prob_sh, "vector", "numeric")
        stopifnot(length(vec_prob_sh) >= 1 & length(vec_prob_sh) <= ncol(list_phenos$phenos))
      } else {
        check_structure_(vec_prob_sh, "vector", "numeric", c(1, length(ind_p0)))
      }
      check_zero_one_(vec_prob_sh)

    } else {

      stopifnot(is.null(ind_d0) & is.null(ind_p0) & is.null(vec_prob_sh))

      ind_p0 <- which(rowSums(pat)>0)
      ind_d0 <- which(colSums(pat)>0)

      if (length(ind_p0) == 0 | length(ind_d0) == 0) {
        stop("The pattern supplied must involve at least one active SNP/phenotype.")
      }
    }

    ind_d0 <- sort(unique(ind_d0))
    if (!all(ind_d0 %in% 1:d))
      stop("All indices provided in ind_d0 must be integers between 1 and d.")

    ind_p0 <- sort(unique(ind_p0))
    if (!all(ind_p0 %in% 1:p))
      stop("All indices provided in ind_p0 must be integers between 1 and p.")

    if (n != nrow(phenos))
      stop("The numbers of observations used for list_snps and for list_phenos do not match.")

    if (n == 1)
      stop("The number of observations must be greater than 1 in order to generate associations.")

    phenos_act <- phenos[, ind_d0, drop = FALSE]
    snps_act <- snps[, ind_p0, drop = FALSE]
    list_eff <- generate_eff_sizes_(d, phenos_act, snps_act, ind_d0, ind_p0, pat,
                                    vec_prob_sh, vec_maf, pve_per_snp,
                                    max_tot_pve, var_err, block_phenos,
                                    chunks_ph = ind_bl)
    with(list_eff, {
      phenos <- phenos + snps %*% beta

      if (family == "binomial")
        phenos <- ifelse(phenos > 0, 1, 0)

      list_data <- create_named_list_(phenos, snps, beta, pat, pve_per_snp)
      class(list_data) <- "sim_data"
      list_data
    })
  })
}
