# This file is part of the `echoseq` R package:
#     https://github.com/hruffieux/echoseq
#
#' Generate epigenetic marks, SNPs and phenotypes with epigenome-driven genetic
#' associations
#'
#' @param n Number of observations.
#' @param n_loci Number of loci.
#' @param mean_locus_size Mean locus size (drawn from a Poisson distribution).
#' @param p0 Minimum number of active SNPs (i.e., associated with at least one
#' phenotype).
#' @param rho_min_x Minimum autocorrelation value for blocks of SNPs in
#' linkage-disequilibrium.
#' @param rho_max_x Maximum autocorrelation value for blocks of SNPs in
#' linkage-disequilibrium.
#' @param n_modules Number of modules of phenotypes.
#' @param mean_module_size Mean module size (drawn from a Poisson distribution).
#' @param rho_min_y Minimum equicorrelation value for the phenotypes in a given
#' module. If \code{NULL}, independent SNPs simulated.
#' @param rho_max_y Minimum equicorrelation value for the phenotypes in a given
#' module. If \code{NULL}, independent SNPs simulated.
#' @param r Total number of epigenetic annotations.
#' @param r0 Number of epigenetic annotations which trigger genetic associations.
#' @param prop_act Approximate proportion of associated SNP-phenotype pairs.
#' @param max_tot_pve Maximum variance explained by the SNPs for a given
#' phenotype.
#' @param annots_vs_indep Proportion of active SNPs whose effects are triggered
#' by epigenetic marks. Default is 1, for all effects triggered by the epigenome.
#' @param min_dist Minimum distance between each pair of loci (in terms of
#' number of SNPs). Default is 0 for no distance enforced.
#' @param maf_thres Minor allele frequency threshold (applied for both supplied
#' and simulated SNPs). Default is 0.05.
#' @param max_nb_act_snps_per_locus Maximum number of active SNPs per locus.
#' Default is 3.
#' @param vec_q Exact module sizes. Either mean_module_size or vec_q must be
#' \code{NULL}. Default is \code{NULL}.
#' @param real_snp_mat Matrix of real SNPs supplied by the user. Default is
#' \code{NULL} for simulated SNPs under the Hardy-Weinberg assumption.
#' @param real_annot_mat Matrix of real epigenetic annotations supplied by the
#' user. Default is \code{NULL} for simulated binary annotations.
#' @param sd_act_beta Standard deviation of the simulated QTL effects. Either
#' sd_act_beta or max_tot_pve must be \code{NULL}. Default is \code{NULL}.
#' @param q_pres_annot_loci Quantile for selecting annotations which concern
#' most loci (i.e., at least one SNP in each locus). Should be large so enough
#' candidate active SNPs are available when annots_vs_indep is large. Default is
#' \code{NULL}.
#' @param bin_annot_freq Minimum frequency of SNPs concerned by a given
#' annotation. Default is 0.05.
#' @param candidate_modules_annots The subset of module ids where all
#' associations are triggered by annotations. The complement are the modules
#' where associations are independent of the annotations. Default is \code{NULL}
#' for all modules used as active modules. If \code{n_modules} is large, specify
#' a smaller subset of modules, as the mapping may fail otherwise.
#' @param tpois_lam_act_annots_mm Zero-truncated Poisson parameter for drawing
#' the number of active annots per module. Default is 1.
#' @param sd_act_prob Standard deviation for the effects of SNPs and annotations.
#'  Default is 1.
#' @param sd_pat Standard deviation for the randomness of the SNP-trait
#'  pattern. Default is 1.
#' @param sd_err Response error standard deviation. Default is 1.
#' @param rbeta_sh1_rr  Beta distribution shape2 parameter for the proportion of
#' responses associated with an active SNP (in a given module) rbeta_sh2_rr = 1
#' (default), so right skewed if rbeta_sh1_rr > 1.
#' @param n_cpus number of CPUs to be used. Default is 1.
#' @param maxit Maximum number of iterations for the repeat loops. Default is
#' 1e4.
#' @param module_specific Boolean specifying whether the epigenome activation is
#' module-specific or not. Default is \code{FALSE}
#' @param user_seed Seed set for reproducibility. Default is \code{NULL}, no
#' seed set.
#' @param return_patterns Boolean specifying whether the simulated SNP-phenotype
#' association pattern and active annotation variables.
#'
#' @return A list containing matrices of
#'  \item{snps}{Matrix containing the simulated or supplied SNP data.}
#'  \item{annots}{Matrix containing the simulated or supplied epiegenetic
#'  annotation data.}
#'  \item{phenos}{Matrix containing the simulated phenotypic data.}
#'  \item{pat}{If \code{return_patterns} is \code{TRUE}, simulated SNP-phenotype
#'   association pattern.}
#'  \item{beta}{If \code{return_patterns} is \code{TRUE}, simulated SNP-phenotype
#'   regression coefficients.}
#'  \item{active_annots}{If \code{return_patterns} is \code{TRUE}, active annotation
#'  variables.}
#'
#' @examples
#' user_seed <- 123; set.seed(user_seed)
#'
#' # Number of samples
#' #
#' n <- 500
#'
#' # Loci
#' #
#' n_loci <- 20
#' mean_locus_size <- 100
#' p0 <- 10
#'
#' # Modules of traits
#' #
#' n_modules <- 5
#' mean_module_size <- 50
#'
#' # Autocorrelation within loci and equicorrelation within trait modules
#' #
#' rho_min_x <- rho_min_y <- 0.5
#' rho_max_x <- rho_max_y <- 0.9
#'
#' # Annotations
#' #
#' r <- 200
#' r0 <- 10
#'
#' # Association pattern
#' #
#' prop_act <- 0.1
#' max_tot_pve <- 0.5
#'
#' list_assoc <- generate_dependence_from_annots(n, n_loci, mean_locus_size, p0,
#'                                               rho_min_x, rho_max_x,
#'                                               n_modules, mean_module_size,
#'                                               rho_min_y, rho_max_y, r, r0,
#'                                               prop_act, max_tot_pve,
#'                                               user_seed = user_seed)
#'
#' @export

generate_dependence_from_annots <- function(n,
                                            n_loci,
                                            mean_locus_size,
                                            p0,
                                            rho_min_x,
                                            rho_max_x,
                                            n_modules,
                                            mean_module_size,
                                            rho_min_y,
                                            rho_max_y,
                                            r, r0,
                                            prop_act,
                                            max_tot_pve,
                                            annots_vs_indep = 1,
                                            min_dist = 0,
                                            maf_thres = 0.05,
                                            max_nb_act_snps_per_locus = 3,
                                            vec_q = NULL,
                                            real_snp_mat = NULL,
                                            real_annot_mat = NULL,
                                            sd_act_beta = NULL,
                                            q_pres_annot_loci = NULL,
                                            bin_annot_freq = 0.05,
                                            candidate_modules_annots = NULL,
                                            tpois_lam_act_annots_mm = 1,
                                            sd_act_prob = 1,
                                            sd_pat = 1,
                                            sd_err = 1,
                                            rbeta_sh1_rr = 1,
                                            n_cpus = 1,
                                            maxit = 1e4,
                                            module_specific = FALSE,
                                            user_seed = NULL,
                                            return_patterns = FALSE) {


  check_structure_(user_seed, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(user_seed)){
    RNGkind("L'Ecuyer-CMRG") # ensure reproducibility when using mclapply
    set.seed(user_seed)
  }

  # SOME SANITY CHECKS
  #
  #
  check_natural_(n)
  check_natural_(n_loci)
  check_natural_(p0)
  check_natural_(r0)
  check_natural_(max_nb_act_snps_per_locus)
  check_natural_(n_modules)
  check_natural_(min_dist, zero_ok = TRUE)
  check_zero_one_(annots_vs_indep)

  stopifnot(annots_vs_indep > 0) # we enforce some annotation-triggered effects


  check_structure_(module_specific, "vector", "logical")

  if (module_specific) {
    stopifnot(is.null(tpois_lam_act_annots_mm))
  } else {
    check_positive_(tpois_lam_act_annots_mm)
  }

  check_structure_(real_annot_mat, "matrix", "numeric", null_ok = TRUE, na_ok = TRUE)
  check_structure_(real_snp_mat, "matrix", "numeric", null_ok = TRUE)

  if (is.null(candidate_modules_annots)) { # means: all modules considered.
    candidate_modules_annots <- 1:n_modules
    if (n_modules > 5) {
      warning(paste0("n_modules is large and the mapping may fail if ",
                     "candidate_modules_annots is NULL. If so, please specify ",
                     "a smaller subset of active modules using ",
                     "candidate_modules_annots."))
    }
  }

  stopifnot(candidate_modules_annots %in% 1:n_modules)
  candidate_modules_annots <- sort(unique(candidate_modules_annots))

  if (length(candidate_modules_annots) == n_modules) {
    stopifnot(all.equal(annots_vs_indep, 1)) # there is no module left in which associations can be simulated independently of the marks setdiff(1:n_modules, candidate_modules_annots)
  }

  if (!is.null(real_annot_mat)) {

    stopifnot(!is.null(r))

    r <- ncol(real_annot_mat)

    tb_act_annots <- NULL

    r0_list_map <- r0

  } else {

    check_natural_(r)

    stopifnot(is.null(q_pres_annot_loci))

    ind_r0 <- sort(sample(1:r, r0)) # Indices of active annots. Must be null if real_annots_mat non-NULL, and must be NULL if r0 non-NULL and vice-versa.

    # module_specific <- FALSE # for now don't give the option of module-specific activation.
    if (module_specific) {

      tb_act_annots <- set_act_annots_module_specific_(ind_r0, n_modules)

    } else {

      tb_act_annots <- set_act_annots(candidate_modules_annots,
                                      ind_r0,
                                      tpois_lam_act_annots_mm,
                                      force_ind = ind_r0)

    }

    r0_list_map <- NULL

  }


  check_natural_(min_dist,zero_ok = TRUE) # we should have an assertion for checking that sufficient number of SNPs to enforce this distance.

  if (is.null(real_snp_mat)) {

    check_zero_one_(rho_min_x)
    check_zero_one_(rho_max_x)

    stopifnot(all.equal(min_dist, 0))

  } else {

    stopifnot(is.null(rho_min_x) & is.null(rho_max_x))

    stopifnot(n <= nrow(real_snp_mat))

    stopifnot(real_snp_mat %in% 0:2)

  }

  check_positive_(mean_locus_size)

  stopifnot(xor(mean_module_size, vec_q))

  if (!is.null(mean_module_size)) {
    check_positive_(mean_module_size)
  }
  if (!is.null(vec_q)) {
    check_natural_(vec_q)
    stopifnot(all.equal(length(vec_q), n_modules))
  }

  check_positive_(rbeta_sh1_rr)

  check_zero_one_(annots_vs_indep)
  check_zero_one_(maf_thres)

  if (!is.null(r0_list_map)) {

    stopifnot(r0_list_map %in% 1:r) # (necessary but not sufficient, as some non-binary annots may be dropped)
    stopifnot(is.null(tb_act_annots))

  }

  if (!is.null(tb_act_annots)) {

    stopifnot(names(tb_act_annots) == c("module", "act_annot"))
    stopifnot(tb_act_annots$module %in% 1:n_modules)
    stopifnot(tb_act_annots$act_annot %in% 1:r) # (necessary but not sufficient, as some non-binary annots may be dropped)
    stopifnot(is.null(real_annot_mat)) # the chosen annots will depend on their locus coverage
    stopifnot(is.null(r0_list_map))

    tpois_lam_act_annots_mm <- NULL

  }


  check_zero_one_(prop_act)
  check_positive_(sd_act_prob)
  check_positive_(sd_err)
  check_zero_one_(rho_min_y)
  check_zero_one_(rho_max_y)

  if (!xor(is.null(max_tot_pve), is.null(sd_act_beta))) {

    stop("Either sd_act_beta or max_tot_pve must be provided.")

  }

  if (!is.null(sd_act_beta)) {
    check_positive_(sd_act_beta)
  } else {
    check_zero_one_(max_tot_pve)
  }

  check_natural_(maxit)
  check_natural_(n_cpus)
  stopifnot(n_cpus <= parallel::detectCores())

  it <- 1
  repeat { # regenerate if to catch !(all(colSums(map_act_snps_modules[, candidate_modules, drop = FALSE]) > 0))
    # i.e. case where not all modules have at least one active SNPs for the annots
    # can occur when !is.null(real_annot_mat)

  list_map <- tryCatch({
      set_locus_and_module_pattern_from_annots_(n,
                                               real_annot_mat,
                                               q_pres_annot_loci,
                                               r,
                                               r0_list_map, # if !bool_real_annots, tb_act_annots is provided so no need to provide r0 (must be NULL)
                                               tb_act_annots,
                                               tpois_lam_act_annots_mm,
                                               bin_annot_freq,
                                               real_snp_mat,
                                               rho_min_x,
                                               rho_max_x,
                                               n_loci,
                                               mean_locus_size,
                                               min_dist,
                                               p0,
                                               max_nb_act_snps_per_locus,
                                               maf_thres,
                                               n_modules,
                                               mean_module_size,
                                               vec_q,
                                               candidate_modules_annots, # all modules are considered equally for annot-triggered and indep associations
                                               annots_vs_indep,
                                               rbeta_sh1_rr,
                                               n_cpus)
    }, error=function(e){})

    it <- it + 1

    if (!is.null(list_map) | it > maxit) {
      break
    }
  }


  stopifnot(!is.null(list_map))

  annots <- list_map$V
  snps <- list_map$X

  list_Y <- generate_phenos_from_annots_(list_map,
                                        prop_act,
                                        sd_act_prob,
                                        sd_pat,
                                        sd_act_beta,
                                        sd_err,
                                        max_tot_pve,
                                        rho_min_y,
                                        rho_max_y,
                                        n_cpus)

  phenos <- list_Y$Y

  if (return_patterns) {
    pat <- list_Y$pat
    beta <- list_Y$beta
    active_annots <- list_map$list_map_annots$tb_act_annots
  } else {
    pat <- beta <- active_annots <- NULL
  }
  create_named_list_(snps, annots, phenos, pat, beta, active_annots)
}

set_locus_and_module_pattern_from_annots_ <- function(n, # sample size

                                                     # ANNOTS
                                                     #
                                                     real_annot_mat,          # if NULL, binary annots simulated
                                                     q_pres_annot_loci,       # if !is.null(real_annot_mat), quantile for selecting annots which concern most loci (i.e., at least one SNP in each locus) - should be large so enough candidate active SNPs are available when annots_vs_indep is large
                                                     r,                      # total number of annots, can be NULL if !is.null(real_annot_mat) in which case all annots are used
                                                     r0,                     # total number of candidate active annots (if NULL, all r annots; if small, the same annots will tend to be active for several modules); if 1, mostly 1-2 active annots in each module, 5, mean number of annots approx 5
                                                     tb_act_annots,           # table mapping the modules with their corresponding active annots. Must be obtained using the set_act_annots function. tpois_lam_act_annots_mm must be NULL as unused. tb_act_annots must be NULL if !is.null(real_annot_mat), and r0 and tb_act_annots cannots both be non-NULL.
                                                     tpois_lam_act_annots_mm, # zero-truncated Poisson parameter for drawing the number of active annots per module
                                                     bin_annot_freq,          # minimal frequency of SNPs concerned by a given annotation

                                                     # PREDICTORS (SNPS)
                                                     #
                                                     real_snp_mat,              # if NULL, SNPs simulated under the Hardy-Weinberg assumption
                                                     rho_min_x,                 # if is.null(real_snp_mat), minimum correlation for blocks of simulated autocorrelated SNPs
                                                     rho_max_x,                 # if is.null(real_snp_mat), minimum correlation for blocks of simulated autocorrelated SNPs
                                                     n_loci,                    # number of loci
                                                     mean_locus_size,           # mean locus size (drawn from a Poisson distribution)
                                                     min_dist,                  # minimum number of SNPs separating each pair of loci
                                                     p0,                        # minimum number of active SNPs (the actual number will also depend on the prop_act argument of the generate_phenos_from_annots_ function)
                                                     max_nb_act_snps_per_locus, # maximum number of active SNPs allowed per locus (typically 4-5?)
                                                     maf_thres,                 # minor allele frequency threshold (applied for both real and simulated SNPs)

                                                     # (SIMULATED) RESPONSES
                                                     #
                                                     n_modules,               # number of modules
                                                     mean_module_size,        # mean module size (drawn from a Poisson distribution)
                                                     vec_q,                   # exact module sizes provided. Either mean_module_size or vec_q must be NULL.
                                                     candidate_modules_annots, # the subset of module ids where all associations are triggered by annots.
                                                     # the complement, setdiff(1:n_modules, candidate_modules_annots) are the modules where associations are independent of annots
                                                     #
                                                     annots_vs_indep, # proportion of hotspots driven/triggered by the annots
                                                     rbeta_sh1_rr,   # Beta distribution shape2 parameter for the proportion of responses associated with an active SNP (in a given module)
                                                     # rbeta_sh2_rr = 1, so right skewed if rbeta_sh1_rr > 1

                                                     # OTHER SETTINGS
                                                     #
                                                     n_cpus = 1, # number of CPUs to be used for the VBEM algorithm (if n_modules > 1)
                                                     maxit = 1e4 # maximum number of iterations for the repeat loops
) {  #


  # BUILD LOCI PATTERN
  #
  #
  # Poisson distribution with user-specified mean locus size,
  # but making sure that the minimum size generated is larger than the
  # user-specified maximum number of active SNPs per locus
  # (usually verified if the mean locus size is reasonably large and the maximum
  # number of active SNPs per locus reasonably small).
  #
  it <- 1

  repeat {

    vec_p <- rpois(n_loci, lambda = mean_locus_size) # E(X) = lambda

    it <- it + 1

    if (min(vec_p) >= max_nb_act_snps_per_locus | it > maxit) {
      break
    }

  }

  # Stop if condition still not verified after maxit iterations
  #
  stopifnot(min(vec_p) >= max_nb_act_snps_per_locus)

  p <- sum(vec_p)
  n_loci <- length(vec_p)

  pos_loci <- c(1, cumsum(vec_p[-n_loci]) + 1)


  # List of loci, positions etc
  #
  list_loci <- set_blocks_(p, pos_loci, n_cpus = n_cpus)


  # SOME FURTHER SANITY CHECKS
  #
  stopifnot(list_loci$n_var_blocks == p)
  stopifnot(list_loci$n_bl== n_loci)

  if (is.null(real_annot_mat)) {

    stopifnot(is.null(q_pres_annot_loci))

  } else if (is.null(real_snp_mat)){

    stop("If real annots supplied, the corresponding real SNP data must be supplied as well.")

  } else { # both real_annot_mat and real_snp_mat non-NULL

    stopifnot(ncol(real_annot_mat) >= r) # (necessary but not sufficient, as some non-binary annots may be dropped)
    check_zero_one_(q_pres_annot_loci)

    stopifnot(r*(1-q_pres_annot_loci) > 1) # at least one candidate annot to draw from (necessary but not sufficient)

    snps_with_annots <- intersect(colnames(real_snp_mat), rownames(real_annot_mat))

    stopifnot(length(snps_with_annots) >= p) # necessary condition for drawing X and V (not sufficient)

  }

  # BUILD MODULE PATTERN
  #
  #
  # Poisson distribution with user-specified mean module size,
  # but making sure that the minimum size generated is larger than 10, so there
  # is enough information to estimate the module-specific hotspot propensity
  # variance (usually verified if the mean module size is reasonably large, e.g.,
  # mean_module_size = 40)
  #
  it <- 1

  if (is.null(vec_q)) { # exact module sizes not provided, only mean module size
    repeat {

      vec_q <- rpois(n_modules, lambda = mean_module_size)

      it <- it + 1

      if (min(vec_q) >= 10 | it > maxit) {
        break
      }

    }
  }

  # stop if condition still not verified after maxit iterations
  #
  stopifnot(min(vec_q) >= 10)

  q <- sum(vec_q)
  n_modules <- length(vec_q)

  pos_modules <- c(1, cumsum(vec_q[-n_modules]) + 1)


  # list to be given to the epispot function. Partition into modules. No partition of the SNPs.
  #
  list_partition <- set_blocks_(c(p, q), list(1, pos_modules), n_cpus = n_cpus) # n_cpus used for the VBEM algo

  if (n_modules > 1) {
    stopifnot(list_partition$bl_y$n_var_blocks == q)
    stopifnot(list_partition$bl_y$n_bl == n_modules)
  }

  if (!is.null(real_annot_mat)) {

    # BUILD SNP MATRIX AND CORRESPONDING ANNOT MATRIX IF BOTH FROM REAL DATA
    #
    #
    # Exclude non-binary annots (such as the minimum distance to TSS) if present
    # (binary annot assumption needed below when creating the active SNP pattern).
    #
    bool_bin <- apply(real_annot_mat, 2, function(annot) all(annot[!is.na(annot)] %in% 0:1))
    real_annot_mat <- real_annot_mat[, bool_bin, drop = FALSE]

    it <- 1

    repeat {

      # build matrix of SNPs
      #
      X <- generate_snps_from_loci_(n,
                                   list_loci,
                                   real_snp_mat = real_snp_mat,
                                   maf_thres = maf_thres,
                                   min_dist = min_dist) # 150 SNPs corresponds to a median distance of approx 1Mb for chr1 of our monocyte data (see save_chr_1_unstim.R)


      # build matrix of annots
      #
      V <- real_annot_mat[match(colnames(X), rownames(real_annot_mat)), , drop = FALSE] # if a SNP is absent from in real_annot_mat, NA row produced

      # remove annots whose frequency is too low to be informative
      #
      V <- rm_bin_annot_freq_(V, bin_annot_freq)

      if (r < ncol(V)) {
        V <- V[, sort(sample(1:ncol(V), size = r)), drop = FALSE]
      }

      it <- it + 1

      if (sum(is.na(V)) == 0 | it > maxit) { # make sure that all SNPs have all annots available
        break
      }
    }

    if (sum(is.na(V))>0) {
      stop("No complete annot matrix for the selected loci.")
    }

    stopifnot(all.equal(rownames(V), colnames(X)))

    if (r > ncol(V)) {

      r <- ncol(V)

      stopifnot(r*(1-q_pres_annot_loci) > 1) # at least one candidate annot to draw from

      if (!is.null(r0)) {
        stopifnot(r0 %in% 1:r)
      }

      if (!is.null(tb_act_annots)) {
        stopifnot(tb_act_annots$act_annot %in% 1:r)
      }

    }
  }

  if (annots_vs_indep > 0) { # at least one active SNP will be triggered by the annots

    # SELECT ACTIVE ANNOTS
    #
    if (!is.null(real_annot_mat)) {

      # for each annot, proportion of loci in which it concerns at least one SNP
      #
      prop_V_in_loci <- colMeans(apply(V, 2, function(vv) {

        sapply(levels(list_loci$vec_fac_bl), function(ll) sum(vv[list_loci$vec_fac_bl == ll])>0)

      }))

      # threshold on the proportion of loci with at least one SNP concerned by the annot under consideration
      # (needs to be high enough so that an active annot can concern sufficiently many SNPs - and hence one has sufficient info to estimate its effect)
      #
      thres_pr_locus_has_a_snp_in_a_annot <- quantile(prop_V_in_loci, probs = q_pres_annot_loci)

      # the candidate annots are those which have, in each locus,
      # a probability > thres_pr_locus_has_a_snp_in_a_annot to concern at least one SNP
      #
      candidate_act_annots <- which(prop_V_in_loci > thres_pr_locus_has_a_snp_in_a_annot)

    } else {

      if (is.null(tb_act_annots)) {
        candidate_act_annots <- 1:r # simulated annot matrix --> all annots can be taken as candidate as V will be generated accordingly late
      }

      V <- NULL # V to be simulated later, based on the annot/snp association pattern

    }

    if (is.null(tb_act_annots)) {

      if (!is.null(r0) && r0 < length(candidate_act_annots)) { # choose among a restricted subset of r0 annots

        candidate_act_annots <- sort(sample(candidate_act_annots, r0))

      }

      tb_act_annots <- set_act_annots(candidate_modules_annots,
                                      candidate_act_annots,
                                      tpois_lam_act_annots_mm)

    }

    # START BY CHOOSING SNPS WHOSE ACTIVITY IS DRIVEN BY THE ANNOTS (only if annots_vs_indep > 0)
    #
    p0_annots <- ceiling(p0 * annots_vs_indep) # number of active SNPs driven by the annots

    list_map_annots <- choose_act_snps_("annots",
                                        q, # needed for the case n_modules = 1, and cannots be retrieved from list_partition
                                        p0_annots,
                                        list_loci,
                                        list_partition,
                                        max_nb_act_snps_per_locus,
                                        rbeta_sh1_rr,
                                        n_modules,
                                        candidate_modules = candidate_modules_annots,
                                        excluded_act_snps = NULL,
                                        tb_act_annots = tb_act_annots,
                                        user_locus_ids = NULL,
                                        V = V,
                                        maxit = maxit)

    list_map_annots$tb_act_annots <- tb_act_annots

    # THEN CHOOSE SNPS WHOSE ACTIVITY IS NOT DRIVEN BY THE ANNOTS
    #
    locus_ids_annots <- sort_loci_by_number_of_act_snps_(list_map_annots$tb_act_snps, list_loci$n_bl)

    p0_indep <- p0 - p0_annots
    excluded_act_snps <- list_map_annots$tb_act_snps$act_snp
    user_locus_ids <- locus_ids_annots

  } else {

    p0_indep <- p0
    excluded_act_snps <- NULL
    user_locus_ids <- NULL

    list_map_annots <- NULL

  }

  if (annots_vs_indep < 1) {

    candidate_modules_indep <- setdiff(1:n_modules, candidate_modules_annots)

    list_map_indep <- choose_act_snps_("indep",
                                       q,
                                       p0_indep,
                                       list_loci,
                                       list_partition,
                                       max_nb_act_snps_per_locus,
                                       rbeta_sh1_rr,
                                       n_modules,
                                       candidate_modules = candidate_modules_indep,
                                       excluded_act_snps = excluded_act_snps,
                                       tb_act_annots = list_map_annots$tb_act_annots,
                                       user_locus_ids = user_locus_ids,
                                       V = NULL,
                                       maxit = maxit)

  } else {

    list_map_indep <- NULL

  }

  if (is.null(real_annot_mat)) {

    V <- generate_binary_annots_(p, r, list_map_annots)

    X <- generate_snps_from_loci_(n,
                                 list_loci,
                                 real_snp_mat = real_snp_mat,
                                 rho_min = rho_min_x,
                                 rho_max = rho_max_x,
                                 maf_thres = maf_thres,
                                 min_dist = min_dist) # minimal distance between each locus (in # of SNPs)

  }


  list_map <- create_named_list_(X, V, q,
                                 pos_loci,
                                 pos_modules,
                                 list_partition,
                                 list_loci,
                                 list_map_annots,
                                 list_map_indep)

  class(list_map) <- "sim_map"

  list_map

}




set_act_annots_module_specific_ <- function(ind_r0, n_modules) {

  r0 <- length(ind_r0)

  # module-specific annots # we force the modules to be concerned by different annots
  #
  stopifnot(r0 >= n_modules)

  ind_r0 <- sort(ind_r0)

  # Assign a different annot to each module
  #
  tb_act_annots <- as.data.frame(cbind(1:n_modules, sample(ind_r0, size = n_modules)))

  # Assign the remaining active annots
  #
  if (r0 > n_modules) {
    tb_act_annots <- rbind(tb_act_annots,
                           cbind(sample(1:n_modules, r0 - n_modules, replace = TRUE),
                                 setdiff(ind_r0, tb_act_annots[,2])))
  }

  tb_act_annots <- as.data.frame(tb_act_annots)
  names(tb_act_annots) <- c("module", "act_annot")

  tb_act_annots <- tb_act_annots[order(tb_act_annots$act_annot),]
  tb_act_annots <- tb_act_annots[order(tb_act_annots$module),]

  tb_act_annots

}

set_act_annots <- function(candidate_modules,
                           candidate_ind,      # these annot ids COULD be assigned as active annot
                           tpois_lam_act_annots_mm,
                           force_ind = NULL) { # these annot ids WILL be assigned to at least one module (we make sure that they are)


  if (!is.null(force_ind)) {

    place_all_force_ind <- sample(candidate_modules, size = length(force_ind), replace = TRUE)

  }

  tb_act_annots <- lapply(candidate_modules, function(mm) {

    if (length(candidate_ind) > 1) {

      if (!is.null(force_ind) && any(place_all_force_ind %in% mm)) {

        act_annots_mm <- force_ind[place_all_force_ind %in% mm]

      } else {

        act_annots_mm <- NULL

      }

      # Zero-truncated Poisson distribution
      #
      act_annots_mm <- sort(unique(c(act_annots_mm, sort(sample(candidate_ind,
                                                                size = min(
                                                                  rtpois(1, lambda = tpois_lam_act_annots_mm, a = 0),
                                                                  length(candidate_ind)
                                                                )
      )
      ))))

    } else { # single candidate annot

      act_annots_mm <- candidate_ind

    }

    cbind(rep(mm, length(act_annots_mm)), act_annots_mm)

  })

  tb_act_annots <- as.data.frame(do.call(rbind, tb_act_annots))
  names(tb_act_annots) <- c("module", "act_annot")

  tb_act_annots <- tb_act_annots[order(tb_act_annots$act_annot),]
  tb_act_annots <- tb_act_annots[order(tb_act_annots$module),]

  tb_act_annots

}



# Sort the loci based on the number of already assigned active SNPs in these loci:
# start with the loci with few active SNPs and finish with those with many active SNPs
#
sort_loci_by_number_of_act_snps_ <- function(tb_act_snps, n_loci) {

  # the "sample()" is to shuffle the loci within groups of loci having a given number of active SNPs
  #
  tb_act_snps <- as.data.frame(tb_act_snps)
  names(tb_act_snps) <- c("locus", "act_snp")


  nb_act_snps_per_locus <- tb_act_snps[sample(1:nrow(tb_act_snps)), ] %>% group_by(.data$locus) %>% mutate(nb_act_snps = length(.data$act_snp))
  nb_act_snps_per_locus <- unique(nb_act_snps_per_locus[,c("locus", "nb_act_snps")])
  nb_act_snps_per_locus <- nb_act_snps_per_locus[order(nb_act_snps_per_locus$nb_act_snps), ]

  # the "sample()" is to shuffle the loci within groups of loci having no active SNP
  #
  c(setdiff(sample(1:n_loci), nb_act_snps_per_locus$locus), # first the loci which do not have any active SNP yet
    nb_act_snps_per_locus$locus) # then the loci with only 1 active SNPs, then those with 2, up to max_nb_act_snps_per_locus

}




choose_act_snps_ <- function(type,                     # flag "annots" or "indep" depending on whether the chosen active SNPs are driven by the annots or not
                             q,                         # total number of responses
                             p0_sub,                    # number of active SNPs to choose
                             list_loci,                 # locus partition
                             list_partition,            # module partition
                             max_nb_act_snps_per_locus, # maximum number of active SNPs per locus allowed
                             rbeta_sh1_rr,              # Beta distribution shape2 parameter for the proportion of responses associated with an active SNP (in a given module), rbeta_sh1_rr = 1, so left skewed
                             n_modules,
                             candidate_modules,         # modules to be considered. Must match the modules in tb_act_annots if type == "annots"
                             excluded_act_snps = NULL,  # SNP ids excluded from the candidate active SNPs
                             tb_act_annots = NULL,       # dataframe of active annots for each module
                             user_locus_ids = NULL,     # locus ids from which to select the active SNPs (by order of appearance)
                             V = NULL,                  # matrix of annots
                             maxit = 1e4                # maximum number of iterations to be used in the repeat loops
) {


  stopifnot(type %in% c("annots", "indep"))

  stopifnot(candidate_modules %in% 1:n_modules)

  if (type == "annots") {

    stopifnot(!is.null(tb_act_annots))
    candidate_modules <- sort(unique(candidate_modules))

    stopifnot(all.equal(candidate_modules, sort(unique(tb_act_annots$module))))

  }

  check_natural_(q)
  check_natural_(p0_sub)
  check_natural_(max_nb_act_snps_per_locus)
  check_natural_(maxit)
  check_positive_(rbeta_sh1_rr)

  if (!is.null(excluded_act_snps)) {
    check_natural_(excluded_act_snps)
  }

  if (!is.null(user_locus_ids)) {
    check_natural_(user_locus_ids)
  }


  if (!is.null(V)) {

    if (is.null(tb_act_annots)) {

      stop("Object tb_act_annots must be non-NULL if V is non-NULL.")

    }

    vec_act_annots <- sort(unique(tb_act_annots$act_annot))
    stopifnot(all(V %in% 0:1))

    stopifnot(sum(rowSums(V[, vec_act_annots, drop = FALSE]) > 0) >= p0_sub) # the number of SNPs having at least one annot
    # must be larger than the number of active SNPs driven by annots
  }

  n_loci <- list_loci$n_bl
  n_modules <- ifelse(is.null(list_partition$bl_y), 1, list_partition$bl_y$n_bl)

  tb_act_snps <-  map_act_snps_modules <- map_act_snps_responses <- NULL
  p0_current <- 0
  it <- 1

  repeat{

    if (it > 1) { # If we still have p0_current < p0_sub after walking through all the loci, do a second (third, fourth, etc if needed) round

      # Exclude from the list of candidate active SNPs, those that have already been selected in the previous iteration(s)
      #
      excluded_act_snps <- unique(c(excluded_act_snps, tb_act_snps[,2]))

      # Choose the order in which we will walk through the loci based on the number of already assigned active SNPs in these loci: start with the loci with few active SNPs
      #
      locus_ids <-  sort_loci_by_number_of_act_snps_(tb_act_snps, n_loci)

      if (!is.null(user_locus_ids)) {

        stopifnot(length(user_locus_ids) == n_loci)

        # Combine the info from the two locus_ids vectors in order to start walking through the loci with the fewest already assigned active SNPs, overall
        # i.e., sum the ranks of the loci in both user_locus_ids and locus_ids and sort the loci from that with the smaller summed rank to the largest summed rank
        #
        locus_ids <- order(order(user_locus_ids) + order(locus_ids))

      }

    } else if (is.null(user_locus_ids)) {

      locus_ids <- sample(1:n_loci)

    } else {

      locus_ids <- user_locus_ids

    }

    # Assign active SNPs locus by locus
    #
    for (ll in locus_ids) {

      if (p0_current < p0_sub) {

        if (!is.null(V)) { # Real matrix of annots used and active SNPs triggered by annots (list_map_annots)

          # The candidate SNPs in the locus are those which are concerned by at least one active annot
          #
          V_ll_act_annots <- V[list_loci$vec_fac_bl == ll, vec_act_annots, drop = FALSE]
          V_ll_act_annots_nonzero <- V_ll_act_annots[rowSums(V_ll_act_annots)>0,, drop = FALSE]

          candidate_act_snps_ll <- which(rownames(V) %in% rownames(V_ll_act_annots_nonzero))

        } else { # Simulated matrix of annots used (will be generated below, based on the generated active SNP pattern --> so no constraint)
          # and/or active SNPs not triggered by annots (list_map_indep)

          candidate_act_snps_ll <- which(list_loci$vec_fac_bl == ll)

        }

        if (!is.null(excluded_act_snps)) {

          # Number of SNPs that are already active in the locus
          #
          nb_act_snps_excluded_ll <- length(intersect(candidate_act_snps_ll, excluded_act_snps))
          candidate_act_snps_ll <- setdiff(candidate_act_snps_ll, excluded_act_snps)

        } else {

          nb_act_snps_excluded_ll <- 0

        }


        if (length(candidate_act_snps_ll) > 0 & # At least one candidate active SNP in the locus
            max_nb_act_snps_per_locus - nb_act_snps_excluded_ll > 0) { # The maximum number of active SNPs per locus has not been reached yet

          if (length(candidate_act_snps_ll) == 1) {

            p0_ll <- 1

            act_snps_ll <- candidate_act_snps_ll

          } else {

            p0_ll <- sample(1:min(c(max_nb_act_snps_per_locus - nb_act_snps_excluded_ll, # substract the number of SNPs already active in the locus from the total number of active SNPs allowed per locus
                                    length(candidate_act_snps_ll),
                                    p0_sub - p0_current)),
                            size = 1)

            act_snps_ll <- sort(sample(candidate_act_snps_ll,
                                       size = p0_ll)) # choose up to max_nb_act_snps_per_locus active SNPs for each locus ll
            # among the candidate active SNPs in ll
          }

          tb_act_snps <- rbind(tb_act_snps, cbind(rep(ll, length(act_snps_ll)), act_snps_ll))

          for(ss in act_snps_ll) {


            # ASSIGN ACTIVE MODULES TO EACH ACTIVE SNP SS
            #
            if (is.null(V)) { # Simulated matrix of annots and/or active SNPs not triggered by annots (list_map_indep)

              candidate_act_mm_for_ss <- candidate_modules

            } else { # Real matrix of annots and active SNPs triggered by annots (list_map_annots)

              # The candidate module(s) for SNP ss are those whose active annot(s) concern ss
              #
              ss_name <- rownames(V)[ss]
              candidate_act_mm_for_ss <- sort(unique(tb_act_annots$module[tb_act_annots$act_annot %in% vec_act_annots[which(V_ll_act_annots_nonzero[ss_name, ] > 0)]]))

            }

            if (length(candidate_act_mm_for_ss) > 1 & !is.null(tb_act_annots)) { # Choose the number of modules associated with a given SNP
              # based on the number of modules concerned by a same active annot

              pr_module_assoc_ss <- sample(table(tb_act_annots$act_annot) / length(unique(tb_act_annots$module)), 1)

              mm_assoc_ss <- sort(sample(candidate_act_mm_for_ss,
                                         size = ceiling(pr_module_assoc_ss * length(candidate_act_mm_for_ss))))

            } else { # Only enters here when annots_vs_indep = 0 (and test ranges of annots_vs_indep only for n_modules = 1, so ok)

              mm_assoc_ss <- candidate_act_mm_for_ss

            }

            map_act_ss_modules <- rep(FALSE, n_modules)
            map_act_ss_modules[mm_assoc_ss] <- TRUE
            map_act_snps_modules <- rbind(map_act_snps_modules, map_act_ss_modules)


            # ASSIGN RESPONSE(S) TO EACH ACTIVE SNP SS (WITHIN EACH MODULE ASSOCIATED WITH SS)
            #
            rr_assoc_ss <- rep(NA, n_modules)

            rr_assoc_ss[mm_assoc_ss] <- sapply(mm_assoc_ss, function(mm) {

              # Its probability to be associated with each response in the active module mm, right-skewed distribution if rbeta_sh1_rr > 1
              #
              pr_within_module_assoc_ss <- rbeta(1, shape1 = rbeta_sh1_rr, shape2 = 1)

              if (n_modules > 1) {

                # Retrieve all the responses within module mm
                #
                rr_in_mm <- which(list_partition$bl_y$vec_fac_bl == mm)
                n_rr_in_mm <- length(rr_in_mm)

              } else { # list_partition$bl_y does not exist

                rr_in_mm <- 1:q
                n_rr_in_mm <- q

              }

              list(sort(sample(rr_in_mm, size = ceiling(n_rr_in_mm * pr_within_module_assoc_ss))))

            } )

            map_act_snps_responses <- rbind(map_act_snps_responses, rr_assoc_ss)

          }

          # Bump the counter for the number of active SNPs currently assigned
          #
          p0_current <- p0_current + p0_ll

        }

      }

    }

    it <- it + 1

    if (p0_current == p0_sub | it > maxit) { # the total number of active SNPs to be assigned has been reach.

      if ((type == "annots" & all(colSums(map_act_snps_modules[, candidate_modules, drop = FALSE]) > 0)) |  # make sure that all modules are concerned by at least one active SNP,
          # as each module is assigned to at least one active annot
          # (in the real_annots_mat case, does not mean that all active annots
          # are concerned by at least one SNP, unfortunately...
          # But this is ensured when the annots are simulated, even when using real SNPs).
          type == "indep" |
          it > maxit) {

        break

      }

    }

  }

  tb_act_snps <- as.data.frame(tb_act_snps)
  names(tb_act_snps) <- c("locus", "act_snp")

  stopifnot(p0_current == p0_sub)

  if (type == "annots") {

    stopifnot(all(colSums(map_act_snps_modules[, candidate_modules, drop = FALSE]) > 0)) # each module must be associated with at least one active SNP (because to each module is attributed at least one active annot)

  }

  rownames(map_act_snps_modules) <- rownames(map_act_snps_responses) <- paste0("act_snp_", tb_act_snps$act_snp)
  colnames(map_act_snps_modules) <- colnames(map_act_snps_responses) <- paste0("module_", 1:n_modules)

  stopifnot(all(rowSums(map_act_snps_modules) > 0))                              # each active SNP must be associated with at least one module
  stopifnot(identical(!is.na(map_act_snps_responses), map_act_snps_modules))     # each pair of active SNP / active module must have at least one SNP-response association
  # and there must be no association for inactive pairs.

  stopifnot(nrow(tb_act_snps) == p0_sub)

  # Reorder dataframes so they look nicer
  #
  ord_snps <- order(tb_act_snps$act_snp)

  map_act_snps_modules <- map_act_snps_modules[ord_snps, , drop = FALSE]
  map_act_snps_responses<- map_act_snps_responses[ord_snps, , drop = FALSE]

  tb_act_snps <- tb_act_snps[ord_snps, , drop = FALSE]
  rownames(tb_act_snps) <- NULL

  create_named_list_(tb_act_snps,
                     map_act_snps_modules,
                     map_act_snps_responses)

}

generate_phenos_from_annots_ <- function(list_map,           # object supplied by the set_locus_and_module_pattern_from_annots_ function
                                        prop_act,           # approximate proportion of non-zero predictor-response effects
                                        sd_act_prob,        # standard deviation for the effects of hotspot propensities and annots (within probit link - Gaussian effects)
                                        sd_pat,             # standard deviation for the randomness of the predictor-response pattern
                                        sd_act_beta,        # standard deviation for the regression effects of between predictors and responses (Gaussian effects)
                                        sd_err,             # response error standard deviation
                                        max_tot_pve = NULL, # maximum proportion of response variance explained by the SNPs, for each response
                                        rho_min = 0,        # minimum block equicorrelation level for the responses
                                        rho_max = 0.5,       # maximum block equicorrelation level for the responses
                                        n_cpus = 1
) {


  if (!inherits(list_map, "sim_map"))
    stop("The provided list_map must be an object of class ``sim_map''. \n")

  check_zero_one_(prop_act)
  check_positive_(sd_act_prob)
  check_positive_(sd_err)
  check_zero_one_(rho_min)
  check_zero_one_(rho_max)

  if (!xor(is.null(max_tot_pve), is.null(sd_act_beta))) {

    stop("Either sd_act_beta or max_tot_pve must be provided.")

  }

  if (!is.null(sd_act_beta)) {
    check_positive_(sd_act_beta)
  } else {
    check_zero_one_(max_tot_pve)
  }


  V <- list_map$V
  stopifnot(all(V %in% 0:1))

  X <- list_map$X

  n <- nrow(X)
  p <- ncol(X)
  r <- ncol(V)
  q <- list_map$q

  stopifnot(nrow(V) == p)

  if (!is.null(list_map$list_partition$bl_y)) { # more than a single module

    list_modules <- list_map$list_partition$bl_y
    n_modules <- list_modules$n_bl
    stopifnot(list_modules$n_var_blocks == q)

  } else {

    n_modules <- 1

  }

  # Gather patterns for effects driven by the annots
  #
  if (!is.null(list_map$list_map_annots)) {

    tb_act_annots <- list_map$list_map_annots$tb_act_annots
    tb_act_snps_annots <- list_map$list_map_annots$tb_act_snps
    map_act_snps_responses_annots <- list_map$list_map_annots$map_act_snps_responses
    map_act_snps_modules_annots <- list_map$list_map_annots$map_act_snps_modules

  }

  # Gather patterns for effects not driven by the annots
  #
  if (!is.null(list_map$list_map_indep)) {

    tb_act_snps_indep <- list_map$list_map_indep$tb_act_snps
    map_act_snps_responses_indep <- list_map$list_map_indep$map_act_snps_responses
    map_act_snps_modules_indep <- list_map$list_map_indep$map_act_snps_modules

  }


  # Set the sparsity level
  #
  n0 <- -2.5
  W_arg <- matrix(n0, nrow = p, ncol = q)

  for (mm in 1:n_modules) {

    # ADD EFFECTS DRIVEN BY THE ANNOTS
    #
    if (!is.null(list_map$list_map_annots)) {

      xi_mm <- rep(0, r)
      act_annots_mm <- tb_act_annots$act_annot[tb_act_annots$module == mm]

      # log normal: positive effect --> SNP within an annot favours its potential associations with transcripts
      #
      xi_mm[act_annots_mm] <-  rlnorm(length(act_annots_mm), meanlog = 1, sdlog = sd_act_prob)

      for (ss in which(map_act_snps_modules_annots[, mm])) {

        act_snp_ss <- tb_act_snps_annots$act_snp[ss]
        act_resp_ss_mm <- unlist(map_act_snps_responses_annots[ss, mm])

        # Divide the effects by the number of active annots that concern SNP ss
        # in order to prevent SNPs concerned by many annots to dominate over the
        # others (SNPs concerned by many annots should not necessarily be larger
        # hotspots than the others)
        # sum(V[act_snp_ss, act_annots_mm]): since V is binary OK
        #
        W_arg[act_snp_ss, act_resp_ss_mm] <- W_arg[act_snp_ss, act_resp_ss_mm] + as.numeric(V[act_snp_ss, , drop = FALSE] %*% xi_mm) / sum(V[act_snp_ss, act_annots_mm])

      }

    }


    # ADD EFFECTS NOT DRIVEN BY THE ANNOTS
    #
    if (!is.null(list_map$list_map_indep)) {

      for (ss in which(map_act_snps_modules_indep[, mm])) {

        act_snp_ss <- tb_act_snps_indep$act_snp[ss]
        act_resp_ss_mm <- unlist(map_act_snps_responses_indep[ss, mm])

        # SNP-specific effect
        #
        theta_ss_mm <- rlnorm(1, meanlog = 1, sdlog = sd_act_prob)

        W_arg[act_snp_ss, act_resp_ss_mm] <- W_arg[act_snp_ss, act_resp_ss_mm] + theta_ss_mm

      }

    }

  }

  W <- apply(W_arg, 2, function(W_arg_col) rnorm(p, mean = W_arg_col, sd = sd_pat))

  # prop_act is the desired proportion of associations
  #
  pat <- beta <- matrix(0, nrow = p, ncol = q)
  pat[W >= quantile(W, probs = 1 - prop_act)] <- 1

  # Make sure that at least one association for each pair of active SNP/active module
  #
  add_missing_activations <- function(pat,
                                      list_map,
                                      n_modules) {

    if (!is.null(list_map)) {

      for (mm in 1:n_modules) {

        for (ss in which(list_map$map_act_snps_modules[, mm])) {

          act_snp_ss <- list_map$tb_act_snps$act_snp[ss]
          act_resp_ss_mm <- unlist(list_map$map_act_snps_responses[ss, mm])

          if (sum(pat[act_snp_ss, act_resp_ss_mm]) == 0) {

            pat[act_snp_ss, sample(act_resp_ss_mm, 1)] <- 1

          }
        }

      }

    }

    pat

  }

  pat <- add_missing_activations(pat, list_map$list_map_annots, n_modules)
  pat <- add_missing_activations(pat, list_map$list_map_indep, n_modules)

  vec_rho <- runif(n_modules, min = rho_min, max = rho_max)

  if (n_modules > 1) {

    # for the case the modules do not have the same number of responses
    # (not handled by generate_phenos, which would split by blocks of equal size)
    #
    list_phenos <- lapply(1:n_modules, function(mm) {
      q_m <- sum(list_modules$vec_fac_bl == mm)
      generate_phenos(n, q_m, cor_type = "equicorrelated",
                      vec_rho = vec_rho[mm],
                      n_cpus = n_cpus, var_err = sd_err^2)$phenos
    })

    phenos <- do.call(cbind, list_phenos)

  } else {

    phenos <- generate_phenos(n, q, cor_type = "equicorrelated",
                              vec_rho = vec_rho,
                              n_cpus = n_cpus, var_err = sd_err^2)$phenos

  }

  ind_p0 <- which(rowSums(pat)>0)
  ind_q0 <- which(colSums(pat)>0)

  if(is.null(max_tot_pve)) {

    beta[pat == 1] <- rnorm(sum(pat == 1), sd = sd_act_beta)

  } else {

    vec_maf <- apply(X, 2, mean) / 2
    var_err <- apply(phenos, 2, var) # empirical error variance

    max_per_resp <- max(colSums(pat))

    pve_per_snp <- max_tot_pve / max_per_resp

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

    beta[, ind_q0] <- sapply(ind_q0, function(k) {

      p0_k <- sum(pat[,k])
      vec_pve_per_snp <- rbeta(p0_k, shape1 = 1, shape2 = 5) # left-skewed Beta distribution,
      # to give more weight to smaller effect sizes
      vec_pve_per_snp <- vec_pve_per_snp / sum(vec_pve_per_snp) * pve_per_snp * p0_k

      tot_var_expl_k <- pve_per_snp * p0_k * var_err[k] / (1 - pve_per_snp * p0_k)

      vec_maf_act <- vec_maf[pat[,k] == 1]
      vec_var_act <- 2 * vec_maf_act * (1 - vec_maf_act)

      beta_k <- rep(0, p)
      beta_k[pat[,k] == 1] <- sqrt((tot_var_expl_k + var_err[k]) * vec_pve_per_snp / vec_var_act)

      # switches signs with probabilty 0.5
      beta_k[pat[,k] == 1] <- sample(c(1, -1), p0_k, replace = TRUE) * beta_k[pat[,k] == 1]

      beta_k
    })

  }

  # Y <- phenos + scale(X) %*% beta
  Y <- phenos + X %*% beta #######################

  stopifnot(identical(pat>0, abs(beta)>0))

  create_named_list_(Y,
                     pat,
                     beta,
                     ind_p0, # all the variables involved in associations, can be larger than p0
                     ind_q0)


}

generate_binary_annots_ <- function(p,
                                   r,
                                   list_map_annots = NULL # object obtained from the set_locus_and_module_pattern_from_annots_ function, if NULL, annots drawn from noise
) {

  # Each SNP has a given probability of falling in annots
  #
  prob_annots <- rbeta(p, shape1 = 1, shape2 = 25) # left-skewed (a few SNPs will have the tendency to fall in many annots)

  # Simulate a matrix of binary annots
  #
  V <- t(sapply(prob_annots, function(pm) sample(c(0, 1), replace = TRUE, size = r, prob = c(1-pm, pm))))

  # Make sure that each annot concerns at least two SNPs
  # (otherwise not enough information to be acquired from the annot)
  #
  V <- apply(V, 2, function(V_l) {

    if (sum(V_l) < 2) {
      V_l[V_l == 0][sample(1:(p-sum(V_l)), size = 2 - sum(V_l))] <- 1
    }
    V_l
  })

  if (!is.null(list_map_annots)) {

    # Add a one for each active entry of SNP / annot
    #
    map_act_snps_modules <- list_map_annots$map_act_snps_modules
    tb_act_annots <- list_map_annots$tb_act_annots
    tb_act_snps <- list_map_annots$tb_act_snps

    n_modules <- ncol(map_act_snps_modules)

    for (mm in 1:n_modules) {

      act_annots_mm <- tb_act_annots$act_annot[tb_act_annots$module == mm]

      act_snps_mm <- tb_act_snps$act_snp[map_act_snps_modules[, mm]]

      V[act_snps_mm, act_annots_mm] <- 1
    }

  }

  colnames(V) <- paste0("Annot_z_", 1:r)
  V

}


generate_snps_from_loci_ <- function(n,
                                    list_loci,
                                    real_snp_mat = NULL,
                                    rho_min = NULL,
                                    rho_max = NULL,
                                    maf_thres = 0.05,
                                    min_dist = 0) {

  check_natural_(n)
  check_natural_(min_dist, zero_ok = TRUE)
  check_zero_one_(maf_thres)

  p <- list_loci$n_var_blocks
  n_loci <- list_loci$n_bl

  if (is.null(real_snp_mat)) {

    stopifnot(min_dist == 0) # simulated SNPs so no minimal distance to be specified

    check_zero_one_(rho_min)
    check_zero_one_(rho_max)

    vec_rho <- runif(n_loci, min = rho_min, max = rho_max)

    # echoseq::generate_snps(n = n, p = p, vec_rho = vec_rho, cor_type = "autocorrelated")$snps

    # for the case the loci do not have the same number of SNPs
    # (not handled by generate_snps, which would split by blocks of equal size)
    #
    list_snps <- lapply(1:n_loci, function(ll) {
      p_l <- sum(list_loci$vec_fac_bl == ll)
      generate_snps(n = n, p = p_l, vec_rho = vec_rho[ll],
                    cor_type = "autocorrelated",
                    vec_maf = runif(p_l, min = maf_thres, max = 0.5))$snps
    })

    X <- do.call(cbind, list_snps)

  } else {

    stopifnot(is.null(rho_min) & is.null(rho_max))

    real_n <- nrow(real_snp_mat)
    stopifnot(real_n >= n)
    stopifnot(all(as.vector(real_snp_mat) %in% c(0:2)))


    # Sample loci from the real snp matrix
    #
    real_snp_mat <- real_snp_mat[sample(1:real_n, n), , drop = FALSE]

    # Exclude rare variants after restricting to the considered samples
    #
    maf <- apply(real_snp_mat, 2, function(s) mean(s) / 2)
    real_snp_mat <-  real_snp_mat[, maf > maf_thres, drop = FALSE]

    real_p <- ncol(real_snp_mat)
    stopifnot(real_p > p + min_dist * n_loci) # in order to have sufficient SNPs to draw loci from (randomly)

    vec_ch <- make_chunks_(1:real_p, n_loci)
    len_loci <- table(list_loci$vec_fac_bl)

    stopifnot(all(sapply(vec_ch, length) >= len_loci)) # each loci must fit within its corresponding chunk
    stopifnot(min_dist %in% 0:floor((real_p - sum(len_loci)) / n_loci))

    vec_ind <- unlist(lapply(1:n_loci, function(ll) {
      ch <- vec_ch[[ll]]
      len <- len_loci[ll]

      start_choices <- ch[1:(length(ch) - len - min_dist)]

      if (length(start_choices) > 1) {
        start <- sample(start_choices, size = 1)
      } else {
        start <- start_choices
      }

      start:(start + len - 1)

    }))

    X <- real_snp_mat[, vec_ind, drop = FALSE]

  }

  rownames(X) <- paste0("ind_", 1:nrow(X)) # in order to have the same sample
  # names as in the simulated
  # response matrix

  X

}

rm_bin_annot_freq_ <- function(mat, bin_annot_freq) {

  if (!is.null(bin_annot_freq)) {

    bool_bin_annot_freq <- apply(mat, 2, function(x) {
      if (check_binary_(x) &
          (sum(x) < bin_annot_freq * length(x) | sum(1-x) < bin_annot_freq * length(x))) {
        TRUE
      } else {
        FALSE
      }
    })

    mat[, !bool_bin_annot_freq, drop = FALSE]

  } else {

    mat

  }

}


