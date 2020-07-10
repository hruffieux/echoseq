# This file is part of the `echoseq` R package:
#     https://github.com/hruffieux/echoseq
#

# Diverse utility functions implementing sanity checks and basic operations.
#
check_natural_ <- function(x, zero_ok = FALSE, eps = .Machine$double.eps^0.75){

  stopifnot(!is.null(x))

  fac <- ifelse(zero_ok, -1, 1)

  if (any(x < fac * eps | abs(x - round(x)) > eps)) {
    stop(paste0(deparse(substitute(x)),
                " must be natural."))
  }
}


check_binary_ <-function(x) {
  identical(as.vector(x), as.numeric(as.logical(x)))
}


check_positive_ <- function(x, eps = .Machine$double.eps^0.75){
  if (any(x < eps)) {
    err_mess <- paste(deparse(substitute(x)), " must be positive.", sep="")
    if (length(x) > 1) err_mess <- paste("All entries of ", err_mess, sep="")
    stop(err_mess)
  }
}

check_zero_one_ <- function(x){
  if (any(x < 0) | any(x > 1)) {
    err_mess <- paste(deparse(substitute(x)), " must lie between 0 and 1.", sep="")
    if (length(x) > 1) err_mess <- paste("All entries of ", err_mess, sep="")
    stop(err_mess)
  }
}

check_structure_ <- function(x, struct, type, size = NULL,
                             null_ok = FALSE,  inf_ok = FALSE, na_ok = FALSE) {
  if (type == "double") {
    bool_type <-  is.double(x)
    type_mess <- "a double-precision "
  } else if (type == "integer") {
    bool_type <- is.integer(x)
    type_mess <- "an integer "
  } else if (type == "numeric") {
    bool_type <- is.numeric(x)
    type_mess <- "a numeric "
  } else if (type == "logical") {
    bool_type <- is.logical(x)
    type_mess <- "a boolean "
  } else if (type == "string") {
    bool_type <- is.character(x)
    type_mess <- "string "
  }

  bool_size <- TRUE # for case size = NULL (no assertion on the size/dimension)
  size_mess <- ""
  if (struct == "vector") {
    bool_struct <- is.vector(x) & (length(x) > 0) # not an empty vector
    if (!is.null(size)) {
      bool_size <- length(x) %in% size
      size_mess <- paste(" of length ", paste(size, collapse=" or "), sep = "")
    }
  } else if (struct == "matrix") {
    bool_struct <- is.matrix(x) & (length(x) > 0) # not an empty matrix
    if (!is.null(size)) {
      bool_size <- all(dim(x) == size)
      size_mess <- paste(" of dimension ", size[1], " x ", size[2], sep = "")
    }
  }

  correct_obj <- bool_struct & bool_type & bool_size

  bool_null <- is.null(x)

  if (!is.list(x) & type != "string") {
    na_mess <- ""
    if (!na_ok) {
      if (!bool_null) correct_obj <- correct_obj & !any(is.na(x))
      na_mess <- " without missing value"
    }

    inf_mess <- ""
    if (!inf_ok) {
      if (!bool_null) correct_obj <- correct_obj & all(is.finite(x[!is.na(x)]))
      inf_mess <- ", finite"
    }
  } else {
    na_mess <- ""
    inf_mess <- ""
  }

  null_mess <- ""
  if (null_ok) {
    correct_obj <- correct_obj | bool_null
    null_mess <- " or must be NULL"
  }

  if(!(correct_obj)) {
    stop(paste(deparse(substitute(x)), " must be a non-empty ", type_mess, struct,
               size_mess, inf_mess, na_mess, null_mess, ".", sep = ""))
  }
}


create_named_list_ <- function(...) {
  setNames(list(...), as.character(match.call()[-1]))
}


make_chunks_ <- function(x, n_g) split(x, factor(sort(rank(x) %% n_g)))



cbind_fill_matrix <- function(...) {
  tr <- lapply(..., as.matrix)
  tr <- lapply(..., t)
  t(as.matrix(plyr::rbind.fill.matrix(tr)))
}


set_blocks_ <- function(tot, pos_bl, n_cpus, verbose = TRUE) {

  check_structure_(n_cpus, "vector", "numeric", 1)
  check_natural_(n_cpus)

  check_structure_(verbose, "vector", "logical", 1)

  list_blocks <- lapply(seq_along(tot), function(ii) {

    tt <- tot[ii]

    if (is.list(pos_bl)) {
      pb <- pos_bl[[ii]]
    } else {
      pb <- pos_bl
    }

    check_structure_(tt, "vector", "numeric", 1)
    check_natural_(tt)

    check_structure_(pb, "vector", "numeric")
    check_natural_(pb)

    if (any(pb < 1) | any(pb > tt))
      stop("The positions provided in pos_bl must range between 1 and total number of variables given in tot.")

    if (any(duplicated(pb)))
      stop("The positions provided in pos_bl must be unique.")

    if (any(pb != cummax(pb)))
      stop("The positions provided in pos_bl must be monotonically increasing.")

    vec_fac_bl <- as.factor(cumsum(seq_along(1:tt) %in% pb))

    n_bl <- length(unique(vec_fac_bl))

    n_var_blocks <- tt

    create_named_list_(n_var_blocks, n_bl, vec_fac_bl)

  })

  if (length(list_blocks) > 1 && list_blocks[[2]]$n_bl != 1) {
    names(list_blocks) <- c("bl_x", "bl_y")
    tot_n_bl <- list_blocks$bl_x$n_bl * list_blocks$bl_y$n_bl
  } else {
    list_blocks <- list_blocks[[1]]
    tot_n_bl <- list_blocks$n_bl
  }

  if (n_cpus > 1) {

    n_cpus_avail <- parallel::detectCores()
    if (n_cpus > n_cpus_avail) {
      n_cpus <- n_cpus_avail
      warning(paste("The number of CPUs specified exceeds the number of CPUs ",
                    "available on the machine. The latter has been used instead.",
                    sep=""))
    }

    if (n_cpus > tot_n_bl){
      message <- paste("The number of cpus in use is at most equal to the number of blocks.",
                       "n_cpus is therefore set to ", tot_n_bl, ". \n", sep ="")
      if(verbose) cat(message)
      else warning(message)
      n_cpus <- tot_n_bl
    }


  }

  list_blocks$n_cpus <- n_cpus

  class(list_blocks) <- "blocks"

  list_blocks
}
