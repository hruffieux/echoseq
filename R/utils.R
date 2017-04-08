# This file is part of the `mimeo` R package:
#     https://github.com/hruffieux/mimeo
#

# Diverse utility functions implementing sanity checks and basic operations.
#

check_natural_ <- function(x, eps = .Machine$double.eps^0.75){
  if (any(x < eps | abs(x - round(x)) > eps)) {
    stop(paste(deparse(substitute(x)),
               " must be natural.", sep=""))
  }
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
