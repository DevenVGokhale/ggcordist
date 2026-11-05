#' Convert Stan output (or arrays) to a `cov_draws` object
#'
#' `as_cov_draws()` extracts draws of a correlation or covariance matrix
#' parameter from a fitted model and returns a [cov_draws][new_cov_draws]
#' object that can be passed to the plotting functions.
#'
#' @param x One of:
#'   * a `stanfit` object from \pkg{rstan},
#'   * a `CmdStanMCMC` object from \pkg{cmdstanr},
#'   * a `brmsfit` object from \pkg{brms},
#'   * a `draws` object from \pkg{posterior},
#'   * a numeric `S x K x K` array of matrices,
#'   * a `cov_draws` object (returned unchanged).
#' @param par Name of the matrix parameter in the Stan program, e.g.
#'   `"Omega"` or `"Sigma"`. Not used for `brmsfit` or array input.
#' @param type What `par` contains: a correlation matrix (`"cor"`), a
#'   covariance matrix (`"cov"`), the lower Cholesky factor of a correlation
#'   matrix (`"chol_cor"`) or of a covariance matrix (`"chol_cov"`).
#' @param sd_par Optional name of a vector parameter holding the standard
#'   deviations, used when `type` is `"cor"` or `"chol_cor"`. It is ignored
#'   for covariance types because the standard deviations are implied.
#' @param group For `brmsfit` objects: the grouping factor whose
#'   random-effect correlations (`cor_<group>__a__b`) and standard deviations
#'   (`sd_<group>__a`) should be used. Use `NULL` together with
#'   `residual = TRUE` for the residual correlations of a multivariate model.
#' @param residual For `brmsfit` objects: use the residual correlations
#'   (`rescor__a__b`) and residual standard deviations (`sigma_a`) instead.
#' @param sd For array input and `type` `"cor"`/`"chol_cor"`: optional `S x K`
#'   matrix of standard deviations.
#' @param labels Optional variable labels (length `K`).
#' @param source Label describing the draws, e.g. `"posterior"` or `"prior"`.
#' @param ... Passed to methods.
#'
#' @return A [cov_draws][new_cov_draws] object.
#' @export
#' @examples
#' # An S x K x K array of covariance matrices
#' arr <- cov_array(r_inv_wishart(50, k = 3))
#' as_cov_draws(arr, type = "cov")
as_cov_draws <- function(x, ...) {
  UseMethod("as_cov_draws")
}

#' @rdname as_cov_draws
#' @export
as_cov_draws.cov_draws <- function(x, ...) {
  x
}

#' @rdname as_cov_draws
#' @export
as_cov_draws.default <- function(x, ...) {
  cli::cli_abort(
    "Don't know how to convert {.obj_type_friendly {x}} to {.cls cov_draws}."
  )
}

#' @rdname as_cov_draws
#' @export
as_cov_draws.array <- function(
  x,
  type = c("cor", "cov", "chol_cor", "chol_cov"),
  sd = NULL,
  labels = NULL,
  source = "draws",
  ...
) {
  type <- rlang::arg_match(type)
  if (length(dim(x)) != 3L) {
    cli::cli_abort("Array input must have dimensions {.code S x K x K}.")
  }
  build_cov_draws(x, type = type, sd = sd, labels = labels, source = source)
}

#' @rdname as_cov_draws
#' @export
as_cov_draws.draws <- function(
  x,
  par,
  type = c("cor", "cov", "chol_cor", "chol_cov"),
  sd_par = NULL,
  labels = NULL,
  source = "posterior",
  ...
) {
  type <- rlang::arg_match(type)
  if (missing(par)) {
    cli::cli_abort("{.arg par} is required for {.cls draws} objects.")
  }
  rv <- posterior::as_draws_rvars(x)
  if (!par %in% names(rv)) {
    cli::cli_abort(c(
      "Parameter {.val {par}} not found in the draws.",
      i = "Available variables: {.val {names(rv)}}."
    ))
  }
  arr <- posterior::draws_of(rv[[par]])
  if (length(dim(arr)) != 3L) {
    cli::cli_abort(
      "{.val {par}} must be a matrix parameter, but has dimension {.val {dim(arr)[-1]}}."
    )
  }
  sd <- NULL
  if (!is.null(sd_par)) {
    if (!sd_par %in% names(rv)) {
      cli::cli_abort("Parameter {.val {sd_par}} not found in the draws.")
    }
    sd <- posterior::draws_of(rv[[sd_par]])
  }
  build_cov_draws(arr, type = type, sd = sd, labels = labels, source = source)
}

#' @rdname as_cov_draws
#' @export
as_cov_draws.stanfit <- function(
  x,
  par,
  type = c("cor", "cov", "chol_cor", "chol_cov"),
  sd_par = NULL,
  labels = NULL,
  source = "posterior",
  ...
) {
  rlang::check_installed("rstan")
  type <- rlang::arg_match(type)
  if (missing(par)) {
    cli::cli_abort("{.arg par} is required for {.cls stanfit} objects.")
  }
  pars <- c(par, sd_par)
  ex <- rstan::extract(x, pars = pars, permuted = TRUE)
  arr <- ex[[par]]
  sd <- if (!is.null(sd_par)) ex[[sd_par]] else NULL
  build_cov_draws(arr, type = type, sd = sd, labels = labels, source = source)
}

#' @rdname as_cov_draws
#' @export
as_cov_draws.CmdStanMCMC <- function(
  x,
  par,
  type = c("cor", "cov", "chol_cor", "chol_cov"),
  sd_par = NULL,
  labels = NULL,
  source = "posterior",
  ...
) {
  type <- rlang::arg_match(type)
  if (missing(par)) {
    cli::cli_abort("{.arg par} is required for {.cls CmdStanMCMC} objects.")
  }
  d <- x$draws(variables = c(par, sd_par), format = "draws_array")
  as_cov_draws(
    d,
    par = par,
    type = type,
    sd_par = sd_par,
    labels = labels,
    source = source
  )
}

#' @rdname as_cov_draws
#' @export
as_cov_draws.brmsfit <- function(
  x,
  group = NULL,
  residual = FALSE,
  labels = NULL,
  source = "posterior",
  ...
) {
  if (isTRUE(residual)) {
    cor_rx <- "^rescor__"
    sd_rx <- "^sigma_"
    var_regex <- "^(rescor__|sigma_)"
  } else if (!is.null(group)) {
    cor_rx <- paste0("^cor_", group, "__")
    sd_rx <- paste0("^sd_", group, "__")
    var_regex <- paste0("^(cor|sd)_", group, "__")
  } else {
    cli::cli_abort("Supply {.arg group} or set {.code residual = TRUE}.")
  }
  m <- posterior::as_draws_matrix(x, variable = var_regex, regex = TRUE)
  brms_matrix_to_cov_draws(m, cor_rx, sd_rx, labels, source)
}

# Parse brms-style `cor_<g>__a__b` / `sd_<g>__a` columns of a draws matrix
brms_matrix_to_cov_draws <- function(m, cor_rx, sd_rx, labels, source) {
  nms <- colnames(m)
  sd_nms <- grep(sd_rx, nms, value = TRUE)
  cor_nms <- grep(cor_rx, nms, value = TRUE)
  if (length(sd_nms) < 2L) {
    cli::cli_abort(
      "Fewer than two standard deviation parameters found; need at least a 2 x 2 matrix."
    )
  }
  terms <- sub(sd_rx, "", sd_nms)
  K <- length(terms)
  S <- nrow(m)

  cor_arr <- array(0, dim = c(S, K, K))
  for (k in seq_len(K)) {
    cor_arr[, k, k] <- 1
  }
  for (nm in cor_nms) {
    pair <- strsplit(sub(cor_rx, "", nm), "__", fixed = TRUE)[[1]]
    if (length(pair) != 2L || !all(pair %in% terms)) {
      cli::cli_abort("Could not parse correlation parameter {.val {nm}}.")
    }
    i <- match(pair[1], terms)
    j <- match(pair[2], terms)
    cor_arr[, i, j] <- as.numeric(m[, nm])
    cor_arr[, j, i] <- as.numeric(m[, nm])
  }
  sd <- matrix(as.numeric(m[, sd_nms]), nrow = S)
  new_cov_draws(cor_arr, sd = sd, labels = labels %||% terms, source = source)
}

# Internal helpers ---------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

build_cov_draws <- function(arr, type, sd, labels, source) {
  S <- dim(arr)[1]
  K <- dim(arr)[2]
  if (dim(arr)[3] != K) {
    cli::cli_abort("Matrix draws must be square.")
  }

  if (type %in% c("chol_cor", "chol_cov")) {
    full <- array(0, dim = dim(arr))
    for (s in seq_len(S)) {
      L <- matrix(arr[s, , ], K, K)
      L[upper.tri(L)] <- 0
      full[s, , ] <- tcrossprod(L)
    }
    arr <- full
    type <- sub("chol_", "", type, fixed = TRUE)
  }

  if (type == "cov") {
    sd <- matrix(NA_real_, S, K)
    for (k in seq_len(K)) {
      sd[, k] <- sqrt(arr[, k, k])
    }
    cor <- arr
    for (i in seq_len(K)) {
      for (j in seq_len(K)) {
        cor[, i, j] <- arr[, i, j] / (sd[, i] * sd[, j])
      }
    }
  } else {
    cor <- arr
  }

  new_cov_draws(cor, sd = sd, labels = labels, source = source)
}
