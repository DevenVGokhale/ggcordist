#' Draws of covariance / correlation matrices
#'
#' `cov_draws` is the common container used by all plotting functions in
#' ggcordist. It stores `S` draws of a `K x K` correlation matrix and,
#' optionally, the `K` standard deviations for each draw (so that the full
#' covariance matrix `Sigma = diag(sd) R diag(sd)` can be reconstructed).
#'
#' Use [as_cov_draws()] to create one from a Stan fit, and the `r_*()`
#' samplers (e.g. [r_lkj()]) to create reference draws.
#'
#' @param cor A numeric `S x K x K` array of correlation matrices.
#' @param sd Optional numeric `S x K` matrix of standard deviations.
#' @param labels Optional character vector of `K` variable labels.
#' @param source A single string describing where the draws come from, for
#'   example `"prior"` or `"posterior"`. Used as the legend label by
#'   [plot_compare()].
#' @param tol Numerical tolerance for the symmetry and unit-diagonal checks.
#'
#' @return An object of class `cov_draws`: a list with elements `cor`, `sd`,
#'   `labels` and `source`.
#' @export
#' @examples
#' x <- r_lkj(200, k = 3, eta = 2)
#' x
new_cov_draws <- function(
  cor,
  sd = NULL,
  labels = NULL,
  source = "draws",
  tol = 1e-4
) {
  if (!is.array(cor) || length(dim(cor)) != 3L || dim(cor)[2] != dim(cor)[3]) {
    cli::cli_abort("{.arg cor} must be a numeric {.code S x K x K} array.")
  }
  S <- dim(cor)[1]
  K <- dim(cor)[2]
  if (K < 2L) {
    cli::cli_abort("Matrices must be at least 2 x 2.")
  }
  if (anyNA(cor)) {
    cli::cli_abort("{.arg cor} contains missing values.")
  }

  if (!is.null(sd)) {
    sd <- as.matrix(sd)
    if (!all(dim(sd) == c(S, K))) {
      cli::cli_abort("{.arg sd} must be an {.code S x K} matrix ({S} x {K}).")
    }
    if (any(sd <= 0)) {
      cli::cli_abort("{.arg sd} must be strictly positive.")
    }
    dimnames(sd) <- NULL
  }

  if (is.null(labels)) {
    labels <- as.character(seq_len(K))
  }
  if (length(labels) != K) {
    cli::cli_abort("{.arg labels} must have length {K}.")
  }

  asym <- max(abs(cor - aperm(cor, c(1L, 3L, 2L))))
  if (asym > tol) {
    cli::cli_abort(
      "Correlation draws are not symmetric (max asymmetry {signif(asym, 3)})."
    )
  }
  dev <- max(vapply(
    seq_len(K),
    function(k) max(abs(cor[, k, k] - 1)),
    numeric(1)
  ))
  if (dev > tol) {
    cli::cli_abort(
      "Correlation draws do not have a unit diagonal (max deviation {signif(dev, 3)})."
    )
  }

  # Remove tiny numerical noise (e.g. from 6-digit CSV output)
  cor <- (cor + aperm(cor, c(1L, 3L, 2L))) / 2
  for (k in seq_len(K)) {
    cor[, k, k] <- 1
  }
  dimnames(cor) <- NULL

  structure(
    list(
      cor = cor,
      sd = sd,
      labels = as.character(labels),
      source = as.character(source)[1]
    ),
    class = "cov_draws"
  )
}

#' @export
print.cov_draws <- function(x, ...) {
  d <- dim(x$cor)
  cli::cli_text(
    "{.cls cov_draws}: {d[1]} draws of {d[2]} x {d[3]} matrices (source: {x$source})"
  )
  cli::cli_text("Variables: {paste(x$labels, collapse = ', ')}")
  cli::cli_text(
    "Standard deviations: {if (is.null(x$sd)) 'not available' else 'available'}"
  )
  invisible(x)
}

#' @export
length.cov_draws <- function(x) {
  dim(x$cor)[1]
}

#' Number of draws and dimension of a `cov_draws` object
#'
#' @param x A [cov_draws][new_cov_draws] object.
#' @return An integer.
#' @name cov_draws_dims
NULL

#' @rdname cov_draws_dims
#' @export
n_draws <- function(x) {
  check_cov_draws(x)
  dim(x$cor)[1]
}

#' @rdname cov_draws_dims
#' @export
n_vars <- function(x) {
  check_cov_draws(x)
  dim(x$cor)[2]
}

check_cov_draws <- function(
  x,
  arg = rlang::caller_arg(x),
  call = rlang::caller_env()
) {
  if (!inherits(x, "cov_draws")) {
    cli::cli_abort(
      "{.arg {arg}} must be a {.cls cov_draws} object, not {.obj_type_friendly {x}}. Use {.fn as_cov_draws}.",
      call = call
    )
  }
  invisible(x)
}

#' Reconstruct covariance matrices from a `cov_draws` object
#'
#' @param x A [cov_draws][new_cov_draws] object with standard deviations.
#' @return A numeric `S x K x K` array of covariance matrices.
#' @export
#' @examples
#' cov_array(r_inv_wishart(10, k = 3))
cov_array <- function(x) {
  check_cov_draws(x)
  if (is.null(x$sd)) {
    cli::cli_abort(
      "{.arg x} has no standard deviations, so covariances cannot be reconstructed."
    )
  }
  K <- n_vars(x)
  out <- x$cor
  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      out[, i, j] <- x$cor[, i, j] * x$sd[, i] * x$sd[, j]
    }
  }
  out
}
