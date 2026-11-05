#' Effective variance and effective dependence
#'
#' Scalar summaries of an entire covariance matrix (Pena and Rodriguez, 2003;
#' layer 4 of Tokuda et al.). For a `K x K` covariance matrix `Sigma` with
#' correlation matrix `R`,
#'
#' * effective variance: `|Sigma|^(1/K)`,
#' * effective dependence: `1 - |R|^(1/K)`.
#'
#' @param x A [cov_draws][new_cov_draws] object. `effective_variance()`
#'   requires standard deviations.
#' @return A numeric vector with one value per draw.
#' @export
#' @examples
#' x <- r_lkj(100, k = 4)
#' summary(effective_dependence(x))
effective_dependence <- function(x) {
  check_cov_draws(x)
  path <- effective_dependence_path(x)
  path[, ncol(path)]
}

#' @rdname effective_dependence
#' @export
effective_variance <- function(x) {
  check_cov_draws(x)
  if (is.null(x$sd)) {
    cli::cli_abort(
      "Effective variance needs standard deviations, which {.arg x} does not have."
    )
  }
  K <- n_vars(x)
  logdet_R <- log_det_path(x)[, K]
  exp((2 * rowSums(log(x$sd)) + logdet_R) / K)
}

#' Effective dependence of leading principal submatrices
#'
#' For each draw, the effective dependence `1 - |R_p|^(1/p)` of the leading
#' `p x p` submatrix of the correlation matrix, for `p = 1, ..., K`.
#'
#' @inheritParams effective_dependence
#' @return A numeric `S x K` matrix; column `p` holds the effective dependence
#'   of the leading `p x p` submatrix.
#' @export
#' @examples
#' dim(effective_dependence_path(r_lkj(20, k = 5)))
effective_dependence_path <- function(x) {
  check_cov_draws(x)
  K <- n_vars(x)
  ld <- log_det_path(x)
  1 - exp(sweep(ld, 2, seq_len(K), "/"))
}

# log |R_p| for p = 1..K via a single Cholesky factor per draw
log_det_path <- function(x) {
  S <- n_draws(x)
  K <- n_vars(x)
  out <- matrix(NA_real_, S, K)
  for (s in seq_len(S)) {
    U <- tryCatch(chol(matrix(x$cor[s, , ], K, K)), error = function(e) NULL)
    if (!is.null(U)) {
      out[s, ] <- cumsum(2 * log(diag(U)))
    }
  }
  out
}

#' Eigenvalues of the draws
#'
#' @inheritParams effective_dependence
#' @param scale Eigenvalues of the correlation matrices (`"cor"`) or of the
#'   covariance matrices (`"cov"`, requires standard deviations).
#' @return A numeric `S x K` matrix of eigenvalues in decreasing order.
#' @export
#' @examples
#' dim(cov_eigen(r_lkj(20, k = 4)))
cov_eigen <- function(x, scale = c("cor", "cov")) {
  check_cov_draws(x)
  scale <- rlang::arg_match(scale)
  arr <- if (scale == "cor") x$cor else cov_array(x)
  S <- dim(arr)[1]
  K <- dim(arr)[2]
  out <- matrix(NA_real_, S, K)
  for (s in seq_len(S)) {
    out[s, ] <- eigen(
      matrix(arr[s, , ], K, K),
      symmetric = TRUE,
      only.values = TRUE
    )$values
  }
  out
}

# Classify draws by effective dependence: lowest / highest `p` are flagged
extreme_group <- function(ed, p) {
  g <- rep("mid", length(ed))
  if (p > 0) {
    lo <- stats::quantile(ed, p, names = FALSE, na.rm = TRUE)
    hi <- stats::quantile(ed, 1 - p, names = FALSE, na.rm = TRUE)
    g[!is.na(ed) & ed <= lo] <- "low"
    g[!is.na(ed) & ed >= hi] <- "high"
  }
  factor(g, levels = c("mid", "low", "high"))
}

# Evenly thinned draw indices (deterministic, no RNG use)
thin_index <- function(S, n) {
  if (n >= S) {
    return(seq_len(S))
  }
  unique(round(seq(1, S, length.out = n)))
}

# Data frame used by the panels: variables for index set `idx`
cov_frame <- function(x, idx, extreme = 0.05) {
  K <- n_vars(x)
  idx <- as.integer(idx)
  if (anyNA(idx) || any(idx < 1L) || any(idx > K) || anyDuplicated(idx)) {
    cli::cli_abort("{.arg idx} must be distinct integers between 1 and {K}.")
  }
  df <- data.frame(draw = seq_len(n_draws(x)))
  if (!is.null(x$sd)) {
    df$lsd1 <- log(x$sd[, idx[1]])
    if (length(idx) >= 2L) df$lsd2 <- log(x$sd[, idx[2]])
  }
  if (length(idx) >= 2L) {
    df$r12 <- x$cor[, idx[1], idx[2]]
  }
  if (length(idx) >= 3L) {
    df$r13 <- x$cor[, idx[1], idx[3]]
    df$r23 <- x$cor[, idx[2], idx[3]]
  }
  if (length(idx) >= 4L) {
    df$r34 <- x$cor[, idx[3], idx[4]]
  }
  df$ed <- effective_dependence(x)
  df$group <- extreme_group(df$ed, extreme)
  df
}

default_idx <- function(x) {
  seq_len(min(n_vars(x), 4L))
}
