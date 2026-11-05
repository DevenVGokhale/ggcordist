#' Reference samplers for covariance matrix distributions
#'
#' Draw from the distributions discussed by Tokuda et al., returned as
#' [cov_draws][new_cov_draws] objects so they can be plotted directly or
#' compared against posterior draws with [plot_compare()].
#'
#' Parameterisations follow the paper:
#'
#' * `r_wishart()`: `Sigma ~ Wishart(df, scale)`.
#' * `r_inv_wishart()`: `Sigma ~ inv-Wishart(df, scale^-1)`, with density
#'   proportional to `|Sigma|^-(df+k+1)/2 exp(-tr(scale Sigma^-1)/2)` and mean
#'   `scale / (df - k - 1)`.
#' * `r_scaled_inv_wishart()`: `Sigma = diag(xi) Q diag(xi)` with
#'   `Q ~ inv-Wishart(df, I)`; `xi` standard half-normal by default.
#' * `r_lkj()`: correlation matrices from the LKJ distribution with shape
#'   `eta` (uniform on correlation matrices when `eta = 1`), drawn with the
#'   C-vine method of Lewandowski, Kurowicka and Joe (2009). Standard
#'   deviations are not included.
#' * `r_separation()`: separation strategy of Barnard, McCulloch and Meng
#'   (2000), `Sigma = diag(sd) R diag(sd)`, with `R` either LKJ or the
#'   correlation matrix of an inverse-Wishart draw, and standard half-normal
#'   standard deviations by default.
#'
#' @param n Number of draws.
#' @param k Dimension of the matrices.
#' @param df Degrees of freedom. For [r_inv_wishart()] and
#'   [r_scaled_inv_wishart()] the default is `k + 1`, which makes every
#'   marginal correlation uniform on (-1, 1).
#' @param scale Positive definite `k x k` scale matrix (default identity).
#' @param eta LKJ shape parameter, `eta > 0`.
#' @param sd_sampler A function of one argument `n` returning `n` positive
#'   standard deviations (or scale parameters). Defaults to a standard
#'   half-normal.
#' @param cor_dist For `r_separation()`: distribution of the correlation
#'   matrix, `"lkj"` or `"inv_wishart"`.
#' @param labels Optional variable labels.
#'
#' @return A [cov_draws][new_cov_draws] object (with `sd` unless noted).
#' @name samplers
#' @examples
#' set.seed(1)
#' r_inv_wishart(100, k = 4)
#' r_lkj(100, k = 4, eta = 1)
NULL

half_normal <- function(n) abs(stats::rnorm(n))

#' @rdname samplers
#' @export
r_wishart <- function(n, k, df = k, scale = diag(k), labels = NULL) {
  check_sampler_args(n, k, df, scale)
  if (df < k) {
    cli::cli_abort(
      "{.arg df} must be at least {.arg k} for a non-singular Wishart."
    )
  }
  w <- stats::rWishart(n, df = df, Sigma = scale)
  cov_from_array(aperm(w, c(3L, 1L, 2L)), labels, "Wishart")
}

#' @rdname samplers
#' @export
r_inv_wishart <- function(n, k, df = k + 1, scale = diag(k), labels = NULL) {
  check_sampler_args(n, k, df, scale)
  if (df < k) {
    cli::cli_abort(
      "{.arg df} must be at least {.arg k} for a proper inverse-Wishart."
    )
  }
  cov_from_array(inv_wishart_array(n, df, scale), labels, "inverse-Wishart")
}

#' @rdname samplers
#' @export
r_scaled_inv_wishart <- function(
  n,
  k,
  df = k + 1,
  scale = diag(k),
  sd_sampler = half_normal,
  labels = NULL
) {
  check_sampler_args(n, k, df, scale)
  q <- inv_wishart_array(n, df, scale)
  xi <- vapply(seq_len(k), function(i) sd_sampler(n), numeric(n))
  xi <- matrix(xi, nrow = n, ncol = k)
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      q[, i, j] <- q[, i, j] * xi[, i] * xi[, j]
    }
  }
  cov_from_array(q, labels, "scaled inverse-Wishart")
}

#' @rdname samplers
#' @export
r_lkj <- function(n, k, eta = 1, labels = NULL) {
  check_sampler_args(n, k, df = k, scale = diag(k))
  if (!is.numeric(eta) || length(eta) != 1L || eta <= 0) {
    cli::cli_abort("{.arg eta} must be a single positive number.")
  }
  cor <- lkj_array(n, k, eta)
  new_cov_draws(
    cor,
    sd = NULL,
    labels = labels,
    source = paste0("LKJ(", eta, ")")
  )
}

#' @rdname samplers
#' @export
r_separation <- function(
  n,
  k,
  cor_dist = c("lkj", "inv_wishart"),
  eta = 1,
  df = k + 1,
  sd_sampler = half_normal,
  labels = NULL
) {
  cor_dist <- rlang::arg_match(cor_dist)
  check_sampler_args(n, k, df, diag(k))
  cor <- if (cor_dist == "lkj") {
    lkj_array(n, k, eta)
  } else {
    arr <- inv_wishart_array(n, df, diag(k))
    cov_from_array(arr, NULL, "")$cor
  }
  sd <- matrix(
    vapply(seq_len(k), function(i) sd_sampler(n), numeric(n)),
    nrow = n,
    ncol = k
  )
  src <- if (cor_dist == "lkj") {
    paste0("separation (LKJ(", eta, "))")
  } else {
    "separation (inverse-Wishart)"
  }
  new_cov_draws(cor, sd = sd, labels = labels, source = src)
}

# Internal -----------------------------------------------------------------

check_sampler_args <- function(n, k, df, scale) {
  if (!is.numeric(n) || length(n) != 1L || n < 1) {
    cli::cli_abort("{.arg n} must be a positive integer.")
  }
  if (!is.numeric(k) || length(k) != 1L || k < 2) {
    cli::cli_abort("{.arg k} must be an integer >= 2.")
  }
  if (!is.matrix(scale) || !all(dim(scale) == k)) {
    cli::cli_abort("{.arg scale} must be a {k} x {k} matrix.")
  }
  invisible(TRUE)
}

cov_from_array <- function(arr, labels, source) {
  S <- dim(arr)[1]
  K <- dim(arr)[2]
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
  new_cov_draws(cor, sd = sd, labels = labels, source = source)
}

# Sigma^-1 ~ Wishart(df, solve(scale))  =>  Sigma ~ inv-Wishart(df, scale^-1)
inv_wishart_array <- function(n, df, scale) {
  k <- nrow(scale)
  w <- stats::rWishart(n, df = df, Sigma = solve(scale))
  out <- array(NA_real_, dim = c(n, k, k))
  for (s in seq_len(n)) {
    out[s, , ] <- chol2inv(chol(w[,, s]))
  }
  out
}

# C-vine method, Lewandowski, Kurowicka and Joe (2009), Section 2.4
lkj_array <- function(n, k, eta) {
  out <- array(NA_real_, dim = c(n, k, k))
  for (s in seq_len(n)) {
    out[s, , ] <- lkj_one(k, eta)
  }
  out
}

lkj_one <- function(k, eta) {
  beta <- eta + (k - 1) / 2
  P <- matrix(0, k, k)
  R <- diag(k)
  for (i in seq_len(k - 1L)) {
    beta <- beta - 1 / 2
    for (j in (i + 1L):k) {
      P[i, j] <- 2 * stats::rbeta(1, beta, beta) - 1
      p <- P[i, j]
      if (i > 1L) {
        for (l in (i - 1L):1L) {
          p <- p * sqrt((1 - P[l, i]^2) * (1 - P[l, j]^2)) + P[l, i] * P[l, j]
        }
      }
      R[i, j] <- p
      R[j, i] <- p
    }
  }
  R
}
