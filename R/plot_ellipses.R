#' Equiprobability ellipses of a 2 x 2 submatrix
#'
#' For each (thinned) draw, plots the equiprobability ellipse of a zero-mean
#' bivariate normal whose covariance is the `(i, j)` submatrix of the draw
#' (layer 3 of Tokuda et al.). Ellipses are coloured by effective dependence
#' as in [plot_tableau()]. If `x` has no standard deviations, the correlation
#' submatrix is used, i.e. unit variances.
#'
#' @inheritParams plot_tableau
#' @param i,j Indices of the two variables.
#' @param n Number of ellipses to draw (draws are thinned evenly).
#' @return A ggplot object.
#' @export
#' @examples
#' plot_ellipses(r_inv_wishart(500, k = 4))
plot_ellipses <- function(
  x,
  i = 1,
  j = 2,
  n = 100,
  level = 0.5,
  extreme = 0.05,
  ...
) {
  x <- as_cov_draws(x, ...)
  K <- n_vars(x)
  if (i == j || !all(c(i, j) %in% seq_len(K))) {
    cli::cli_abort(
      "{.arg i} and {.arg j} must be two different indices between 1 and {K}."
    )
  }
  if (!is.numeric(level) || level <= 0 || level >= 1) {
    cli::cli_abort("{.arg level} must be strictly between 0 and 1.")
  }
  group <- extreme_group(effective_dependence(x), extreme)
  draws <- thin_index(n_draws(x), n)

  r <- sqrt(stats::qchisq(level, df = 2))
  theta <- seq(0, 2 * pi, length.out = 101)
  circle <- rbind(cos(theta), sin(theta))

  pieces <- lapply(draws, function(s) {
    si <- if (is.null(x$sd)) 1 else x$sd[s, i]
    sj <- if (is.null(x$sd)) 1 else x$sd[s, j]
    rho <- x$cor[s, i, j]
    S2 <- matrix(c(si^2, si * sj * rho, si * sj * rho, sj^2), 2L, 2L)
    e <- eigen(S2, symmetric = TRUE)
    A <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), 2L)
    p <- r * (A %*% circle)
    data.frame(draw = s, x = p[1, ], y = p[2, ], group = group[s])
  })
  d <- do.call(rbind, pieces)
  d <- d[order(d$group), , drop = FALSE]
  d$draw <- factor(d$draw, levels = unique(d$draw))

  ggplot2::ggplot(
    d,
    ggplot2::aes(.data$x, .data$y, group = .data$draw, colour = .data$group)
  ) +
    ggplot2::geom_path(alpha = 0.5, linewidth = 0.3) +
    ggplot2::scale_colour_manual(values = cov_palette, guide = "none") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = sprintf("X[%s]", x$labels[i]),
      y = sprintf("X[%s]", x$labels[j]),
      subtitle = sprintf(
        "%d %g%% equiprobability ellipses",
        length(draws),
        100 * level
      )
    ) +
    ggplot2::theme_bw()
}

#' Three-dimensional scatterplot of three correlations
#'
#' Plots `(rho_ij, rho_ik, rho_jk)` for three variables in an interactive
#' \pkg{rgl} window (layer 3 of Tokuda et al.). The support of three
#' correlations is an elliptical tetrahedron. Draws are coloured by
#' effective dependence as in [plot_tableau()].
#'
#' @inheritParams plot_tableau
#' @param idx Three distinct variable indices.
#' @return Invisibly, the result of `rgl::plot3d()`. Requires the \pkg{rgl}
#'   package.
#' @export
plot_cor3d <- function(x, idx = 1:3, extreme = 0.05, ...) {
  rlang::check_installed("rgl", reason = "to draw 3D scatterplots.")
  x <- as_cov_draws(x, ...)
  if (length(idx) != 3L) {
    cli::cli_abort("{.arg idx} must have length 3.")
  }
  df <- cov_frame(x, idx, extreme)
  lab <- var_labels(x, idx)
  cols <- unname(cov_palette[as.character(df$group)])
  invisible(rgl::plot3d(
    df$r12,
    df$r13,
    df$r23,
    col = cols,
    size = 3,
    type = "p",
    xlab = lab$r12,
    ylab = lab$r13,
    zlab = lab$r23
  ))
}
