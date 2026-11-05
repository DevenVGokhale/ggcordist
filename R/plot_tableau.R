#' Four-layer visualization of a distribution of covariance matrices
#'
#' Reproduces the display of Tokuda et al. (Figure 1) for a
#' [cov_draws][new_cov_draws] object. Each column is one layer:
#'
#' 1. histograms of `log(sd_i)` and of a correlation `rho_ij`,
#' 2. scatterplots of pairs of log standard deviations and correlations,
#' 3. equiprobability ellipses of a zero-mean bivariate normal for the
#'    `(i, j)` submatrix (see [plot_ellipses()]),
#' 4. effective variance, effective dependence, and effective dependence of
#'    growing leading submatrices (see [plot_effective()]).
#'
#' Draws with the lowest (highest) effective dependence are coloured blue
#' (red) in every panel.
#'
#' The paper relies on exchangeability of the variables, so that one panel
#' per *type* of pair is enough. Posterior draws are generally not
#' exchangeable: choose the variables to show with `idx`, or use
#' [plot_pairs()] to see every pair.
#'
#' @param x A [cov_draws][new_cov_draws] object, or anything accepted by
#'   [as_cov_draws()].
#' @param idx Integer vector of up to four variable indices
#'   `(i, j, k, l)` used for the panels: `log(sd_i)`, `log(sd_j)`,
#'   `rho_ij`, `rho_jk` and `rho_kl`. Defaults to the first four variables.
#' @param extreme Proportion of draws coloured blue (lowest effective
#'   dependence) and red (highest) in each tail. Use `0` to switch colouring
#'   off.
#' @param n_ellipses Number of ellipses to draw.
#' @param level Probability content of each ellipse (default 50%).
#' @param bins Number of histogram bins.
#' @param ... Passed to [as_cov_draws()] when `x` is not a `cov_draws`.
#'
#' @return A \pkg{patchwork} object (a ggplot).
#' @export
#' @examples
#' x <- r_inv_wishart(1000, k = 4)
#' plot_tableau(x)
plot_tableau <- function(
  x,
  idx = NULL,
  extreme = 0.05,
  n_ellipses = 100,
  level = 0.5,
  bins = 30,
  ...
) {
  x <- as_cov_draws(x, ...)
  idx <- idx %||% default_idx(x)
  df <- cov_frame(x, idx, extreme)
  lab <- var_labels(x, idx)

  col1 <- list()
  if (!is.null(x$sd)) {
    col1$sd <- hist_panel(df, "lsd1", lab$lsd[1], bins)
  }
  col1$r12 <- hist_panel(df, "r12", lab$r12, bins)

  col2 <- scatter_panels(df, lab)

  col3 <- list(plot_ellipses(
    x,
    i = idx[1],
    j = idx[2],
    n = n_ellipses,
    level = level,
    extreme = extreme
  ))

  col4 <- effective_panels(x, extreme = extreme, bins = bins)

  cols <- Filter(length, list(col1, col2, col3, col4))
  cols <- lapply(cols, function(p) patchwork::wrap_plots(p, ncol = 1))
  patchwork::wrap_plots(cols, nrow = 1)
}

#' @rdname plot_tableau
#' @param object A [cov_draws][new_cov_draws] object.
#' @importFrom ggplot2 autoplot
#' @exportS3Method ggplot2::autoplot
autoplot.cov_draws <- function(object, ...) {
  plot_tableau(object, ...)
}

# Panels -------------------------------------------------------------------

cov_palette <- c(mid = "grey55", low = "#2C7BB6", high = "#D7191C")

var_labels <- function(x, idx) {
  l <- x$labels[idx]
  list(
    lsd = sprintf("log(\u03c3[%s])", l),
    r12 = if (length(idx) >= 2L) sprintf("\u03c1[%s,%s]", l[1], l[2]),
    r23 = if (length(idx) >= 3L) sprintf("\u03c1[%s,%s]", l[2], l[3]),
    r13 = if (length(idx) >= 3L) sprintf("\u03c1[%s,%s]", l[1], l[3]),
    r34 = if (length(idx) >= 4L) sprintf("\u03c1[%s,%s]", l[3], l[4])
  )
}

hist_panel <- function(df, var, xlab, bins = 30) {
  ext <- df[df$group != "mid", , drop = FALSE]
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[var]])) +
    ggplot2::geom_histogram(
      bins = bins,
      fill = "grey85",
      colour = "grey40",
      linewidth = 0.2
    )
  if (nrow(ext) > 0L) {
    p <- p +
      ggplot2::geom_rug(
        data = ext,
        ggplot2::aes(x = .data[[var]], colour = .data$group),
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }
  p +
    ggplot2::scale_colour_manual(values = cov_palette) +
    ggplot2::labs(x = xlab, y = "Frequency") +
    ggplot2::theme_bw()
}

scatter_panel <- function(df, xvar, yvar, xlab, ylab) {
  df <- df[order(df$group), , drop = FALSE]
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data[[xvar]], y = .data[[yvar]], colour = .data$group)
  ) +
    ggplot2::geom_point(size = 0.8, alpha = 0.6) +
    ggplot2::scale_colour_manual(values = cov_palette, guide = "none") +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_bw()
}

# Layer 2: the pair types needed under exchangeability
scatter_panels <- function(df, lab) {
  out <- list()
  has_sd <- "lsd1" %in% names(df)
  if (has_sd) {
    out$sd_sd <- scatter_panel(df, "lsd2", "lsd1", lab$lsd[2], lab$lsd[1])
    out$sd_r12 <- scatter_panel(df, "r12", "lsd1", lab$r12, lab$lsd[1])
  }
  if ("r23" %in% names(df)) {
    if (has_sd) {
      out$sd_r23 <- scatter_panel(df, "r23", "lsd1", lab$r23, lab$lsd[1])
    }
    out$r12_r23 <- scatter_panel(df, "r23", "r12", lab$r23, lab$r12)
  }
  if ("r34" %in% names(df)) {
    out$r12_r34 <- scatter_panel(df, "r34", "r12", lab$r34, lab$r12)
  }
  out
}
