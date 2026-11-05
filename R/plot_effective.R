#' Effective variance and effective dependence plots
#'
#' Layer 4 of Tokuda et al.: histograms of the effective variance (if
#' standard deviations are available) and effective dependence, and the
#' effective dependence of the leading `p x p` submatrix as `p` grows. Each
#' draw gives one line; the shaded band is the pointwise 5%-95% range over all
#' draws.
#'
#' @inheritParams plot_tableau
#' @param n_lines Number of draws shown as lines in the third panel.
#' @return A \pkg{patchwork} object (a ggplot).
#' @seealso [effective_variance()], [effective_dependence()]
#' @export
#' @examples
#' plot_effective(r_lkj(500, k = 10))
plot_effective <- function(x, extreme = 0.05, n_lines = 200, bins = 30, ...) {
  x <- as_cov_draws(x, ...)
  patchwork::wrap_plots(effective_panels(x, extreme, bins, n_lines), ncol = 1)
}

effective_panels <- function(x, extreme = 0.05, bins = 30, n_lines = 200) {
  ed <- effective_dependence(x)
  group <- extreme_group(ed, extreme)
  out <- list()

  if (!is.null(x$sd)) {
    dv <- data.frame(ev = effective_variance(x), group = group)
    out$ev <- hist_panel(dv, "ev", "Effective variance", bins)
  }
  dd <- data.frame(ed = ed, group = group)
  out$ed <- hist_panel(dd, "ed", "Effective dependence", bins)

  path <- effective_dependence_path(x)
  K <- ncol(path)
  band <- data.frame(
    p = seq_len(K),
    lo = apply(
      path,
      2,
      stats::quantile,
      probs = 0.05,
      na.rm = TRUE,
      names = FALSE
    ),
    hi = apply(
      path,
      2,
      stats::quantile,
      probs = 0.95,
      na.rm = TRUE,
      names = FALSE
    )
  )
  s <- thin_index(nrow(path), n_lines)
  lines <- data.frame(
    draw = rep(s, each = K),
    p = rep(seq_len(K), times = length(s)),
    ed = as.vector(t(path[s, , drop = FALSE])),
    group = rep(group[s], each = K)
  )
  lines <- lines[order(lines$group), , drop = FALSE]
  lines$draw <- factor(lines$draw, levels = unique(lines$draw))

  out$path <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = band,
      ggplot2::aes(x = .data$p, ymin = .data$lo, ymax = .data$hi),
      fill = "lightblue",
      alpha = 0.5
    ) +
    ggplot2::geom_line(
      data = lines,
      ggplot2::aes(
        x = .data$p,
        y = .data$ed,
        group = .data$draw,
        colour = .data$group
      ),
      alpha = 0.3,
      linewidth = 0.3
    ) +
    ggplot2::scale_colour_manual(values = cov_palette, guide = "none") +
    ggplot2::labs(
      x = "p (size of leading submatrix)",
      y = "Effective dependence"
    ) +
    ggplot2::theme_bw()
  out
}

#' Eigenvalue distributions
#'
#' Distribution of the ordered log eigenvalues, and of the log condition
#' number, of the draws.
#'
#' @inheritParams plot_tableau
#' @inheritParams cov_eigen
#' @return A \pkg{patchwork} object (a ggplot).
#' @export
#' @examples
#' plot_eigen(r_lkj(300, k = 5))
plot_eigen <- function(x, scale = c("cor", "cov"), bins = 30, ...) {
  x <- as_cov_draws(x, ...)
  scale <- rlang::arg_match(scale)
  ev <- cov_eigen(x, scale)
  K <- ncol(ev)
  ev[ev <= 0 & !is.na(ev)] <- NA_real_
  long <- data.frame(
    rank = factor(rep(seq_len(K), each = nrow(ev))),
    log_eigenvalue = log(as.vector(ev))
  )
  cond <- data.frame(log_cond = log(ev[, 1] / ev[, K]))

  p1 <- ggplot2::ggplot(
    long,
    ggplot2::aes(x = .data$rank, y = .data$log_eigenvalue)
  ) +
    ggplot2::geom_boxplot(outlier.size = 0.5, na.rm = TRUE) +
    ggplot2::labs(x = "Eigenvalue rank", y = "log eigenvalue") +
    ggplot2::theme_bw()
  p2 <- ggplot2::ggplot(cond, ggplot2::aes(x = .data$log_cond)) +
    ggplot2::geom_histogram(
      bins = bins,
      fill = "grey85",
      colour = "grey40",
      linewidth = 0.2,
      na.rm = TRUE
    ) +
    ggplot2::labs(x = "log condition number", y = "Frequency") +
    ggplot2::theme_bw()
  patchwork::wrap_plots(p1, p2, nrow = 1)
}
