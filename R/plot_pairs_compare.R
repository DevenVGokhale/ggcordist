#' Pairs plot of all log standard deviations and correlations
#'
#' A scatterplot-matrix with histograms on the diagonal for every
#' `log(sd_i)` and every correlation `rho_ij`. Unlike [plot_tableau()] this
#' does not assume exchangeability, so it is suited to posterior draws of
#' small matrices.
#'
#' @inheritParams plot_tableau
#' @param max_vars Maximum number of variables (`K` log standard deviations
#'   plus `K(K-1)/2` correlations) to plot, as the number of panels grows
#'   with its square.
#' @return A \pkg{patchwork} object (a ggplot).
#' @export
#' @examples
#' plot_pairs(r_inv_wishart(300, k = 3))
plot_pairs <- function(x, extreme = 0.05, bins = 30, max_vars = 12, ...) {
  x <- as_cov_draws(x, ...)
  K <- n_vars(x)
  lower <- which(lower.tri(diag(K)), arr.ind = TRUE)
  lower <- lower[order(lower[, 2], lower[, 1]), , drop = FALSE]

  vars <- list()
  if (!is.null(x$sd)) {
    for (i in seq_len(K)) {
      vars[[length(vars) + 1L]] <- list(
        name = sprintf("log(\u03c3[%s])", x$labels[i]),
        values = log(x$sd[, i])
      )
    }
  }
  for (r in seq_len(nrow(lower))) {
    i <- lower[r, 1]
    j <- lower[r, 2]
    vars[[length(vars) + 1L]] <- list(
      name = sprintf("\u03c1[%s,%s]", x$labels[j], x$labels[i]),
      values = x$cor[, i, j]
    )
  }
  m <- length(vars)
  if (m > max_vars) {
    cli::cli_abort(c(
      "{m} variables to plot exceeds {.arg max_vars} = {max_vars}.",
      i = "Reduce the dimension or raise {.arg max_vars}."
    ))
  }

  group <- extreme_group(effective_dependence(x), extreme)
  df <- data.frame(group = group)
  for (v in seq_len(m)) {
    df[[paste0("v", v)]] <- vars[[v]]$values
  }

  panels <- vector("list", m * m)
  for (r in seq_len(m)) {
    for (c in seq_len(m)) {
      k <- (r - 1L) * m + c
      panels[[k]] <- if (r == c) {
        hist_panel(df, paste0("v", r), vars[[r]]$name, bins)
      } else if (r > c) {
        scatter_panel(
          df,
          paste0("v", c),
          paste0("v", r),
          vars[[c]]$name,
          vars[[r]]$name
        )
      } else {
        patchwork::plot_spacer()
      }
    }
  }
  patchwork::wrap_plots(panels, ncol = m)
}

#' Compare several distributions of covariance matrices
#'
#' Overlays marginal densities and scatterplots for several
#' [cov_draws][new_cov_draws] objects, e.g. prior versus posterior draws, or a
#' posterior against a reference prior from [r_lkj()].
#'
#' @param x A named list of [cov_draws][new_cov_draws] objects (or objects
#'   accepted by [as_cov_draws()] without further arguments). Names are used as
#'   legend labels; unnamed elements use their `source`.
#' @param idx Up to four variable indices as in [plot_tableau()].
#' @return A \pkg{patchwork} object (a ggplot).
#' @export
#' @examples
#' prior <- r_lkj(500, k = 3, eta = 1)
#' tight <- r_lkj(500, k = 3, eta = 10)
#' plot_compare(list(`LKJ(1)` = prior, `LKJ(10)` = tight))
plot_compare <- function(x, idx = NULL) {
  if (inherits(x, "cov_draws") || !is.list(x) || length(x) < 2L) {
    cli::cli_abort(
      "{.arg x} must be a list of at least two {.cls cov_draws} objects."
    )
  }
  x <- lapply(x, as_cov_draws)
  nms <- names(x)
  if (is.null(nms)) {
    nms <- rep("", length(x))
  }
  src <- vapply(x, function(d) d$source, character(1))
  nms[nms == ""] <- src[nms == ""]
  nms <- make.unique(nms)

  K <- min(vapply(x, n_vars, integer(1)))
  idx <- idx %||% seq_len(min(K, 4L))
  frames <- Map(
    function(d, nm) {
      f <- cov_frame(d, idx, extreme = 0)
      f$source <- nm
      f
    },
    x,
    nms
  )
  common <- Reduce(intersect, lapply(frames, names))
  df <- do.call(rbind, lapply(frames, function(f) f[, common, drop = FALSE]))
  df$source <- factor(df$source, levels = nms)
  lab <- var_labels(x[[1]], idx)

  dens <- function(var, xlab) {
    ggplot2::ggplot(
      df,
      ggplot2::aes(x = .data[[var]], colour = .data$source, fill = .data$source)
    ) +
      ggplot2::geom_density(alpha = 0.25) +
      ggplot2::labs(x = xlab, y = "Density", colour = NULL, fill = NULL) +
      ggplot2::theme_bw()
  }
  sc <- function(xv, yv, xlab, ylab) {
    ggplot2::ggplot(
      df,
      ggplot2::aes(x = .data[[xv]], y = .data[[yv]], colour = .data$source)
    ) +
      ggplot2::geom_point(size = 0.7, alpha = 0.4) +
      ggplot2::labs(x = xlab, y = ylab, colour = NULL) +
      ggplot2::theme_bw()
  }

  panels <- list()
  if ("lsd1" %in% common) {
    panels$lsd1 <- dens("lsd1", lab$lsd[1])
  }
  panels$r12 <- dens("r12", lab$r12)
  panels$ed <- dens("ed", "Effective dependence")
  if ("lsd1" %in% common) {
    panels$sd_r12 <- sc("r12", "lsd1", lab$r12, lab$lsd[1])
  }
  if ("r23" %in% common) {
    panels$r12_r23 <- sc("r23", "r12", lab$r23, lab$r12)
  }
  patchwork::wrap_plots(panels, guides = "collect")
}
