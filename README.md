# ggcordist

[![R-CMD-check](https://github.com/DevenVGokhale/ggcordist/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DevenVGokhale/ggcordist/actions/workflows/R-CMD-check.yaml)

Visualize distributions of covariance and correlation matrices, such as
posterior draws from a Stan model, following the four-layer method of Tokuda,
Goodrich, Van Mechelen, Gelman and Tuerlinckx, *Visualizing Distributions of
Covariance Matrices*.

## Installation

``` r
# install.packages("pak")
pak::pak("DevenVGokhale/ggcordist")
```

## Usage

``` r
library(ggcordist)

# Draws from a Stan fit (rstan, cmdstanr, brms, or posterior draws)
post <- as_cov_draws(fit, par = "Omega", type = "cor", sd_par = "tau")

plot_tableau(post)             # four-layer display
plot_pairs(post)               # every pair, no exchangeability assumption
plot_ellipses(post, i = 1, j = 2)
plot_effective(post)

# Compare with a reference prior
prior <- r_lkj(2000, k = 4, eta = 2)
plot_compare(list(prior = prior, posterior = post))
```

See `vignette("ggcordist")` for details.
