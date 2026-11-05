# ggcordist: plan and progress

Last updated: 2026-11-05

Goal: an R package that visualizes draws of correlation / covariance matrices
(e.g. from a Stan fit) using the four-layer method of Tokuda, Goodrich, Van
Mechelen, Gelman and Tuerlinckx, *Visualizing Distributions of Covariance
Matrices* (`reference/428c03ac-7fa3-4916-81bc-5d3ad9aa5aa7.pdf`).

Status key: `[x]` done, `[~]` done but not fully verified, `[ ]` not started.

## Decisions

- Inputs: rstan, cmdstanr, brms, posterior `draws_*`, plain arrays; all
  converted to one internal `cov_draws` object.
- Graphics: ggplot2 + patchwork; rgl only for the 3D plot.
- The paper's grid assumes exchangeable variables. Posterior draws usually are
  not, so `plot_tableau()` takes `idx` and `plot_pairs()` shows every pair.
- Draws in the lowest / highest 5% of effective dependence are coloured
  blue / red (as in the paper).

## Milestones

### 0. Method verification

- [x] Extract paper text with pdftools (`reference/tokuda_text.txt`, git-ignored)
- [x] Confirm the four layers, ellipse level (50%), ~100 ellipses, L = 1000 draws
- [x] Confirm effective variance / dependence definitions (eqs. 4, 5)
- [ ] Read the rest of the paper (Sections 3.5 onward: Wishart, appendix) and
      check against the samplers and plots
- [ ] Compare each vignette figure side by side with the paper's Figures 1-8

### 1. Package scaffold

- [x] DESCRIPTION, MIT licence, NAMESPACE via roxygen2, tests/, vignettes/
- [x] `.Rbuildignore`, `.gitignore`, NEWS.md, README.md
- [x] GitHub Actions workflow for R CMD check (green on all 5 jobs, 2026-11-05)
- [ ] `usethis::use_pkgdown()` site (optional)
- [x] Replace `<github-user>` placeholder in README install instructions
- [x] Initial git commit and push to https://github.com/DevenVGokhale/ggcordist

### 2. Data container and inputs

- [x] `new_cov_draws()` with validation; `print()`, `length()`, `n_draws()`, `n_vars()`
- [x] `cov_array()` reconstructs covariances
- [x] `as_cov_draws()` for arrays, posterior draws; types `cor`, `cov`, `chol_cor`, `chol_cov`
- [~] `as_cov_draws()` for `stanfit` (rstan): code written, never run on a real fit
- [~] `as_cov_draws()` for `CmdStanMCMC` (cmdstanr): code written, never run on a real fit
- [~] `as_cov_draws()` for `brmsfit`: name parsing tested on mocked columns only
- [ ] Tests with small real Stan / brms fits (skipped when Stan is unavailable)
- [ ] Support for multiple chains / `chain` column in `cov_draws`
- [ ] Decide whether to keep `residual = TRUE` brms path (untested)

### 3. Summaries

- [x] `effective_dependence()`, `effective_variance()`, `effective_dependence_path()`
- [x] `cov_eigen()`
- [ ] Optional Fisher z transform for correlations
- [ ] Tidy long-format helper (`cov_quantities()`) from the original plan

### 4. Plots

- [x] `plot_tableau()` and `autoplot()` method (four layers)
- [x] `plot_ellipses()`
- [x] `plot_effective()`
- [x] `plot_eigen()`
- [x] `plot_pairs()`
- [x] `plot_compare()` (prior vs posterior overlay)
- [~] `plot_cor3d()`: needs `rgl`, which is not installed; untested
- [ ] Check `plot_tableau()` layout at paper size (14 x 9 in); tidy axis labels and
      panel proportions (panels are cramped in a small viewer)
- [ ] Add the reference line for a uniform correlation density (green line in paper Fig. 1)
- [ ] Add the `"pooled"` mode that stacks all index pairs of the same type for
      exchangeable priors (in the original plan, not implemented)
- [ ] Faceted `plot_compare()` for the full tableau, not only marginals

### 5. Reference samplers

- [x] `r_wishart()`, `r_inv_wishart()`, `r_scaled_inv_wishart()`
- [x] `r_lkj()` (C-vine method), `r_separation()` (LKJ or inverse-Wishart correlation)
- [x] Statistical checks: uniform marginal for inv-Wishart with df = k + 1;
      Beta(k/2, k/2) marginal for LKJ(1); inv-Wishart mean
- [ ] Check the correlation marginal of `r_separation(cor_dist = "inv_wishart")`
      against the paper's Eq. 7 kernel

### 6. Tests

- [x] 117 expectations across 5 test files, all passing
- [ ] vdiffr snapshot tests for each plot (vdiffr is not installed)
- [ ] Coverage report (`covr`)

### 7. Documentation

- [x] roxygen2 docs and examples for every exported function
- [x] Vignette `ggcordist.Rmd` (reference distributions, Stan usage, comparison)
- [ ] Run the Stan / brms vignette chunks (currently `eval = FALSE`) on a real model
- [ ] Second vignette: hierarchical model with varying slopes, prior vs posterior

### 8. Release checks

- [x] `R CMD check --no-manual` on Windows: Status OK (2026-11-05), run with
      `_R_CHECK_FORCE_SUGGESTS_=false`
- [ ] Check with Suggests installed (needs `rgl`, `vdiffr`)
- [x] Check on Linux / macOS (GitHub Actions: macOS, Windows, Ubuntu release/devel/oldrel-1 all pass)
- [ ] `devtools::spell_check()`, `urlchecker`

## Known issues / ideas

- Axis labels must stay ASCII: Greek characters fail on CI (C locale, pdf device).
  A plotmath-based label (`expression(rho[12])`) would restore Greek letters safely.
- `Rscript` is not on the PATH in the integrated terminal; run checks from the R
  session or add R's `bin` folder to PATH.
- The paper's package, VisCov, uses base graphics; ggcordist output will look different.
- Possible extension: plot a posterior against its prior directly from a Stan fit
  with `sample_prior`-style draws.

## Log

- 2026-11-05: Scaffolded package, implemented all functions above, 117 tests
  pass, `R CMD check` Status OK. Plan saved to `dev/PLAN.md`.
- 2026-11-05: Published to GitHub (`main`). `reference/` (paper PDF and text)
  is git-ignored for copyright reasons.
- 2026-11-05: Added `.github/workflows/R-CMD-check.yaml` (macOS, Windows, Ubuntu
  release/devel/oldrel-1; hard dependencies only, so Stan packages and rgl are
  not installed on CI).
- 2026-11-05: First CI run failed on all platforms: Greek letters in axis labels
  could not be drawn in the CI's non-UTF-8 locale. Switched labels to ASCII
  (`rho[1,2]`, `log(sd[1])`); all 5 jobs now pass.
