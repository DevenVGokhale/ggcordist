build_ok <- function(p) {
  expect_s3_class(p, "ggplot")
  expect_no_error(suppressWarnings(ggplot2::ggplot_build(p)))
}

test_that("plot_tableau works with and without standard deviations", {
  set.seed(1)
  build_ok(plot_tableau(r_inv_wishart(200, 4)))
  build_ok(plot_tableau(r_lkj(200, 4)))
  build_ok(plot_tableau(r_inv_wishart(200, 3)))
  build_ok(plot_tableau(r_inv_wishart(200, 2)))
  build_ok(plot_tableau(r_inv_wishart(200, 6), idx = c(2, 4, 5, 6)))
  expect_s3_class(ggplot2::autoplot(r_lkj(50, 3)), "ggplot")
  expect_error(plot_tableau(r_inv_wishart(20, 4), idx = c(1, 1, 2)), "distinct")
  expect_error(plot_tableau(r_inv_wishart(20, 4), idx = c(1, 9)), "between")
})

test_that("plot_ellipses validates arguments", {
  set.seed(1)
  x <- r_inv_wishart(100, 4)
  build_ok(plot_ellipses(x, i = 2, j = 3, n = 20, level = 0.9))
  build_ok(plot_ellipses(r_lkj(50, 3)))
  expect_error(plot_ellipses(x, i = 1, j = 1), "different")
  expect_error(plot_ellipses(x, level = 1), "between")
})

test_that("plot_effective and plot_eigen work", {
  set.seed(1)
  build_ok(plot_effective(r_inv_wishart(100, 5)))
  build_ok(plot_effective(r_lkj(100, 5)))
  build_ok(plot_eigen(r_lkj(100, 5)))
  build_ok(plot_eigen(r_inv_wishart(100, 5), scale = "cov"))
})

test_that("plot_pairs works and respects max_vars", {
  set.seed(1)
  build_ok(plot_pairs(r_inv_wishart(100, 3)))
  expect_error(plot_pairs(r_inv_wishart(100, 6)), "max_vars")
})

test_that("plot_compare overlays several sources", {
  set.seed(1)
  a <- r_lkj(100, 3, eta = 1)
  b <- r_lkj(100, 3, eta = 10)
  build_ok(plot_compare(list(a = a, b = b)))
  build_ok(plot_compare(list(r_inv_wishart(100, 4), r_separation(100, 4))))
  expect_error(plot_compare(a), "list")
  expect_error(plot_compare(list(a)), "list")
})

test_that("plot_cor3d requires rgl", {
  skip_if(requireNamespace("rgl", quietly = TRUE))
  expect_error(plot_cor3d(r_lkj(10, 3)), "rgl")
})
