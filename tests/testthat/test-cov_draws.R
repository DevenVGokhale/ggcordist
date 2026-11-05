test_that("new_cov_draws validates input", {
  arr <- array(diag(3), dim = c(1, 3, 3))
  expect_s3_class(new_cov_draws(arr), "cov_draws")

  bad <- arr
  bad[1, 1, 2] <- 0.5
  expect_error(new_cov_draws(bad), "symmetric")

  bad <- arr
  bad[1, 2, 2] <- 2
  expect_error(new_cov_draws(bad), "unit diagonal")

  expect_error(new_cov_draws(arr, sd = matrix(-1, 1, 3)), "positive")
  expect_error(new_cov_draws(arr, sd = matrix(1, 2, 3)), "matrix")
  expect_error(new_cov_draws(arr, labels = "a"), "length")
  expect_error(new_cov_draws(matrix(1, 2, 2)), "array")
})

test_that("cov_array reconstructs covariances", {
  x <- r_inv_wishart(5, k = 3)
  cv <- cov_array(x)
  expect_equal(dim(cv), c(5, 3, 3))
  expect_equal(sqrt(cv[, 2, 2]), x$sd[, 2])
  expect_error(cov_array(r_lkj(5, 3)), "no standard deviations")
})

test_that("print and length work", {
  x <- r_lkj(7, 3)
  expect_equal(length(x), 7)
  expect_equal(n_draws(x), 7)
  expect_equal(n_vars(x), 3)
  expect_no_error(capture.output(print(x), type = "message"))
  expect_error(n_draws(1:3), "cov_draws")
})
