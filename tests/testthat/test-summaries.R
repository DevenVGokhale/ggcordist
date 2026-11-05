test_that("effective dependence is 0 for identity and 1 - |R|^(1/k) otherwise", {
  arr <- array(0, c(2, 3, 3))
  arr[1, , ] <- diag(3)
  R <- matrix(c(1, .5, .2, .5, 1, .3, .2, .3, 1), 3)
  arr[2, , ] <- R
  x <- new_cov_draws(arr)
  expect_equal(effective_dependence(x), c(0, 1 - det(R)^(1 / 3)))
})

test_that("effective variance equals det(Sigma)^(1/k)", {
  set.seed(1)
  x <- r_inv_wishart(5, 3)
  cv <- cov_array(x)
  expected <- vapply(1:5, function(s) det(cv[s, , ])^(1 / 3), numeric(1))
  expect_equal(effective_variance(x), expected)
  expect_error(effective_variance(r_lkj(5, 3)), "standard deviations")
})

test_that("effective dependence path matches leading submatrices", {
  set.seed(1)
  x <- r_lkj(4, 5)
  path <- effective_dependence_path(x)
  expect_equal(dim(path), c(4, 5))
  expect_equal(path[, 1], rep(0, 4))
  expect_equal(path[2, 3], 1 - det(x$cor[2, 1:3, 1:3])^(1 / 3))
})

test_that("cov_eigen returns decreasing eigenvalues summing to K for correlations", {
  set.seed(1)
  x <- r_lkj(6, 4)
  ev <- cov_eigen(x)
  expect_true(all(apply(ev, 1, function(r) !is.unsorted(rev(r)))))
  expect_equal(rowSums(ev), rep(4, 6))
  expect_error(cov_eigen(x, "cov"), "standard deviations")
})

test_that("extreme_group flags the requested tails", {
  g <- extreme_group(1:100, 0.1)
  expect_equal(sum(g == "low"), 10)
  expect_equal(sum(g == "high"), 10)
  expect_true(all(extreme_group(1:10, 0) == "mid"))
})

test_that("thin_index is deterministic and bounded", {
  expect_equal(thin_index(5, 10), 1:5)
  i <- thin_index(1000, 100)
  expect_lte(length(i), 100)
  expect_equal(i, thin_index(1000, 100))
})
