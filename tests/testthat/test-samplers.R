test_that("samplers return valid cov_draws", {
  set.seed(1)
  for (x in list(
    r_wishart(30, 4),
    r_inv_wishart(30, 4),
    r_scaled_inv_wishart(30, 4),
    r_separation(30, 4, "lkj"),
    r_separation(30, 4, "inv_wishart")
  )) {
    expect_s3_class(x, "cov_draws")
    expect_equal(dim(x$cor), c(30, 4, 4))
    expect_equal(dim(x$sd), c(30, 4))
  }
  lkj <- r_lkj(30, 4)
  expect_null(lkj$sd)
  expect_equal(dim(lkj$cor), c(30, 4, 4))
})

test_that("LKJ draws are positive definite correlation matrices", {
  set.seed(2)
  x <- r_lkj(100, 6, eta = 0.5)
  ok <- vapply(
    seq_len(100),
    function(s) {
      min(eigen(x$cor[s, , ], symmetric = TRUE, only.values = TRUE)$values) > 0
    },
    logical(1)
  )
  expect_true(all(ok))
})

test_that("LKJ(1) marginals are Beta(K/2, K/2) on (-1, 1)", {
  set.seed(3)
  k <- 4
  x <- r_lkj(3000, k, eta = 1)
  cdf <- function(q) stats::pbeta((q + 1) / 2, k / 2, k / 2)
  for (ij in list(c(1, 2), c(2, 3), c(1, 4), c(3, 4))) {
    p <- stats::ks.test(x$cor[, ij[1], ij[2]], cdf)$p.value
    expect_gt(p, 0.001)
  }
})

test_that("inverse-Wishart with df = k + 1 has uniform correlations", {
  set.seed(4)
  x <- r_inv_wishart(3000, 4, df = 5)
  p <- stats::ks.test(x$cor[, 1, 2], function(q) (q + 1) / 2)$p.value
  expect_gt(p, 0.001)
})

test_that("inverse-Wishart mean is scale / (df - k - 1)", {
  set.seed(5)
  df <- 4 + 10
  sc <- diag(c(1, 2, 3, 4))
  x <- r_inv_wishart(5000, 4, df = df, scale = sc)
  m <- apply(cov_array(x), c(2, 3), mean)
  expect_equal(m, sc / (df - 4 - 1), tolerance = 0.1)
})

test_that("sampler arguments are checked", {
  expect_error(r_lkj(10, 3, eta = -1), "positive")
  expect_error(r_inv_wishart(10, 1), ">= 2")
  expect_error(r_wishart(10, 3, df = 2), "at least")
  expect_error(r_wishart(10, 3, scale = diag(2)), "matrix")
})
