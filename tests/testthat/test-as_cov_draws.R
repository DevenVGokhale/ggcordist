make_cov_arr <- function(S = 20, K = 3) {
  set.seed(42)
  cov_array(r_inv_wishart(S, k = K))
}

test_that("array input: cov, cor, and Cholesky types agree", {
  cv <- make_cov_arr()
  from_cov <- as_cov_draws(cv, type = "cov")
  expect_equal(from_cov$sd^2, t(apply(cv, 1, diag)))

  # correlation + sd route
  from_cor <- as_cov_draws(from_cov$cor, type = "cor", sd = from_cov$sd)
  expect_equal(cov_array(from_cor), cv)

  # Cholesky routes
  S <- dim(cv)[1]
  K <- dim(cv)[2]
  chol_cov <- array(NA_real_, dim(cv))
  chol_cor <- array(NA_real_, dim(cv))
  for (s in seq_len(S)) {
    chol_cov[s, , ] <- t(chol(cv[s, , ]))
    chol_cor[s, , ] <- t(chol(from_cov$cor[s, , ]))
  }
  expect_equal(as_cov_draws(chol_cov, type = "chol_cov")$cor, from_cov$cor)
  expect_equal(as_cov_draws(chol_cov, type = "chol_cov")$sd, from_cov$sd)
  expect_equal(as_cov_draws(chol_cor, type = "chol_cor")$cor, from_cov$cor)
})

test_that("cov_draws passes through and bad input errors", {
  x <- r_lkj(3, 3)
  expect_identical(as_cov_draws(x), x)
  expect_error(as_cov_draws("a"), "Don't know")
  expect_error(
    as_cov_draws(array(0, c(2, 2)), type = "cor"),
    "Don't know|dimensions"
  )
  expect_error(as_cov_draws(array(1, c(2, 2, 3))), "square")
})

test_that("posterior draws objects are supported", {
  cv <- make_cov_arr(S = 10, K = 3)
  S <- dim(cv)[1]
  mat <- matrix(NA_real_, S, 9)
  nms <- character(9)
  m <- 0
  for (j in 1:3) {
    for (i in 1:3) {
      m <- m + 1
      mat[, m] <- cv[, i, j]
      nms[m] <- sprintf("Sigma[%d,%d]", i, j)
    }
  }
  colnames(mat) <- nms
  dm <- posterior::as_draws_matrix(mat)

  x <- as_cov_draws(dm, par = "Sigma", type = "cov")
  expect_s3_class(x, "cov_draws")
  expect_equal(cov_array(x), cv)
  expect_equal(x$source, "posterior")

  expect_error(as_cov_draws(dm, par = "nope", type = "cov"), "not found")
  expect_error(as_cov_draws(dm, type = "cov"), "par")
})

test_that("brms-style cor_ and sd_ columns are parsed", {
  S <- 15
  set.seed(1)
  m <- cbind(
    sd_g__Intercept = runif(S, 0.5, 1.5),
    sd_g__x = runif(S, 0.5, 1.5),
    sd_g__z = runif(S, 0.5, 1.5),
    cor_g__Intercept__x = runif(S, -0.3, 0.3),
    cor_g__Intercept__z = runif(S, -0.3, 0.3),
    cor_g__x__z = runif(S, -0.3, 0.3)
  )
  x <- brms_matrix_to_cov_draws(m, "^cor_g__", "^sd_g__", NULL, "posterior")
  expect_equal(x$labels, c("Intercept", "x", "z"))
  expect_equal(x$cor[, 1, 2], unname(m[, "cor_g__Intercept__x"]))
  expect_equal(x$cor[, 2, 1], unname(m[, "cor_g__Intercept__x"]))
  expect_equal(x$cor[, 3, 2], unname(m[, "cor_g__x__z"]))
  expect_equal(x$sd[, 3], unname(m[, "sd_g__z"]))

  bad <- cbind(m, cor_g__Intercept__w = 0.1)
  expect_error(brms_matrix_to_cov_draws(bad, "^cor_g__", "^sd_g__", NULL, "p"), "parse")
  expect_error(
    brms_matrix_to_cov_draws(m[, 1, drop = FALSE], "^cor_g__", "^sd_g__", NULL, "p"),
    "Fewer than two"
  )
})

test_that("brmsfit method requires group or residual", {
  fake <- structure(list(), class = "brmsfit")
  expect_error(as_cov_draws(fake), "group")
})
