test_that("RemoveSeasonal computes seasonal means and anomalies", {
  y <- c(2, 4, 6, 8, 10, 12)
  t <- c(1, 2, 1, 2, 1, 2)
  out <- rpkg.template:::RemoveSeasonal(y, t)
  expect_equal(out$mu, c(6, 8))
  expect_equal(out$anomaly, c(-4, -4, 0, 0, 4, 4))
})

test_that("RemoveSeasonal orders mu by sorted unique t_mmdd", {
  y <- c(4, 2, 8, 6)
  t <- c(2, 1, 2, 1)
  expect_equal(rpkg.template:::RemoveSeasonal(y, t)$mu, c(4, 6))
})

test_that("RemoveSeasonal handles NA", {
  y <- c(2, NA, 6, 8, 10, 12)
  t <- c(1, 2, 1, 2, 1, 2)
  out <- rpkg.template:::RemoveSeasonal(y, t)
  expect_equal(out$mu, c(6, 10))
  expect_equal(out$anomaly, c(-4, NA, 0, -2, 4, 2))
})

test_that("LS_OneBreak recovers coefficients exactly", {
  x <- 1:10
  y <- 2 + 3 * x + 5 * (x > 4)
  out <- rpkg.template::LS_OneBreak(y, x, Ic = 4)
  expect_equal(as.numeric(out$coef), c(2, 3, 5))
  expect_equal(out$fitted, y)
  expect_equal(out$SSE, 0)
})

test_that("LS_OneBreak without break equals simple lm", {
  x <- 1:5
  y <- 1 + 2 * x
  out <- rpkg.template::LS_OneBreak(y, x, Ic = NA)
  expect_equal(as.numeric(out$coef), c(1, 2))
  expect_equal(out$SSE, 0)
})

test_that("LS_MultiBreaks recovers coefficients exactly", {
  x <- 1:9
  y <- 1 + 2 * (x - mean(x)) + c(0, 0, 0, 3, 3, 3, -2, -2, -2)
  out <- rpkg.template:::LS_MultiBreaks(y, x, Ips = c(3, 6, 9))
  expect_equal(as.numeric(out$coef), c(1, 2, 3, -2))
  expect_equal(out$fitted, y)
  expect_equal(out$SSE, 0)
})

test_that("lm_solve only.SSE matches full SSE", {
  x <- 1:10
  design <- cbind(1, x, c(rep(0, 4), rep(1, 6)))
  y <- 2 + 3 * x + 5 * c(rep(0, 4), rep(1, 6))
  full <- rpkg.template:::lm_solve(design, y)
  sse <- rpkg.template:::lm_solve(design, y, only.SSE = TRUE)
  expect_equal(sse$SSE, full$SSE)
  expect_equal(sse$SSE, sum(full$resi^2))
})
