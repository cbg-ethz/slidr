context("IH test validation")
library(slidr)

set.seed(156)

n1 <- 4
x1 <- runif(n1,0,1)
n2 <- 13
x2 <- runif(n2,0,1)
n3 <- 45
x3 <- runif(n3,0,1)

test_that("IH returns the correct CDF", {
  expect_equal(IH_CDF(sum(x1),n1), 0.076403, tolerance = 1e-3)
  expect_equal(IH_CDF(sum(x2),n2), 0.523384, tolerance = 1e-3)
  expect_equal(IH_CDF(0.394,4), 0.001, tolerance = 1e-3)
  expect_equal(IH_CDF(sum(x3),n3), pnorm(sum(x3), n3/2, sqrt(n3/12), lower.tail = TRUE))
})
