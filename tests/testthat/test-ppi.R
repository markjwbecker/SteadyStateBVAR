test_that("ppi returns correct structure", {
  result <- ppi(l = 1, u = 3)
  expect_type(result, "list")
  expect_named(result, c("mean", "var"))
})

test_that("ppi calculates mean correctly", {
  result <- ppi(l = 1, u = 3, interval = 0.95)
  expect_equal(result$mean, 2)  # (1 + 3) / 2
})

test_that("ppi with symmetric interval", {
  result <- ppi(l = 0, u = 4, interval = 0.95)
  expect_equal(result$mean, 2)
})

test_that("ppi variance is positive", {
  result <- ppi(l = 1, u = 3)
  expect_gt(result$var, 0)
})

test_that("ppi wider interval has larger variance", {
  result_narrow <- ppi(l = 2, u = 3)
  result_wide <- ppi(l = 1, u = 4)
  expect_lt(result_narrow$var, result_wide$var)
})

test_that("ppi with different confidence levels", {
  result_95 <- ppi(l = 1, u = 3, interval = 0.95)
  result_68 <- ppi(l = 1, u = 3, interval = 0.68)
  expect_gt(result_68$var, result_95$var)
})

test_that("ppi annualized growth rate conversion", {
  result_annual <- ppi(l = 4, u = 8, annualized_growthrate = TRUE, freq = 4)
  # Mean should be (4 + 8) / 2 / 4 = 1.5
  expect_equal(result_annual$mean, 1.5)
})

test_that("ppi with custom frequency", {
  result_monthly <- ppi(l = 12, u = 24, annualized_growthrate = TRUE, freq = 12)
  # Mean should be (12 + 24) / 2 / 12 = 1.5
  expect_equal(result_monthly$mean, 1.5)
})
