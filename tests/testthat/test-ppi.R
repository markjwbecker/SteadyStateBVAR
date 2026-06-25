test_that("ppi validates lower < upper", {
  expect_error(ppi(l = 3, u = 1), "lower bound must be less than upper")
})

test_that("ppi validates interval is between 0 and 1", {
  expect_error(ppi(l = 1, u = 3, interval = 1.5), "interval must be between 0 and 1")
  expect_error(ppi(l = 1, u = 3, interval = 0), "interval must be between 0 and 1")
})

test_that("ppi validates freq is positive", {
  expect_error(ppi(l = 1, u = 3, annualized_growthrate = TRUE, freq = -1), 
               "freq must be positive")
})

test_that("ppi returns list with mean and var", {
  result <- ppi(l = 1, u = 3)
  
  expect_type(result, "list")
  expect_named(result, c("mean", "var"))
  expect_type(result$mean, "double")
  expect_type(result$var, "double")
})

test_that("ppi calculates mean correctly", {
  result <- ppi(l = 1, u = 3, interval = 0.95)
  
  expect_equal(result$mean, 2)
})

test_that("ppi variance is always positive", {
  result <- ppi(l = 1, u = 3)
  expect_gt(result$var, 0)
})

test_that("ppi wider interval has larger variance", {
  result_narrow <- ppi(l = 2, u = 3)
  result_wide <- ppi(l = 1, u = 4)
  
  expect_lt(result_narrow$var, result_wide$var)
})

test_that("ppi with different confidence intervals", {
  result_68 <- ppi(l = 1, u = 3, interval = 0.68)
  result_95 <- ppi(l = 1, u = 3, interval = 0.95)
  result_99 <- ppi(l = 1, u = 3, interval = 0.99)
  
  expect_equal(result_68$mean, result_95$mean)
  expect_equal(result_95$mean, result_99$mean)
  
  # Lower confidence = larger variance (for same bounds)
  expect_gt(result_68$var, result_95$var)
  expect_gt(result_95$var, result_99$var)
})

test_that("ppi with annualized growth rate conversion", {
  result_annual <- ppi(l = 4, u = 8, annualized_growthrate = TRUE, freq = 4)
  
  expect_equal(result_annual$mean, 1.5)
})

test_that("ppi annualized conversion scales variance correctly", {
  result_level <- ppi(l = 4, u = 8, annualized_growthrate = FALSE)
  result_annual <- ppi(l = 4, u = 8, annualized_growthrate = TRUE, freq = 4)
  
  expect_equal(result_annual$var, result_level$var / 16)
})

test_that("ppi with custom frequency", {
  result_monthly <- ppi(l = 12, u = 24, annualized_growthrate = TRUE, freq = 12)
  result_quarterly <- ppi(l = 12, u = 24, annualized_growthrate = TRUE, freq = 4)
  
  expect_equal(result_monthly$mean, 1.5)
  expect_equal(result_quarterly$mean, 4.5)
  
  expect_lt(result_monthly$var, result_quarterly$var)
})

test_that("ppi symmetric intervals", {
  result_pos <- ppi(l = 0, u = 4)
  expect_equal(result_pos$mean, 2)
  
  result_shift <- ppi(l = 3, u = 7)
  expect_equal(result_shift$mean, 5)
})

test_that("ppi with very tight interval", {
  result <- ppi(l = 2.99, u = 3.01, interval = 0.95)
  
  expect_equal(result$mean, 3)
  expect_lt(result$var, 0.0001)
})

test_that("ppi with very wide interval", {
  result <- ppi(l = 0, u = 100, interval = 0.95)
  
  expect_equal(result$mean, 50)
  expect_gt(result$var, 100)
})

test_that("ppi realistic inflation prior", {
  result <- ppi(l = 1, u = 3, interval = 0.90)
  
  expect_equal(result$mean, 2)
  expect_gt(result$var, 0)
})

test_that("ppi realistic growth rate prior annualized", {
  result <- ppi(l = 2, u = 4, annualized_growthrate = TRUE, freq = 4)
  
  expect_equal(result$mean, 0.75)
})