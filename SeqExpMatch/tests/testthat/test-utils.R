test_that("logit and inv_logit work correctly", {
  p <- seq(0.1, 0.9, by = 0.1)
  expect_equal(inv_logit(logit(p)), p)
  
  x <- seq(-5, 5, by = 1)
  expect_equal(logit(inv_logit(x)), x)
  
  expect_error(logit(0))
  expect_error(logit(1))
  expect_error(logit(-0.1))
  expect_error(logit(1.1))
})

test_that("sample_mode works correctly", {
  x <- c(1, 1, 2, 3)
  expect_equal(sample_mode(x), 1)
  
  x <- c("a", "b", "b", "c")
  expect_equal(sample_mode(x), "b")
  
  # Ties?
  x <- c(1, 1, 2, 2)
  # Ideally it returns one of them, but check implementation if unsure.
  # Assuming it returns the first mode or similar.
  res <- sample_mode(x)
  expect_true(res %in% c(1, 2))
})
