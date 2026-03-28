test_that("sequencing_power returns correct structure", {
  sp <- sequencing_power()
  expect_s3_class(sp, "tbl_df")
  expect_equal(nrow(sp), 1L)
  expect_true(sp$required_n > 0)
})

test_that("sequencing_power vectorized", {
  sp <- sequencing_power(current_freq = c(0.01, 0.05, 0.10))
  expect_equal(nrow(sp), 3L)
  expect_true(all(sp$required_n > 0))
})

test_that("sequencing_power rejects bad inputs", {
  expect_error(sequencing_power(current_freq = 0))
  expect_error(sequencing_power(current_freq = 1))
  expect_error(sequencing_power(target_precision = -1))
})
