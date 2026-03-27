test_that("lfq_data constructs correctly from minimal input", {
  d <- data.frame(
    date    = rep(as.Date("2024-01-01") + c(0, 7, 14), each = 2),
    lineage = rep(c("A", "B"), 3),
    count   = c(80, 20, 60, 40, 40, 60)
  )
  x <- lfq_data(d, lineage = lineage, date = date, count = count)

  expect_s3_class(x, "lfq_data")
  expect_equal(nrow(x), 6)
  expect_equal(attr(x, "lineages"), c("A", "B"))
  expect_equal(attr(x, "n_timepoints"), 3L)
  expect_true(all(c(".lineage", ".date", ".count", ".total",
                    ".freq", ".reliable") %in% names(x)))
})

test_that("lfq_data computes frequencies correctly", {
  d <- data.frame(
    date    = rep(as.Date("2024-01-01"), 3),
    lineage = c("A", "B", "C"),
    count   = c(50, 30, 20)
  )
  x <- lfq_data(d, lineage = lineage, date = date, count = count)

  expect_equal(x$.total, rep(100L, 3))
  expect_equal(x$.freq,  c(0.5, 0.3, 0.2))
})

test_that("lfq_data parses character dates", {
  d <- data.frame(
    date    = rep("2024-01-15", 2),
    lineage = c("X", "Y"),
    count   = c(10, 90)
  )
  x <- lfq_data(d, lineage = lineage, date = date, count = count)

  expect_s3_class(x$.date, "Date")
  expect_equal(x$.date[1], as.Date("2024-01-15"))
})

test_that("lfq_data rejects non-existent columns", {
  d <- data.frame(date = Sys.Date(), v = "A", n = 5)
  expect_error(
    lfq_data(d, lineage = missing_col, date = date, count = n)
  )
})

test_that("lfq_data rejects negative counts", {
  d <- data.frame(
    date = Sys.Date(), lineage = "A", count = -5
  )
  expect_error(
    lfq_data(d, lineage = lineage, date = date, count = count),
    "negative"
  )
})

test_that("lfq_data replaces NA counts with 0 and warns", {
  d <- data.frame(
    date    = rep(Sys.Date(), 2),
    lineage = c("A", "B"),
    count   = c(10, NA)
  )
  expect_warning(
    x <- lfq_data(d, lineage = lineage, date = date, count = count),
    "NA"
  )
  expect_equal(x$.count[x$.lineage == "B"], 0L)
})

test_that("lfq_data aggregates duplicates with warning", {
  d <- data.frame(
    date    = rep(Sys.Date(), 4),
    lineage = c("A", "A", "B", "B"),
    count   = c(10, 5, 20, 30)
  )
  expect_warning(
    x <- lfq_data(d, lineage = lineage, date = date, count = count),
    "duplicate"
  )
  expect_equal(nrow(x), 2)
  expect_equal(x$.count[x$.lineage == "A"], 15L)
})

test_that("lfq_data flags low-count time points", {
  d <- data.frame(
    date    = rep(as.Date("2024-01-01") + c(0, 7), each = 2),
    lineage = rep(c("A", "B"), 2),
    count   = c(3, 2, 100, 50)
  )
  x <- lfq_data(d, lineage = lineage, date = date, count = count,
                min_total = 10)

  expect_false(x$.reliable[1])
  expect_true(x$.reliable[3])
})

test_that("is_lfq_data returns correct values", {
  d <- data.frame(
    date = Sys.Date(), lineage = "A", count = 10
  )
  x <- lfq_data(d, lineage = lineage, date = date, count = count)

  expect_true(is_lfq_data(x))
  expect_false(is_lfq_data(d))
  expect_false(is_lfq_data(42))
})

test_that("print.lfq_data runs without error", {
  d <- data.frame(
    date    = rep(Sys.Date(), 2),
    lineage = c("A", "B"),
    count   = c(50, 50)
  )
  x <- lfq_data(d, lineage = lineage, date = date, count = count)
  expect_output(print(x))
})
