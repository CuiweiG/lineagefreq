test_that("as_lfq_data.data.frame works", {
  df <- data.frame(
    date    = rep(as.Date("2024-01-01") + c(0, 7), each = 2),
    lineage = rep(c("A", "B"), 2),
    count   = c(80, 20, 60, 40)
  )
  x <- as_lfq_data(df, lineage = lineage, date = date, count = count)
  expect_s3_class(x, "lfq_data")
})

test_that("as_lfq_data.lfq_data is identity", {
  d <- data.frame(date = Sys.Date(), lineage = "A", count = 10)
  x <- lfq_data(d, lineage = lineage, date = date, count = count)
  y <- as_lfq_data(x)
  expect_identical(x, y)
})

test_that("lfq_engines returns tibble with all engines", {
  eng <- lfq_engines()
  expect_s3_class(eng, "tbl_df")
  expect_equal(nrow(eng), 5L)
  expect_true(all(c("mlr", "hier_mlr", "piantham", "fga", "garw")
                  %in% eng$engine))
  expect_true(all(c("engine", "type", "time_varying",
                    "available", "description") %in% names(eng)))
  # Core engines always available
  expect_true(all(eng$available[eng$type == "frequentist"]))
})

test_that("lfq_version returns expected structure", {
  v <- lfq_version()
  expect_type(v, "list")
  expect_true(all(c("version", "r_version", "stan_available",
                    "engines") %in% names(v)))
  expect_true(is.character(v$version))
  expect_true(is.logical(v$stan_available))
  expect_true("mlr" %in% v$engines)
})
