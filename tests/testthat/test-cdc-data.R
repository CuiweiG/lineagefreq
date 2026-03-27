test_that("cdc_sarscov2_jn1 dataset loads correctly", {
  data(cdc_sarscov2_jn1, package = "lineagefreq")

  expect_true(is.data.frame(cdc_sarscov2_jn1))
  expect_true(all(c("date", "lineage", "count") %in%
                    names(cdc_sarscov2_jn1)))
  expect_true(inherits(cdc_sarscov2_jn1$date, "Date"))
  expect_true(is.numeric(cdc_sarscov2_jn1$count))
  expect_true(all(cdc_sarscov2_jn1$count >= 0))
  expect_true("JN.1" %in% cdc_sarscov2_jn1$lineage)
  expect_true("proportion" %in% names(cdc_sarscov2_jn1))
})

test_that("full pipeline runs on CDC data", {
  data(cdc_sarscov2_jn1, package = "lineagefreq")

  vd <- lfq_data(cdc_sarscov2_jn1,
                 date = date, lineage = lineage, count = count)
  expect_s3_class(vd, "lfq_data")

  fit <- fit_model(vd, engine = "mlr", pivot = "XBB.1.5")
  expect_s3_class(fit, "lfq_fit")

  ga <- growth_advantage(fit, type = "relative_Rt",
                         generation_time = 5)
  expect_s3_class(ga, "tbl_df")

  # JN.1 should have relative Rt > 1 (vs XBB.1.5)
  jn1 <- ga[ga$lineage == "JN.1", ]
  if (nrow(jn1) > 0) {
    expect_gt(jn1$estimate[1], 1.0)
  }

  fc <- forecast(fit, horizon = 14, n_sim = 200)
  expect_s3_class(fc, "lfq_forecast")
})
