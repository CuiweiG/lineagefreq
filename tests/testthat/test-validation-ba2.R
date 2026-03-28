# Independent validation: BA.1 → BA.2 transition
# Validates against published literature estimates

test_that("BA.2 growth advantage matches published estimates", {
  data(cdc_ba2_transition, package = "lineagefreq")
  vd  <- lfq_data(cdc_ba2_transition, date = date, lineage = lineage,
                   count = count)
  fit <- fit_model(vd, engine = "mlr", pivot = "BA.1")

  # Generation time 3.2 days for Omicron
  # (Du et al. 2022, doi:10.3201/eid2806.220158)
  ga <- growth_advantage(fit, type = "relative_Rt", generation_time = 3.2)

  ba2_rt <- ga$estimate[ga$lineage == "BA.2"]
  # Published: 1.3-1.5 (Lyngse et al. 2022)
  expect_gt(ba2_rt, 1.1)
  expect_lt(ba2_rt, 2.0)

  # BA.2.12.1 should have higher Rt than BA.2
  ba2121_rt <- ga$estimate[ga$lineage == "BA.2.12.1"]
  expect_gt(ba2121_rt, ba2_rt)

  # Pivot (BA.1) should be exactly 1
  ba1_rt <- ga$estimate[ga$lineage == "BA.1"]
  expect_equal(unname(ba1_rt), 1.0)
})

test_that("BA.2 dataset pipeline is complete", {
  data(cdc_ba2_transition, package = "lineagefreq")
  vd  <- lfq_data(cdc_ba2_transition, date = date, lineage = lineage,
                   count = count)
  fit <- fit_model(vd, engine = "mlr", pivot = "BA.1")

  # All S3 methods work
  expect_true(is.numeric(coef(fit)))
  expect_true(is.matrix(vcov(fit)))
  expect_equal(nrow(glance.lfq_fit(fit)), 1L)

  # Forecast
  fc <- forecast(fit, horizon = 14, n_sim = 200)
  expect_s3_class(fc, "lfq_forecast")

  # MAE < 10% for all lineages (full period)
  resid <- fit$residuals
  for (lin in fit$lineages) {
    r <- resid[resid$.lineage == lin, ]
    mae <- mean(abs(r$.observed - r$.fitted_freq), na.rm = TRUE)
    expect_lt(mae, 0.10, label = paste("MAE for", lin))
  }
})
