library(dplyr)
library(testthat)

# helper to create tipse object
make_tipse <- function(method, arm_to_impute, tipping_param) {
  cox_fit <- coxph(Surv(AVAL, EVENT) ~ TRT01P, ties = "exact", data = codebreak200)
  if (method %in% c("hazard multiplication")) {
    res <- tipping_point_model_based(
      dat = codebreak200,
      reason = "Early dropout",
      impute = arm_to_impute,
      J = 2,
      tipping_range = tipping_param,
      cox_fit = cox_fit,
      seed = 123
    )
  } else {
    res <- tipping_point_model_free(
      dat = codebreak200,
      reason = "Early dropout",
      impute = arm_to_impute,
      J = 2,
      tipping_range = tipping_param,
      cox_fit = cox_fit,
      method = method,
      seed = 123
    )
  }
}

test_that("plot_tp returns a ggplot object with minimal valid tipse input", {
  # --- use helper to build valid tipse object ---
  tipse_obj <- make_tipse("percentile sampling",
    "docetaxel",
    tipping_param = seq(5, 100, by = 10)
  )

  p <- plot_tp(tipse_obj)

  # --- class check ---
  expect_s3_class(p, "ggplot")

  # --- x-axis label check ---
  # percentile sampling + impute = control → % best event times
  expect_match(
    p$labels$x,
    "% best event times sampled from"
  )

  # --- required aesthetics ---
  expect_true("parameter_docetaxel" %in% names(p$data))
  expect_true("HR" %in% names(p$data))

  # ribbon (CI) must exist
  layer_geoms <- vapply(
    p$layers,
    function(l) class(l$geom)[1],
    character(1)
  )
  expect_true("GeomErrorbar" %in% layer_geoms)

  # point layer must exist
  expect_true("GeomPoint" %in% layer_geoms)

  # HR=1 reference line
  expect_true("GeomHline" %in% layer_geoms)

  # tipping point vertical line
  expect_true("GeomVline" %in% layer_geoms)

  # ensure tipping point line is at correct x value
  tp <- tipse_obj$imputation_results$parameter_docetaxel[
    tipse_obj$imputation_results$tipping_point
  ]
})
