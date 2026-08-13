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

test_that("plot_km works with model-free tipse object and draws expected layers", {
  # --- use helper to build valid tipse object ---
  tipse_obj <- make_tipse("percentile sampling",
    "docetaxel",
    tipping_param = seq(5, 100, by = 10)
  )
  p <- plot_km(tipse_obj)

  # Check ggplot object
  expect_s3_class(p, "ggplot")

  # Extract layers
  geoms <- vapply(p$layers, function(x) class(x$geom)[1], character(1))

  # There should be:
  #  - pooled curves (geom_step)
  #  - highlight curve (geom_step, red)
  #  - original curve (geom_step, orange)
  expect_true("GeomStep" %in% geoms)

  # Build plot for layer data inspection
  built <- ggplot_build(p)

  # Identify pooled KM curves: these come from first layer
  pooled_idx <- which(geoms == "GeomStep")[1]
  pooled_data <- built$data[[pooled_idx]]

  # Check x and y exist (instead of original column names)
  expect_true(all(c("x", "y") %in% names(pooled_data)))

  # Highlight curve (red)
  hl_idx <- which(geoms == "GeomStep")[2]
  hl_data <- built$data[[hl_idx]]

  expect_true(nrow(hl_data) > 0)

  # highlight line must be red
  expect_true(unique(hl_data$colour) == "red")


  # Original KM curve (orange dashed)
  orig_idx <- which(geoms == "GeomStep")[3]
  orig_data <- built$data[[orig_idx]]

  expect_true("surv" %in% names(orig_data) || "y" %in% names(orig_data))

  # Must be orange
  expect_true(unique(orig_data$colour) == "orange")
})
