library(testthat)
library(dplyr)

# Cox model for testing
cox_fit <- survival::coxph(
  Surv(AVAL, EVENT) ~ TRT01P,
  ties = "exact",
  data = codebreak200
)

# ----------------------------------------------------------------------
# Create model-free tipse object for testing
# ----------------------------------------------------------------------
tp_free <- tipping_point_model_free(
  dat = codebreak200,
  reason = "Early dropout",
  impute = "docetaxel",
  J = 10,
  tipping_range = seq(5, 95, by = 5),
  cox_fit = cox_fit,
  method = "percentile sampling",
  seed = 123
)

# ----------------------------------------------------------------------
# Create landmark sampling tipse
# ----------------------------------------------------------------------
tp_det <- tipping_point_model_free(
  dat = codebreak200,
  reason = "Early dropout",
  impute = "docetaxel",
  J = 10,
  tipping_range = 1:20,
  cox_fit = cox_fit,
  method = "landmark sampling",
  seed = 123
)

# ----------------------------------------------------------------------
# Create model-based tipse object for testing
# ----------------------------------------------------------------------
tp_model <- tipping_point_model_based(
  dat = codebreak200,
  reason = "Early dropout",
  impute = "docetaxel",
  imputation_model = "weibull",
  J = 10,
  tipping_range = seq(0.5, 1, by = 0.25),
  cox_fit = cox_fit,
  seed = 123
)

# ======================================================================
# TESTS
# ======================================================================

test_that("plot.tipse rejects non-tipse objects", {
  expect_error(plot.tipse(list(a = 1)), "tipse")
})


test_that("plot.tipse returns ggplot for Kaplan–Meier (percentile sampling)", {
  expect_message(
    p <- plot(tp_free, type = "Kaplan-Meier"),
    regexp = "Kaplan-Meier"
  )
  expect_s3_class(p, "ggplot")
})



test_that("plot.tipse works for landmark sampling", {
  # KM
  expect_message(
    p1 <- plot(tp_det, type = "Kaplan-Meier"),
    regexp = "landmark sampling"
  )
  expect_s3_class(p1, "ggplot")

  # TP
  expect_message(
    p2 <- plot(tp_det, type = "Tipping Point"),
    regexp = "landmark sampling"
  )
  expect_s3_class(p2, "ggplot")
})


test_that("plot.tipse works for model-based (hazard multiplication)", {
  # KM
  expect_message(
    p1 <- plot(tp_model, type = "Kaplan-Meier"),
    regexp = "hazard multiplication"
  )
  expect_s3_class(p1, "ggplot")

  # TP
  expect_message(
    p2 <- plot(tp_model, type = "Tipping Point"),
    regexp = "hazard multiplication"
  )
  expect_s3_class(p2, "ggplot")
})


test_that("plot.tipse enforces valid `type` argument", {
  expect_error(plot(tp_free, type = "NotAType"))
})


test_that("plot.tipse messages include correct arm and method", {
  expect_message(
    plot(tp_free, type = "Kaplan-Meier"),
    regexp = "docetaxel"
  )
  expect_message(
    plot(tp_det, type = "Kaplan-Meier"),
    regexp = "landmark sampling"
  )
})
