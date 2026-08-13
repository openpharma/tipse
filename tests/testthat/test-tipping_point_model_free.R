library(testthat)
library(survival)
library(dplyr)

# Cox model for testing
cox_fit <- coxph(Surv(AVAL, EVENT) ~ TRT01P, ties = "exact", data = codebreak200)

test_that("tipping_point_model_free returns correct structure for percentile sampling", {
  expect_warning(res <- tipping_point_model_free(
    dat = codebreak200,
    reason = "Early dropout",
    impute = "sotorasib",
    J = 2,
    tipping_range = c(20, 40),
    cox_fit = cox_fit,
    method = "percentile sampling",
    seed = 123
  ), "Tipping point not found, please check 'tipping_range'.")

  expect_s3_class(res, "tipse")
  expect_true(all(c(
    "original_data", "original_HR", "reason_to_impute", "arm_to_impute",
    "method_to_impute", "imputation_results", "imputation_data"
  ) %in% names(res)))

  # Check imputation_results structure (single-arm uses arm-named parameter column)
  expect_true(all(c("HR", "HR_upperCI", "HR_lowerCI", "parameter_sotorasib", "tipping_point") %in% names(res$imputation_results)))
})


test_that("tipping_point_model_free throws error for invalid cox_fit", {
  expect_error(
    tipping_point_model_free(
      dat = codebreak200,
      reason = "Early dropout",
      impute = "sotorasib",
      J = 2,
      tipping_range = c(20, 40),
      cox_fit = "not_a_coxph",
      method = "percentile sampling"
    ),
    "Argument 'cox_fit' must be a valid cox model object"
  )
})

test_that("tipping_point_model_free throws error for invalid impute arm", {
  expect_error(
    tipping_point_model_free(
      dat = codebreak200,
      reason = "Early dropout",
      impute = "invalid_arm",
      J = 2,
      tipping_range = c(20, 40),
      cox_fit = cox_fit,
      method = "percentile sampling"
    ),
    "Argument 'impute' must be one or both of the arms"
  )
})

test_that("tipping_point_model_free throws error when reason is empty", {
  expect_error(
    tipping_point_model_free(
      dat = codebreak200,
      reason = character(0),
      impute = "sotorasib",
      J = 2,
      tipping_range = c(20, 40),
      cox_fit = cox_fit,
      method = "percentile sampling"
    ),
    "Argument 'reason' must specify at least one censoring reason"
  )
})

test_that("tipping_point_model_free throws error when percentile_range is missing for percentile sampling", {
  expect_error(
    tipping_point_model_free(
      dat = codebreak200,
      reason = "Early dropout",
      impute = "sotorasib",
      J = 2,
      tipping_range = NULL,
      cox_fit = cox_fit,
      method = "percentile sampling"
    ),
    "tipping_range"
  )
})

test_that("tipping_point_model_free throws error when npts_range is missing for landmark sampling", {
  expect_error(
    tipping_point_model_free(
      dat = codebreak200,
      reason = "Early dropout",
      impute = "sotorasib",
      J = 2,
      tipping_range = NULL,
      cox_fit = cox_fit,
      method = "landmark sampling"
    ),
    "tipping_range"
  )
})

test_that("tipping_point_model_free sets tipping_point flag correctly", {
  res <- tipping_point_model_free(
    dat = codebreak200,
    reason = "Early dropout",
    impute = "docetaxel",
    J = 2,
    tipping_range = seq(10, 90, by = 10),
    cox_fit = cox_fit,
    method = "percentile sampling",
    seed = 123
  )

  expect_true("tipping_point" %in% names(res$imputation_results))
  expect_true(any(res$imputation_results$tipping_point %in% c(TRUE, FALSE)))
})

test_that("tipping_point_model_free respects seed for reproducibility", {
  res1 <- tipping_point_model_free(
    dat = codebreak200,
    reason = "Early dropout",
    impute = "docetaxel",
    J = 2,
    tipping_range = seq(10, 90, by = 10),
    cox_fit = cox_fit,
    method = "percentile sampling",
    seed = 123
  )

  res2 <- tipping_point_model_free(
    dat = codebreak200,
    reason = "Early dropout",
    impute = "docetaxel",
    J = 2,
    tipping_range = seq(10, 90, by = 10),
    cox_fit = cox_fit,
    method = "percentile sampling",
    seed = 123
  )

  expect_equal(res1$imputation_results, res2$imputation_results)
})

test_that("tipping_point_model_free two-arm percentile: correct structure and arm-named columns", {
  expect_warning(
    res <- tipping_point_model_free(
      dat = codebreak200,
      reason = "Early dropout",
      impute = c("docetaxel", "sotorasib"),
      J = 2,
      tipping_range = c(20, 40),
      cox_fit = cox_fit,
      method = "percentile sampling",
      seed = 123
    ),
    "Tipping point not found"
  )

  expect_s3_class(res, "tipse")
  expect_equal(res$arm_to_impute, c("docetaxel", "sotorasib"))

  # arm-named parameter columns
  expect_true("parameter_docetaxel" %in% names(res$imputation_results))
  expect_true("parameter_sotorasib" %in% names(res$imputation_results))
  expect_false("parameter" %in% names(res$imputation_results))
  expect_false("parameter_1" %in% names(res$imputation_results))

  # Grid: 2 x 2 = 4 rows
  expect_equal(nrow(res$imputation_results), 4L)
})

test_that("tipping_point_model_free two-arm landmark: correct structure and arm-named columns", {
  expect_warning(
    res <- tipping_point_model_free(
      dat = codebreak200,
      reason = "Early dropout",
      impute = c("docetaxel", "sotorasib"),
      J = 2,
      tipping_range = c(1, 2),
      cox_fit = cox_fit,
      method = "landmark sampling",
      seed = 123
    ),
    "Tipping point not found"
  )

  expect_s3_class(res, "tipse")
  expect_true("parameter_docetaxel" %in% names(res$imputation_results))
  expect_true("parameter_sotorasib" %in% names(res$imputation_results))
  expect_equal(nrow(res$imputation_results), 4L)
})

test_that("tipping_point_model_free two-arm: per-arm list tipping_range produces Cartesian grid", {
  expect_warning(
    res <- tipping_point_model_free(
      dat = codebreak200,
      reason = "Early dropout",
      impute = c("docetaxel", "sotorasib"),
      J = 2,
      tipping_range = list(c(10, 20, 30), c(40, 50)),
      cox_fit = cox_fit,
      method = "percentile sampling",
      seed = 123
    ),
    "Tipping point not found"
  )

  # 3 x 2 = 6 rows
  expect_equal(nrow(res$imputation_results), 6L)
  expect_equal(length(unique(res$imputation_results$parameter_docetaxel)), 3L)
  expect_equal(length(unique(res$imputation_results$parameter_sotorasib)), 2L)
})

test_that("tipping_point_model_free two-arm: km_data_list names use pipe-separated parameters", {
  expect_warning(
    res <- tipping_point_model_free(
      dat = codebreak200,
      reason = "Early dropout",
      impute = c("docetaxel", "sotorasib"),
      J = 2,
      tipping_range = c(20, 40),
      cox_fit = cox_fit,
      method = "percentile sampling",
      seed = 123
    ),
    "Tipping point not found"
  )

  expect_true(all(grepl("|", names(res$imputation_data), fixed = TRUE)))
})

test_that("tipping_point_model_free two-arm: errors on duplicate arms", {
  expect_error(
    tipping_point_model_free(
      dat = codebreak200,
      reason = "Early dropout",
      impute = c("docetaxel", "docetaxel"),
      J = 2,
      tipping_range = c(20, 40),
      cox_fit = cox_fit,
      method = "percentile sampling"
    ),
    "must not contain the same arm twice"
  )
})

test_that("tipping_point_model_free two-arm: list tipping_range rejected for single arm", {
  expect_error(
    tipping_point_model_free(
      dat = codebreak200,
      reason = "Early dropout",
      impute = "docetaxel",
      J = 2,
      tipping_range = list(c(10, 20), c(30, 40)),
      cox_fit = cox_fit,
      method = "percentile sampling"
    ),
    "only allowed when two arms"
  )
})

test_that("tipping_point_model_free two-arm: reproducible with seed", {
  args <- list(
    dat = codebreak200, reason = "Early dropout",
    impute = c("docetaxel", "sotorasib"),
    J = 2, tipping_range = c(20, 40), cox_fit = cox_fit,
    method = "percentile sampling", seed = 42
  )
  expect_warning(res1 <- do.call(tipping_point_model_free, args), "Tipping point not found")
  expect_warning(res2 <- do.call(tipping_point_model_free, args), "Tipping point not found")
  expect_equal(res1$imputation_results, res2$imputation_results)
})
