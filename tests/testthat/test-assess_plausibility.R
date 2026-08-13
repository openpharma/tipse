library(testthat)
library(dplyr)

# helper to create tipse object
make_tipse <- function(dat, method, arm_to_impute, reason, tipping_param) {
  cox_fit <- coxph(Surv(AVAL, EVENT) ~ TRT01P, ties = "exact", data = dat)
  if (method == "hazard multiplication") {
    res <- tipping_point_model_based(
      dat = dat,
      reason = reason,
      impute = arm_to_impute,
      J = 2,
      tipping_range = tipping_param,
      cox_fit = cox_fit,
      seed = 123
    )
  } else {
    res <- tipping_point_model_free(
      dat = dat,
      reason = reason,
      impute = arm_to_impute,
      J = 2,
      tipping_range = tipping_param,
      cox_fit = cox_fit,
      method = method,
      seed = 123
    )
  }
}

test_that("assess_plausibility works for control arm with percentile sampling", {
  tipse <- make_tipse(codebreak200, "percentile sampling", "docetaxel", "Early dropout", tipping_param = seq(5, 100, by = 10))

  expect_output(
    assess_plausibility(tipse),
    regexp = "median duration of follow-up"
  )
})

test_that("assess_plausibility works for control arm with hazard multiplication (deflation)", {
  tipse <- suppressWarnings(make_tipse(codebreak200, "hazard multiplication", "docetaxel", "Early dropout", tipping_param = seq(0.1, 0.9, by = 0.2)))

  expect_output(
    assess_plausibility(tipse),
    regexp = "HR between imputed set in docetaxel arm"
  )
})

test_that("assess_plausibility works for treatment arm with landmark sampling", {
  tipse <- make_tipse(extenet, "landmark sampling", "neratinib", "Lost to follow-up", tipping_param = seq(1, 20, by = 1))

  expect_output(
    assess_plausibility(tipse),
    regexp = "event rate in imputed set in neratinib"
  )
})

test_that("assess_plausibility works for treatment arm with hazard multiplication (inflation)", {
  tipse <- make_tipse(extenet, "hazard multiplication", "neratinib", "Lost to follow-up", tipping_param = seq(1, 4, by = 0.5))

  expect_output(
    assess_plausibility(tipse),
    regexp = "HR between imputed set in neratinib arm"
  )
})

test_that("assess_plausibility errors for two-arm imputation", {
  cox_fit <- coxph(Surv(AVAL, EVENT) ~ TRT01P, ties = "exact", data = codebreak200)
  tipse_both <- suppressWarnings(tipping_point_model_based(
    dat = codebreak200, reason = "Early dropout",
    impute = c("docetaxel", "sotorasib"),
    J = 2, tipping_range = seq(0.5, 0.9, by = 0.2),
    cox_fit = cox_fit, seed = 123
  ))
  expect_error(assess_plausibility(tipse_both), "both arms")
})
