library(testthat)
library(dplyr)

# Cox model for testing
cox_fit <- coxph(Surv(AVAL, EVENT) ~ TRT01P, ties = "exact", data = codebreak200)
# tipse object
results <- tipping_point_model_free(
  dat = codebreak200,
  reason = "Early dropout",
  impute = "docetaxel",
  J = 2,
  tipping_range = seq(10, 80, by = 10),
  cox_fit = cox_fit,
  method = "percentile sampling",
  seed = 123
)

test_that("summary returns correct structure for percentile sampling", {
  res <- summary(results)

  expect_s3_class(res, "data.frame")
  expect_true(all(c("HR", "CONFINT", "METHOD", "ARMIMP", "TIPPT", "TIPUNIT", "DESC") %in% names(res)))
  expect_equal(res$METHOD, "percentile sampling")
  expect_true(grepl("percentile", res$TIPUNIT))
})


test_that("summary returns correct structure for landmark sampling", {
  # tipse object
  results2 <- tipping_point_model_free(
    dat = codebreak200,
    reason = "Early dropout",
    impute = "docetaxel",
    J = 2,
    tipping_range = seq(1, 21, by = 5),
    cox_fit = cox_fit,
    method = "landmark sampling",
    seed = 123
  )
  res <- summary(results2)

  expect_equal(res$METHOD, "landmark sampling")
  expect_true(grepl("number of subjects", res$TIPUNIT))
})


test_that("summary handles hazard multiplication (inflation and deflation)", {
  results3 <- tipping_point_model_based(
    dat = codebreak200,
    reason = "Early dropout",
    impute = "sotorasib",
    J = 2,
    tipping_range = seq(1, 2, by = 0.05),
    cox_fit = cox_fit,
    seed = 123
  )
  res_infl <- summary(results3)
  expect_true(grepl("hazard multiplication", res_infl$TIPUNIT))

  results4 <- tipping_point_model_based(
    dat = codebreak200,
    reason = "Early dropout",
    impute = "docetaxel",
    J = 2,
    tipping_range = seq(0.1, 1, by = 0.05),
    cox_fit = cox_fit,
    seed = 123
  )
  res_defl <- summary(results4)
  expect_true(grepl("hazard multiplication", res_defl$TIPUNIT))
})

test_that("summary handles case with no tipping point detected", {
  expect_warning(
    results5 <- tipping_point_model_free(
      dat = codebreak200,
      reason = "Early dropout",
      impute = "sotorasib",
      J = 2,
      tipping_range = c(20, 40),
      cox_fit = cox_fit,
      method = "percentile sampling",
      seed = 123
    ),
    "Tipping point not found"
  )
  res <- summary(results5)

  expect_true(is.na(res$TIPPT))
  expect_true(grepl("No tipping point detected", res$DESC))
})

test_that("summary formats CI and interpretation correctly", {
  res <- summary(results)

  expect_true(grepl("\\(", res$CONFINT))
  expect_true(grepl("Tipping point", res$DESC))
})

test_that("summary two-arm: returns one row per arm", {
  res_both <- suppressWarnings(tipping_point_model_based(
    dat = codebreak200, reason = "Early dropout",
    impute = c("docetaxel", "sotorasib"),
    J = 2, tipping_range = seq(0.5, 0.9, by = 0.2),
    cox_fit = cox_fit, seed = 123
  ))
  res <- summary(res_both)

  expect_equal(nrow(res), 3L)
  expect_equal(unique(res$ARMIMP), c("docetaxel,sotorasib"))
  expect_true(all(c("HR", "CONFINT", "METHOD", "ARMIMP", "TIPPT", "TIPUNIT", "DESC") %in% names(res)))
  expect_true(all(res$METHOD == "hazard multiplication"))
})
