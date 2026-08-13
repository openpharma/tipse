library(testthat)
library(survival)

# Create multiple imputations (list of data.frames)
dat_list <- list(
  codebreak200,
  codebreak200[sample(1:nrow(codebreak200), nrow(codebreak200), replace = TRUE), ],
  codebreak200[sample(1:nrow(codebreak200), nrow(codebreak200), replace = TRUE), ]
)

# Fit initial Cox model
cox_fit <- coxph(Surv(AVAL, EVENT) ~ TRT01P, ties = "exact", data = codebreak200)

test_that("pool_results returns correct structure and values", {
  res <- pool_results(dat_list, cox_fit, conf.level = 0.95)

  expect_s3_class(res, "data.frame")
  expect_true(all(c("HR", "log_HR_var", "HR_upperCI", "HR_lowerCI") %in% names(res)))

  # HR should be positive
  expect_true(res$HR > 0)
  # Variance should be positive
  expect_true(res$log_HR_var > 0)
  # CI bounds should make sense
  expect_true(res$HR_upperCI > res$HR_lowerCI)
})

test_that("pool_results handles different confidence levels", {
  res_90 <- pool_results(dat_list, cox_fit, conf.level = 0.90)
  res_95 <- pool_results(dat_list, cox_fit, conf.level = 0.95)

  # Wider CI for 95% than 90%
  expect_true((res_95$HR_upperCI - res_95$HR_lowerCI) > (res_90$HR_upperCI - res_90$HR_lowerCI))
})

test_that("pool_results works with single imputation (J=1)", {
  res <- pool_results(list(codebreak200), cox_fit, conf.level = 0.95)

  expect_s3_class(res, "data.frame")
  expect_true(res$HR > 0)
})

test_that("pool_results computes pooled HR correctly (manual check)", {
  # Compute manually for comparison
  results <- lapply(dat_list, function(df) {
    fit <- coxph(Surv(AVAL, EVENT) ~ TRT01P, ties = "exact", data = df)
    c(log_hr_est = fit$coefficients[1], log_hr_var = fit$var)
  })
  results_df <- do.call(rbind, results)

  hr_manual <- exp(mean(results_df[, 1]))
  res <- pool_results(dat_list, cox_fit)

  expect_equal(round(res$HR, 6), round(hr_manual, 6))
})

test_that("pool_results handles edge case with very small conf.level", {
  res <- pool_results(dat_list, cox_fit, conf.level = 0.01)
  expect_true(res$HR_upperCI > res$HR_lowerCI)
})
