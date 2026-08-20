#' Assess Clinical Plausibility of Imputation Results
#'
#' This function facilitates the evaluation of clinical plausibility at the tipping point. It provides a text summary comparing event rates,
#' follow-up duration, or hazard ratios between treatment arms depending on the
#' imputation method and arm specified. NOTE: this function only supports imputation in one arm.
#'
#' @param tipse A `tipse` object returned by one of \link{tipping_point_model_free} or \link{tipping_point_model_based}.
#' @param verbose Logical. If `TRUE`, prints assessment details.
#'
#' @return A character string summarizing the key information to facilitate clinical plausibility assessment based on
#' the imputation scenario.
#'
#' @import dplyr
#' @importFrom stats median
#' @importFrom survival coxph Surv
#' @export
#'
#' @examples
#' cox1 <- survival::coxph(survival::Surv(AVAL, EVENT) ~ TRT01P, data = codebreak200)
#' result <- tipping_point_model_free(
#'   dat = codebreak200,
#'   reason = "Early dropout",
#'   impute = "docetaxel",
#'   cox_fit = cox1,
#'   method = "percentile sampling"
#' )
#'
#' assess_plausibility(result)
#'
assess_plausibility <- function(tipse, verbose = TRUE) {
  sanitize_tipse(tipse)

  if (length(tipse$arm_to_impute) == 2) {
    stop("The tipping point analysis involved imputation in both arms, it is not supported by this function.")
  }
  summary_results <- tipse$imputation_results
  km_data <- tipse$imputation_data

  control <- levels(tipse$original_data[["TRT01P"]])[1]
  trt <- levels(tipse$original_data[["TRT01P"]])[2]

  if (!any(summary_results$tipping_point)) {
    stop("No tipping point found, cannot assess plausibility.")
  }

  # cat("\u2192 Clinical plausibility assessment (Oodally et al, 2025): \n")
  param_col_name <- paste0("parameter_", tipse$arm_to_impute)
  tipping_param_val <- summary_results %>% filter(tipping_point) %>% pull(param_col_name)
  imputed_data <- km_data[[as.character(tipping_param_val)]]

  assessment <- NULL
  metric <- NULL

  if (tipse$arm_to_impute == trt) { # imputation in treatment arm

    if (tipse$method %in% c("percentile sampling", "landmark sampling")) {
      # compare event rate between imputed treatment and control
      rate_control <- tipse$original_data %>%
        filter(TRT01P == control) %>%
        summarise(rate = mean(EVENT)) %>%
        pull(rate)
      rate_trt <- mean(unlist(lapply(imputed_data, function(x) {
        x %>%
          filter(impute) %>%
          summarise(rate = mean(EVENT)) %>%
          pull(rate)
      })))

      assessment <- paste0(
        "the event rate in imputed set in ", trt, " arm was ", round(rate_trt * 100, 2),
        ", compared to ", round(rate_control * 100, 2), " in the ", control, " arm."
      )

      metric <- list(rate_trt = rate_trt, rate_control = rate_control)
    } else {
      # compare HR between imputed treatment vs control

      HR_imputed_vs_nonImputed <- 1 + tipping_param_val # 1 + inflation

      newHR <- HR_imputed_vs_nonImputed * tipse$original_HR

      assessment <- paste0("the HR between imputed set in ", trt, " arm and ", control, " arm was approximately ", newHR, ".")

      metric <- list(HR = newHR)
    }
  } else { # imputation in control arm

    if (tipse$method %in% c("percentile sampling", "landmark sampling")) {
      # compare follow up between imputed control vs treatment

      fup_trt <- tipse$original_data %>%
        filter(TRT01P == trt) %>%
        summarise(fup = median(AVAL)) %>%
        pull(fup)
      fup_control <- mean(unlist(lapply(imputed_data, function(x) {
        x %>%
          filter(impute) %>%
          summarise(fup = median(AVAL)) %>%
          pull(fup)
      })))

      assessment <- paste0(
        "the median duration of follow-up in imputed set in ", control, " arm was ", round(fup_control, 2),
        ", compared to ", round(fup_trt, 2), " in the ", trt, " arm."
      )

      metric <- list(fup_control = fup_control, fup_trt = fup_trt)
    } else {
      HR_imputed_vs_nonImputed <- 1 - tipping_param_val # 1 - deflation

      newHR <- HR_imputed_vs_nonImputed / tipse$original_HR

      # compare HR between imputed control vs treatment
      assessment <- paste0("the HR between imputed set in ", control, " arm and ", trt, " arm was approximately ", round(newHR, 2), ".")

      metric <- list(HR = newHR)
    }
  }
  assessment <- list(
      assessment = assessment,
      metric = metric,
      arm_imputed = tipse$arm_to_impute,
      method = tipse$method
    )
  class(assessment) <- "plausibility_assessment"
  if (verbose == FALSE){
    return(assessment)
  } else {
    print(assessment)
  }
  # cat(paste0("  At the tipping point, ", assessment, "\n", "  Please carefully assess the clinical plausibility based on this comparison together with the imputation method used in your study context."))
}
