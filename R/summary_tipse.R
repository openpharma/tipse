#' Summarize Tipping Point Results (ARD Format)
#'
#' Creates a concise, analysis-results dataset (ARD) from a tipping point
#' analysis. Identifies the tipping point parameter where the upper CL of the hazard ratio crosses 1
#' and summarizes key metrics.
#'
#' @param object A `tipse` object returned by \link{tipping_point_model_free} or
#'   \link{tipping_point_model_based}.
#' @param ... Additional arguments not used.
#'
#' @return A data frame summarizing:
#'   \itemize{
#'   \item \code{HR} - hazard ratio at that tipping point
#'   \item \code{CONFINT} - 95% CI at tipping point
#'     \item \code{METHOD} - sampling type used
#'     \item \code{ARMIMP} - arm imputed
#'     \item \code{TIPPT} - parameter where upper CL first crosses 1
#'     \item \code{TIPUNIT} - parameter meaning
#'     \item \code{DESC} - textual interpretation
#'   }
#' @export
#'
#' @examples
#' \donttest{
#' # Hazard deflation in the control arm
#' cox1 <- survival::coxph(Surv(AVAL, EVENT) ~ TRT01P, data = codebreak200)
#' result1 <- tipping_point_model_based(
#'   dat = codebreak200,
#'   reason = "Early dropout",
#'   impute = "docetaxel",
#'   imputation_model = "weibull",
#'   J = 10,
#'   tipping_range = seq(0.1, 1, by = 0.05),
#'   cox_fit = cox1
#' )
#'
#' summary(result1)
#'
#' # Imputation in both arms
#' result2 <- tipping_point_model_based(
#'   dat = codebreak200,
#'   reason = "Early dropout",
#'   impute = c("docetaxel", "sotorasib"),
#'   imputation_model = "weibull",
#'   J = 10,
#'   tipping_range = list(seq(0.1, 1, by = 0.2), seq(0.5, 1.5, by = 0.2)),
#'   cox_fit = cox1,
#'   verbose = TRUE,
#'   seed = 12345
#' )
#'
#' summary(result2)
#' }
summary.tipse <- function(object, ...) {
  #----------------------------#
  # Validation
  #----------------------------#
  sanitize_tipse(object)

  summary_results <- object$imputation_results
  km_data <- object$imputation_data
  dat <- object$original_data
  arm <- object$arm_to_impute
  method <- object$method_to_impute

  control <- levels(dat[["TRT01P"]])[1]
  trt <- levels(dat[["TRT01P"]])[2]

  param_cols <- paste0("parameter_", arm)
  #----------------------------#
  # Identify tipping point
  #----------------------------#
  tipping_param <- summary_results[summary_results$tipping_point,]

  if (method == "percentile sampling") {
    tipping_unit <- paste(ifelse(arm == control, "best", "worst"), "percentile")
  } else if (method == "landmark sampling") {
    tipping_unit <- paste("number of subjects", ifelse(arm == control, "extended censoring", "set as events"))
  } else if (method == "hazard multiplication") {
    tipping_unit <- "% hazard multiplication"
  } else {
    stop("Unsupported tipping point method: ", method)
  }

  if (nrow(tipping_param) == 0 || is.null(tipping_param)) {
    param_tipping <- NA
    hr_tipping <- NA
    ci_tipping <- NA
    interp <- "No tipping point detected (upper CL < 1 across all parameters)."
  } else {
    hr_tipping <- tipping_param$HR
    ci_tipping <- sprintf("(%.4f-%.4f)", tipping_param$HR_lowerCI, tipping_param$HR_upperCI)
    param_tipping <- tipping_param[,param_cols]
    if (method == "hazard multiplication") {
      if (is.null(dim(param_tipping))){
param_tipping <- param_tipping*100
      } else {
      param_tipping <- apply(param_tipping*100, 1, paste, collapse = ",")
      }
    }
    interp <- paste0(
      "Tipping point (upper CL \u2265 1) reached when imputing ", paste(arm, collapse = ","),
      " arm at ", param_tipping, " ", tipping_unit, "."
    )
  }

  #----------------------------#
  # Summarize table (ARD format)
  #----------------------------#
  summary_tbl <- data.frame(
    HR = round(hr_tipping, 4),
    CONFINT = ci_tipping,
    METHOD = method,
    ARMIMP = paste(arm, collapse = ","),
    TIPPT = param_tipping,
    TIPUNIT = tipping_unit,
    DESC = interp,
    stringsAsFactors = FALSE
  )

  return(summary_tbl)
}
