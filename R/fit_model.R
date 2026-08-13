#' Fit parametric model for selected subjects
#'
#' @param dat data.frame containing at least 5 columns: TRT01P (treatment arm as factor), AVAL (survival time), EVENT (event indicator), CNSRRS (censoring reason) and MAXAVAL (maximum potential survival time, duration between randomization to data cut-off)
#' @param reason a string specifying the censoring reasons which require imputation. It must be one of the reasons from column CNSRRS.
#' @param impute a string specifying the treatment arm(s) which require imputation. It must be one of the arms from column TRT01P.
#' @param imputation_model a string specifying the parametric distribution used for imputation, can be "Weibull" or "exponential".
#'
#' @return data.frame with flags and fitted model parameters to be used for imputation
#' @importFrom survival survreg
#' @import dplyr
#' @importFrom MASS mvrnorm
#'
#' @details
#' The data.frame contains original columns, plus the following columns appended:
#' \tabular{ll}{
#' AVAL4  \tab Placeholder column to keep imputed survival times \cr
#' EVENT4 \tab Placeholder column to keep imputed events \cr
#' impute \tab Flag indicating whether the subject was selected for imputation \cr
#' a      \tab Shape parameter, equal to 1 if exponential \cr
#' b      \tab Scale parameter \cr
#' cdf    \tab Cumulative distribution function \cr
#' ...    \tab Some temporary columns
#' }
#' @keywords internal
fit_model <- function(dat, reason, impute, imputation_model = c("weibull", "exponential")) {
  # sanitize input data
  dat <- sanitize_dataframe(dat)
  imputation_model <- match.arg(imputation_model)
  # obtain treatment arm label
  control <- levels(dat[["TRT01P"]])[1]
  trt <- levels(dat[["TRT01P"]])[2]

  # check impute is one of the arms
  if (length(impute) > 1) {
    stop("Imputation only in one arm is allowed", call. = FALSE)
  }

  if (!(impute %in% c(control, trt))) {
    stop("Argument 'impute' must be one of the arms provided in column TRT01P.")
  }

  # Identify Patients from two arms who were not treated due or discontinued treatment (at any time)

  ids <- which(dat$TRT01P == impute & dat$CNSRRS %in% reason) # pts in control who discontinued
  if (length(ids) == 0) {
    stop("No subject identified meeting imputation criteria.")
  }

  ###  Weibull parameters estimated based on control arm patients excluding ids subset)
  fit_control <- survreg(Surv(AVAL, EVENT) ~ 1, data = dat %>% filter(TRT01P == control), dist = imputation_model)

  ###  Weibull parameters estimated based on trt arm patients excluding ids subset)
  fit_trt <- survreg(Surv(AVAL, EVENT) ~ 1, data = dat %>% filter(TRT01P == trt), dist = imputation_model)


  ab_control <- MASS::mvrnorm(nrow(dat), fit_control$icoef, fit_control$var)
  ab_trt <- MASS::mvrnorm(nrow(dat), fit_trt$icoef, fit_trt$var)

  if (imputation_model == "exponential") {
    ab_control <- cbind(ab_control, 0)
    ab_trt <- cbind(ab_trt, 0)
  }


  a_control <- 1 / exp(ab_control[, 2]) # shape = 1 / exp(log(scale))
  b_control <- exp(ab_control[, 1]) # scale = exp(intercept)

  a_trt <- 1 / exp(ab_trt[, 2]) # shape
  b_trt <- exp(ab_trt[, 1]) # scale


  # In patients to be imputed, the cumulative probability of event for these patient
  dat <- dat %>%
    mutate(AVALo = AVAL, EVENTo = EVENT) %>% # creating a copy of event and censor that will be replaced by imputed values where applicable
    mutate(impute = row_number() %in% ids) %>%
    bind_cols(a_control = a_control, b_control = b_control, a_trt = a_trt, b_trt = b_trt) %>%
    mutate(a = ifelse(TRT01P == control, a_control, a_trt), b = ifelse(TRT01P == control, b_control, b_trt)) %>%
    mutate(cdf = ifelse(impute, 1 - exp(-(AVAL / b)^a), NA))

  return(dat)
}
