#' Pooling results using Rubin's Rule
#'
#' @description
#' Pooling results from multiple imputations using Rubin's Rule
#'
#' @param dat a list of data.frames from multiple imputation using one alpha or kappa parameter
#' @param cox.fit a coxph object which is used to compute HRs for each imputed datasets
#' @param conf.level confidence level for the returned confidence interval, default to be 0.95.
#'
#' @details
#' The Rubin's rule is applied to the Cox PH model results across imputed datasets as:
#'
#' \enumerate{
#'   \item \emph{Compute pooled HR:}
#'   \deqn{\bar{HR}_\lambda = \exp\Bigg(\frac{1}{M} \sum_{m=1}^{M} \log(HR_m)\Bigg)}
#'
#'   \item \emph{Compute pooled variance:}
#'   \deqn{\bar{\sigma}_\lambda^2 = \frac{1}{M} \sum_{m=1}^{M} \sigma_m^2 +
#'   \frac{1 + \frac{1}{M}}{M-1} \sum_{m=1}^{M} \big(\log(HR_m) - \overline{\log(HR_\lambda)}\big)^2}
#'
#'   \item \emph{Compute CI:}
#'   \deqn{\bar{HR}_\lambda \times \exp\big(\pm t_{\alpha/2} \sqrt{\bar{\sigma}_\lambda^2}\big)}
#' }
#' @return a data.frame of pooled hazard ratio and confidence interval estimate using Rubin's Rule
#' @importFrom survival coxph
#' @importFrom stats var qt
#' @export
#'
#' @examples
#' cox_fit <- survival::coxph(Surv(AVAL, EVENT) ~ TRT01P, data = codebreak200)
#' imputed_list <- impute_percentile(
#'   dat        = codebreak200,
#'   reason     = "Early dropout",
#'   impute     = "docetaxel",
#'   percentile = 30,
#'   J          = 5,
#'   seed       = 1
#' )
#' pool_results(imputed_list, cox_fit)
pool_results <- function(dat, cox.fit, conf.level = 0.95) {
  # checks for parent function
  if (conf.level <= 0 | conf.level > 1) {
    stop("Argument 'conf.level' must be between 0 and 1.")
  }

  if (!inherits(cox.fit, "coxph")) {
    stop("Argument 'cox_fit' must be a valid cox model object, e.g. coxph(Surv(AVAL, EVENT) ~ TRT01P, data = codebreak200).")
  }


  J <- length(dat) # number of imputations

  small_cox <- function(dat.j) {
    model_call <- cox.fit$call # Extract original call
    model_call$data <- quote(dat.j) # Change the data argument
    fit.cox <- eval(model_call) # Re-run the model
    return(data.frame(log_hr_est = fit.cox$coefficients[1], log_hr_var = fit.cox$var))
  }

  results <- do.call(rbind, lapply(dat, small_cox))

  var_w <- mean(results$log_hr_var)
  var_b <- var(results$log_hr_est)
  var_t <- var_w + var_b + var_b/J

  lambda <- (1 + 1/J) * var_b/var_t
  q <- qt((1 - conf.level) / 2, df = (J - 1)/lambda^2, lower.tail = F)
  # Rubin's rule from J imputations
  hr_pool <- exp(mean(results$log_hr_est))
  log_hr_var_pool <- var_t # mean(results$log_hr_var) + (J + 1) / J * var(results$log_hr_est)
  hr_pool_upperCI <- exp(mean(results$log_hr_est) + q * sqrt(log_hr_var_pool))
  hr_pool_lowerCI <- exp(mean(results$log_hr_est) - q * sqrt(log_hr_var_pool))

  return(data.frame(HR = hr_pool, log_HR_var = log_hr_var_pool, HR_upperCI = hr_pool_upperCI, HR_lowerCI = hr_pool_lowerCI))
}
