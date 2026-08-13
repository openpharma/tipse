#' Average Kaplan-Meier Curves Across Multiple Imputed Datasets
#'
#' Takes a list of multiply imputed datasets corresponding to a single tipping point
#' parameter and pools the Kaplan-Meier survival curves for a given treatment arm.
#' Uses log(-log) transformation and Rubin's rules for pooling across imputations.
#'
#' @param km_data List of data frames, each containing one multiply imputed dataset for a tipping point.
#'   Each data frame must contain columns `AVAL`, `EVENT`, `TRT01P`, `iter`, and a tipping point parameter depending on the method.
#' @param arm Character string specifying the treatment arm to pool (must match `TRT01P` levels).
#' @param conf_level Numeric. Confidence level for CIs (default = 0.95).
#'
#' @importFrom stats var qnorm
#' @return A data frame with the following columns:
#' \describe{
#'   \item{time}{Time points of the KM curve.}
#'   \item{survival_comb}{Pooled survival probability at each time point.}
#'   \item{survival_lcl_comb}{Lower 95% confidence limit of pooled survival.}
#'   \item{survival_ucl_comb}{Upper 95% confidence limit of pooled survival.}
#'   \item{stderr}{Standard error of the pooled log-log transformed estimate.}
#' }
#' @keywords internal
average_km <- function(km_data, arm, conf_level = 0.95) {
  if (!is.list(km_data) || length(km_data) == 0) stop("'km_data' must be a non-empty list")

  km_surv_probs <- vector("list", length(km_data))

  par <- names(km_data[[1]])[ncol(km_data[[1]])] # get parameter name

  for (i in seq_along(km_data)) {
    df <- km_data[[i]]
    if (!all(c("AVAL", "EVENT", "TRT01P", "iter") %in% names(df))) {
      stop("Each km_data element must contain columns: AVAL, EVENT, TRT01P, iter")
    }

    # subset to the specified arm and iteration
    arm_data <- df[df$TRT01P == arm & df$iter == i, ]

    # fit Kaplan-Meier model
    km_fit <- survfit(Surv(AVAL, EVENT) ~ 1, data = arm_data)
    km_sum <- summary(km_fit)

    # log(-log) transformation
    valid_idx <- km_sum$surv > 0 & km_sum$surv < 1
    km_surv_probs[[i]] <- data.frame(
      time = km_sum$time[valid_idx],
      survival_cll = log(-log(km_sum$surv[valid_idx])),
      stderr_cll = sqrt((1 / log(km_sum$surv[valid_idx])^2) *
        ((km_sum$std.err[valid_idx]^2) / (km_sum$surv[valid_idx]^2))),
      imputation = i,
      parameter = unique(arm_data[[par]])
    )
  }

  combined_data <- do.call(rbind, km_surv_probs)

  # keep time points present in more than one imputation
  time_counts <- table(combined_data$time)
  valid_times <- as.numeric(names(time_counts[time_counts > 1]))
  combined_data <- combined_data[combined_data$time %in% valid_times, ]

  # pool estimates across imputations using Rubin's rules
  pooled_results <- do.call(rbind, lapply(split(combined_data, combined_data$time), function(group) {
    m <- length(unique(group$imputation))
    W <- mean(group$stderr_cll^2)
    B <- var(group$survival_cll)
    T_var <- W + (1 + 1 / m) * B
    pooled_estimate <- mean(group$survival_cll)
    pooled_stderr <- sqrt(T_var)
    data.frame(
      time = unique(group$time),
      estimate = pooled_estimate,
      stderr = pooled_stderr,
      parameter = unique(group$parameter), row.names = F
    )
  }))

  # back-transform to survival scale and compute confidence intervals
  zval <- qnorm(1 - (1 - conf_level) / 2)
  pooled_results$survival_comb <- exp(-exp(pooled_results$estimate))
  pooled_results$survival_lcl_comb <- exp(-exp(pooled_results$estimate + zval * pooled_results$stderr))
  pooled_results$survival_ucl_comb <- exp(-exp(pooled_results$estimate - zval * pooled_results$stderr))

  return(pooled_results[, c("time", "parameter", "survival_comb", "survival_lcl_comb", "survival_ucl_comb", "stderr")])
}
