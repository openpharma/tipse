#' Plot Pooled Kaplan–Meier Curves from Model-Free Tipping Point Analysis
#'
#' Visualizes averaged (pooled) Kaplan-Meier survival curves across multiple tipping
#' point parameters, highlighting the tipping point where the upper CL of the hazard ratio crosses 1.
#'
#' @param tipse An S3 object of class `"tipse"` returned from \code{tipping_point_model_free()} or \code{tipping_point_model_based()}.
#'
#' @importFrom survival survfit Surv
#' @return A \code{ggplot2} object displaying pooled Kaplan-Meier curves, with:
#' \itemize{
#'   \item Colored lines - pooled KM curves for each tipping point parameter
#'   \item Red line - tipping point where HR upper CL crosses 1
#'   \item Orange dashed line - original (unimputed) KM curve
#' }
#' @keywords internal

plot_km <- function(tipse) {
  #----------------------------#
  # Sanitize and validate input
  #----------------------------#

  sanitize_tipse(tipse)
  summary_results <- tipse$imputation_results
  dat <- tipse$original_data
  arm <- tipse$arm_to_impute
  method <- tipse$method_to_impute
  km_data <- tipse$imputation_data

  # Identify control and treatment arms
  control <- levels(dat$TRT01P)[1]
  trt <- levels(dat$TRT01P)[2]
  param_cols <- paste0("parameter_", arm)

  #----------------------------#
  # Average KM curves
  #----------------------------#
  res_all_km <- dplyr::bind_rows(
    lapply(
      km_data,
      FUN = average_km,
      arm = arm
    )
  )

  #----------------------------#
  # Identify tipping point
  #----------------------------#

  key <- paste(arm, method, sep = "_")

  highlight_param <- summary_results[summary_results$tipping_point, param_cols]

  if (is.null(highlight_param)) {
    warning("Could not determine tipping point.")
  }

  if (key == paste(control, "hazard multiplication", sep = "_")) {
    res_all_km$parameter <- res_all_km$parameter * 100
    highlight_param <- highlight_param * 100
  }

  highlight_data <- res_all_km[res_all_km$parameter == highlight_param, , drop = FALSE]

  #----------------------------#
  # Original (unimputed) KM curve
  #----------------------------#
  survfit_original <- survfit(Surv(AVAL, EVENT) ~ 1, data = subset(dat, TRT01P == arm))
  survfit_original_data <- data.frame(
    time = survfit_original$time,
    surv = survfit_original$surv,
    parameter = 0
  )

  #----------------------------#
  # Plot
  #----------------------------#
  p <- ggplot(res_all_km, aes(x = time, y = survival_comb, color = parameter, group = parameter)) +
    geom_step(alpha = 0.7) +
    scale_color_gradient(
      low = "darkblue", high = "lightblue",
      name = switch(method,
        "percentile sampling" = ifelse(arm == control, paste("% best event times sampled from"),
          paste("% worst event times sampled from")
        ),
        "landmark sampling" = ifelse(arm == control, paste("Number of ", arm, "subjects event-free at DCO"),
          paste("Number of", arm, "subjects as events at censoring time")
        ),
        "hazard multiplication" = paste("Hazard multiplication factor (%) in", arm, "arm"),
      ),
      guide = guide_colorbar(title.position = "bottom", title.hjust = 0.5)
    ) +
    geom_step(
      data = highlight_data, aes(x = time, y = survival_comb),
      color = "red", linewidth = 1.2, linetype = "solid"
    ) +
    geom_step(
      data = survfit_original_data, aes(x = time, y = surv),
      color = "orange", linewidth = 1, linetype = "dashed"
    ) +
    labs(
      title = paste("Pooled Kaplan-Meier Curves:", arm, "arm"),
      subtitle = paste("Method:", method, "| Tipping point (HR upper 95% CL \u2265 1):", highlight_param),
      caption = paste("Original KM curve is shown in orange and KM curve at tipping point is shown in red."),
      x = "Time",
      y = "Survival Probability"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray30"),
      plot.caption = element_text(hjust = 0),
      legend.position = "bottom",
      legend.title = element_text(hjust = 0.5)
    )

  return(p)
}
