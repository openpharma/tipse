#' Plot Tipping Point Analysis Results
#'
#' Visualizes hazard ratios and corresponding confidence intervals across tipping point
#' parameters from tipping point analyses (one-dimensional). Supports both model-free and model-based
#' approaches, including percentile sampling, landmark sampling, and hazard-based methods.
#' The tipping point—defined as the parameter value where the upper confidence interval
#' crosses 1—is highlighted on the plot.
#'
#' @param tipse An S3 object of class `"tipse"` returned from
#'   \code{tipping_point_model_free()} or \code{tipping_point_model_based()}.
#' @param title Optional character string specifying the plot title. Set to \code{NULL}
#'   to remove the title.
#' @param subtitle Optional character string specifying the plot subtitle. If \code{NULL},
#'   a default subtitle including the method and tipping point is used.
#' @param xlab Optional character string specifying the x-axis label. If \code{NULL},
#'   a method-specific default is used.
#' @param ylab Optional character string specifying the y-axis label. Defaults to
#'   \code{"Hazard Ratio"}.
#'
#' @return A \code{ggplot2} object showing hazard ratios and confidence intervals across
#'   tipping point parameters. The returned plot can be further modified using \code{+}.
#'
#' @keywords internal
#'
plot_tp <- function(tipse,
                    title = "Tipping Point Analysis",
                    subtitle = NULL,
                    xlab = NULL,
                    ylab = "Hazard Ratio") {
  #----------------------------#
  # Input validation
  #----------------------------#
  summary_results <- tipse$imputation_results
  dat <- tipse$original_data
  arm <- tipse$arm_to_impute
  method <- tipse$method_to_impute

  control <- levels(dat[["TRT01P"]])[1]
  trt <- levels(dat[["TRT01P"]])[2]
  param_cols <- paste0("parameter_", arm)

  # Determine default x-axis label
  x_label <- switch(method,
                    "percentile sampling" = if (arm == control) {
                      "% best event times sampled from"
                    } else {
                      "% worst event times sampled from"
                    },
                    "landmark sampling" = if (arm == control) {
                      sprintf("Number of %s subjects event-free at DCO", arm)
                    } else {
                      sprintf("Number of %s subjects as events at censoring time", arm)
                    },
                    "hazard multiplication" = sprintf("Hazard multiplication factor (%%) in %s arm", arm)
  )

  if (method == "hazard multiplication") {
    summary_results[param_cols] <- summary_results[param_cols] * 100
  }

  #----------------------------#
  # Identify tipping point
  #----------------------------#
  tipping_point <- summary_results[summary_results$tipping_point, param_cols]
    scale_x <- scale_x_continuous(
    breaks = seq(min(summary_results[param_cols]), max(summary_results[param_cols]), by = 10),
    labels = function(x) round(x)
  )

  if (is.null(subtitle)) {
    subtitle <- paste("Method:", method,
                      "| Tipping point (HR upper 95% CL \u2265 1):",
                      tipping_point)
  }
  #----------------------------#
  # Build ggplot
  #----------------------------#
  p <- ggplot(summary_results, aes(x = .data[[param_cols]], y = HR)) +
    geom_point(size = 2) +
    geom_errorbar(
      aes(ymin = HR_lowerCI, ymax = HR_upperCI),
      width = 0.1 # controls horizontal width of the CI bar ends
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
    geom_vline(xintercept = tipping_point, linetype = "dashed", color = "red") +
    labs(
      x = if (is.null(xlab)) x_label else xlab,
      y = ylab,
      title = title,
      subtitle = subtitle
    ) +
    scale_x +
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

