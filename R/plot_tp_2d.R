#' Heatmap of Two-Dimensional Tipping Point Analysis
#'
#' Visualises the upper confidence limit of the hazard ratio across a Cartesian
#' grid of two imputation parameters (one per arm). Cells are shaded by
#' \code{HR_upperCL}, with the tipping line (where \code{HR_upperCI} crosses 1)
#' highlighted by a dashed contour line and flagged cells marked with a cross.
#'
#' @param tipse A \code{tipse} object returned by \link{tipping_point_model_based}
#'   or \link{tipping_point_model_free} with \code{impute} of length 2.
#' @param title Optional character string specifying the plot title. Set to \code{NULL}
#'   to remove the title.
#' @param subtitle Optional character string specifying the plot subtitle. If \code{NULL},
#'   a default subtitle including the method and tipping point is used.
#' @param xlab Optional character string specifying the x-axis label. If \code{NULL},
#'   a method-specific default is used.
#' @param ylab Optional character string specifying the y-axis label. Defaults to
#'   \code{"Hazard Ratio"}.
#' @param cross Optional logical to specifying if tipping points will be marked with a cross.
#' @param line Optional logical to specifying if tipping line will be marked with a contour.
#' @return A \code{ggplot2} object.
#' @keywords internal
#'
plot_tp_2d <- function(tipse,
                       title = "Two-Dimensional Tipping Point Analysis",
                       subtitle = NULL,
                       xlab = NULL,
                       ylab = NULL,
                       cross = TRUE,
                       line = TRUE) {
  #----------------------------#
  # Input validation
  #----------------------------#
  sanitize_tipse(tipse)

  arms <- tipse$arm_to_impute
  if (length(arms) != 2L) {
    stop("plot_tp_2d() requires a two-arm tipse object (length(impute) == 2).")
  }

  sr     <- tipse$imputation_results
  method <- tipse$method_to_impute

  col1 <- paste0("parameter_", arms[1])
  col2 <- paste0("parameter_", arms[2])

  #----------------------------#
  # Axis labels
  #----------------------------#
  dat     <- tipse$original_data
  control <- levels(dat[["TRT01P"]])[1]

  axis_label <- function(arm) {
    if (method == "percentile sampling") {
      paste0(ifelse(arm == control, "% best", "% worst"), " event times (", arm, ")")
    } else if (method == "landmark sampling") {
      paste0("No. subjects ", ifelse(arm == control, "extended censoring", "set as events"),
             " (", arm, ")")
    } else {
      # hazard multiplication: show % deviation from 1
      paste0("Hazard multiplication factor (", arm, ")")
    }
  }

  #----------------------------#
  # Scale parameters for display
  #----------------------------#
  plot_data <- sr

  x_label <- if (is.null(xlab)) axis_label(arms[1])
  y_label <- if (is.null(ylab)) axis_label(arms[2])

  if (is.null(subtitle) & line) {
    subtitle <- paste("Method:", method,
                      "| Dashed line marks tipping line (HR upper 95% CL \u2265 1)")
  } else if (is.null(subtitle)){
    subtitle <- paste("Method:", method)
  }

  #----------------------------#
  # Build heatmap
  #----------------------------#
  p <- ggplot(
    plot_data,
    aes(
      x    = .data[[col1]],
      y    = .data[[col2]]
    )
  ) +
    geom_tile(aes(fill = HR_upperCI), colour = "white", linewidth = 0.4) +
    # diverging palette centred on 1 (the decision boundary)
    scale_fill_gradient2(
      low      = "#2166ac",
      mid      = "#f7f7f7",
      high     = "#d6604d",
      midpoint = 1,
      name     = "HR Upper 95% CL"
    ) +
    labs(
      x        = x_label,
      y        = y_label,
      title    = title,
      subtitle = subtitle
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray30"),
      plot.caption = element_text(hjust = 0),
      legend.position = "bottom",
      legend.title = element_text(hjust = 0.5)
    )

  # mark tipping-line cells
  if (cross) {p <- p + geom_point(
    data  = plot_data[plot_data$tipping_point, , drop = FALSE],
    shape = 4,   # cross
    size  = 3,
    colour = "black",
    stroke = 1.2
  )}
    # HR_upperCI = 1 contour
    if (line) {p <- p + geom_contour(
      aes(z = HR_upperCI),
      breaks   = 1,
      colour   = "black",
      linetype = "dashed",
      linewidth = 0.7
    )}

  return(p)
}
