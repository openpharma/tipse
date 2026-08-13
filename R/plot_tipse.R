#' Plot Pooled Kaplan–Meier Curves from Tipping Point Analysis
#'
#' Visualizes averaged (pooled) Kaplan-Meier survival curves across multiple tipping
#' point parameters, highlighting the tipping point where the upper CL of the hazard ratio crosses 1.
#'
#' @param x An S3 object of class `"tipse"` returned from \link{tipping_point_model_free} or \link{tipping_point_model_based}.
#' @param type Type of plot, either "Kaplan-Meier" or "Tipping Point".
#' @param ... Additional arguments to specify title, subtitle, xlab and ylab.
#'
#' @return A ggplot2 object displaying pooled Kaplan–Meier curves.
#' @import ggplot2
#' @export
#'
#' @details
#' - If `type = Kaplan-Meier`, then the KM curves from multiply imputed datasets were pooled using Rubin’s rules
#' after complementary log-log transformation as described in Marshall et al. (2009).
#' it can be of interest to visually assess the scenario that tips the result and the shift it causes to the original KM curve,
#' although there is no objective measure to assess the robustness of the result.
#' - If `type = Tipping Point`, then the HR estimation across the range of tipping point parameters are plotted.
#'
#' @references Marshall, A., Altman, D.G., Holder, R.L. et al.
#' Combining estimates of interest in prognostic modelling studies after multiple imputation: current practice and guidelines.
#' BMC Med Res Methodol 9, 57 (2009). https://doi.org/10.1186/1471-2288-9-57
#'
#' @examples
#' \donttest{
#' cox1 <- survival::coxph(Surv(AVAL, EVENT) ~ TRT01P, data = codebreak200)
#' result <- tipping_point_model_based(
#'   dat = codebreak200,
#'   reason = "Early dropout",
#'   impute = "docetaxel",
#'   imputation_model = "weibull",
#'   J = 10,
#'   tipping_range = seq(0.1, 1, by = 0.05),
#'   cox_fit = cox1
#' )
#'
#' plot(result, type = "Kaplan-Meier")
#' plot(result, type = "Tipping Point")
#'
#'
#' # Imputation in both arms
#' result2 <- tipping_point_model_based(
#'   dat = codebreak200,
#'   reason = "Early dropout",
#'   impute = c("docetaxel", "sotorasib"),
#'   imputation_model = "weibull",
#'   J = 10,
#'   tipping_range = list(seq(0.4, 0.7, by = 0.1), seq(0.5, 1.5, by = 0.2)),
#'   cox_fit = cox1
#' )
#'
#' plot(result2, type = "Tipping Point")
#' }
#'
plot.tipse <- function(x, type = c("Kaplan-Meier", "Tipping Point"), ...) {
  # Validate class
  sanitize_tipse(x)

  type <- match.arg(type, c("Kaplan-Meier", "Tipping Point"))

  arm <- x$arm_to_impute
  method <- x$method_to_impute

  if (type == "Kaplan-Meier") {
    # Inform the user
    if (length(arm) == 2) {stop("Plotting Kaplan-Meier curves is not supported for imputation done in both arms.")}
    message("Plotting pooled Kaplan-Meier curves for arm '", arm, "' using method '", method, "' approach.")
    p <- plot_km(x, ...)
  } else if (type == "Tipping Point") {
    message("Plotting tipping point results for arm '", paste(arm, collapse = " & "), "' using method '", method, "' approach.")
    if (length(arm) == 1){
      p <- plot_tp(x, ...)
  } else if (length(arm == 2)){
    p <- plot_tp_2d(x, ...)
  }
  }
  return(p)
}
