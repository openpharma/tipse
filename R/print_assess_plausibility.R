#' Print method for plausibility assessment
#'
#' @param x An object of class "plausibility_assessment"
#' @param ... Further arguments passed to or from other methods.
#'
#' @keywords internal
#'
print.plausibility_assessment <- function(x, ...) {
  cat("\u2192 Clinical plausibility assessment (Oodally et al, 2025):\n\n")
  cat("  At the tipping point, ", x$assessment, "\n\n", sep = "")
  cat("  Please carefully assess the clinical plausibility in light of the imputation method and study context.\n")
  invisible(x)
}
