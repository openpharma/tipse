#' Tipping Point Analysis (Model-Based)
#'
#' Performs a model-based tipping point analysis on time-to-event data by
#' repeatedly imputing censored observations under varying assumptions. The model-based
#' framework assumes that censored patients have a multiple of hazard fitted via a parametric survival model
#' compared to the rest of patients in the same arm (Akinson et al, 2019).
#'
#' @details
#' The **model-based tipping point analysis** provides a reproducible and intuitive
#' framework for exploring the robustness of treatment effects in time-to-event
#' (survival) endpoints when censoring may differ between study arms.
#'
#' A parametric survival model is fitted using maximum likelihood.
#' This function applies a hazard deflation on control arm or hazard inflation on treatment
#' arm, and impute survival times based on the parametric model with additional sampling of the parameters from a multivariate normal distribution.
#' This imputation procedure is iterated across a range of
#' tipping point parameters `tipping_range`. For each parameter value:
#' \enumerate{
#'   \item Multiple imputed datasets are generated (`J` replicates), where censored
#'         observations in the selected arm are reassigned
#'         event times according to the imputation method.
#'   \item A Cox proportional hazards model is fitted to each imputed dataset.
#'   \item Model estimates are pooled using **Rubin’s rules** to obtain a combined
#'         hazard ratio and confidence interval for that tipping point parameter.
#' }
#'
#' The process yields a series of results showing how the treatment effect changes
#' as increasingly conservative or optimistic assumptions are made about censored
#' observations. The *tipping point* is defined as the smallest value (hazard inflation)
#' or biggest value (hazard deflation) of the sensitivity parameter for which the upper
#' bound of the hazard ratio confidence interval crosses 1 - i.e., where the
#' apparent treatment benefit is lost.
#'
#' @param dat data.frame containing at least 5 columns: TRT01P (treatment arm as factor), AVAL (survival time), EVENT (event indicator), CNSRRS (censoring reason) and MAXAVAL (maximum potential survival time, duration between randomization to data cut-off)
#' @param reason Vector specifying censoring reasons to be imputed.
#' @param impute a character vector specifying the arm(s) to impute. Can be one arm or both arms (a length-2 vector). Each value must be one of the arms from variable TRT01P. When both arms are supplied, imputation is applied independently.
#' @param imputation_model used to fit model to observed data (should be "Weibull" or "exponential")
#' @param J numeric indicating number of imputations.
#' @param tipping_range Numeric vector, or when `length(impute) == 2` optionally a named or unnamed list of two numeric vectors (one per arm). When a list is supplied, all combinations of the two vectors are evaluated.
#' @param verbose Logical. If `TRUE`, prints progress and analysis details.
#' @param seed Integer, default as NULL. Random seed for reproducibility.
#' @param cox_fit A Cox model that will be used to calculate HRs on imputed datasets.
#'   In case of inclusion of stratification factors or covariates, conditional HR will be used.
#'
#' @return A `tipse` object containing:
#' \describe{
#'   \item{original data}{Input argument from 'data'.}
#'   \item{imputation_results}{A data frame of combined pooled model results across tipping points}
#'   \item{original_HR}{The original hazard ratio.}
#'   \item{reason_to_impute}{Input argument from 'reason'.}
#'   \item{arm_to_impute}{Input argument from 'impute'.}
#'   \item{method_to_impute}{Input argument from 'method'.}
#'   \item{imputation_data}{A list of imputed datasets for each tipping point value.}
#'   \item{seed}{Random seed.}
#' }
#' @import dplyr
#' @importFrom utils head tail
#' @export
#'
#' @references Atkinson, A., Kenward, M. G., Clayton, T., & Carpenter, J. R. (2019).
#' Reference‐based sensitivity analysis for time‐to‐event data. Pharmaceutical statistics, 18(6), 645-658.
#'
#' @examples
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
tipping_point_model_based <- function(dat,
                                      reason,
                                      impute,
                                      imputation_model = "weibull",
                                      J = 10,
                                      tipping_range = seq(0.05, 1, by = 0.05),
                                      cox_fit = NULL,
                                      verbose = FALSE,
                                      seed = NULL) {
  #----------------------------#
  # Setup and validation
  #----------------------------#
  # handle seed
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }

  dat <- sanitize_dataframe(dat)

  control <- levels(dat[["TRT01P"]])[1]
  trt <- levels(dat[["TRT01P"]])[2]

  if (!inherits(cox_fit, "coxph")) {
    stop("Argument 'cox_fit' must be a valid cox model object, e.g. coxph(Surv(AVAL, EVENT) ~ TRT01P, data = codebreak200).")
  }

  if (length(impute) > 2 || !all(impute %in% c(control, trt))) {
    stop("Argument 'impute' must be one or both of the arms in column TRT01P.")
  }
  if (length(impute) == 2 && impute[1] == impute[2]) {
    stop("Argument 'impute' must not contain the same arm twice.")
  }

  if (length(reason) == 0) {
    stop("Argument 'reason' must specify at least one censoring reason to impute.")
  }

  # Normalise tipping_range: accept list (per-arm) or plain vector (shared)
  two_arms <- length(impute) == 2
  if (is.list(tipping_range)) {
    if (!two_arms) stop("A list 'tipping_range' is only allowed when two arms are supplied in 'impute'.")
    if (length(tipping_range) != 2) stop("When 'tipping_range' is a list it must have exactly two elements, one per arm.")
    range1 <- sanitize_tipping_range(tipping_range[[1]])
    range2 <- sanitize_tipping_range(tipping_range[[2]])
  } else {
    range1 <- sanitize_tipping_range(tipping_range)
    range2 <- range1
  }

  # determine tipping point search direction from the first arm's range
  deflation <- all(range1 <= 1)
  method_to_impute <- "hazard multiplication"

  # column names for imputation results
  col1 <- paste0("parameter_", impute[1])
  col2 <- if (two_arms) paste0("parameter_", impute[2]) else NULL

  # Build the grid of parameter combinations to evaluate
  if (two_arms) {
    param_grid <- expand.grid(param1 = range1, param2 = range2)
  } else {
    param_grid <- data.frame(param1 = range1, param2 = NA_real_)
  }

  HR <- exp(cox_fit$coefficients[paste0("TRT01P", trt)])

  #----------------------------#
  # Print setup info
  #----------------------------#
  if (verbose) {
    cat("\u2192 Detected arms:\n")
    cat("   Control arm   :", control, "\n")
    cat("   Treatment arm :", trt, "\n")
    cat("   Imputing arm  :", paste(impute, collapse = " & "), "\n\n")
    cat("Starting tipping point analysis using method:", method_to_impute, "\n")
    cat("Replicates per tipping point parameter:", J, "\n\n")
  }

  #----------------------------#
  # Main computation function
  #----------------------------#
  run_imputation <- function(param1, param2, cox_fit, verbose, dat, reason, impute, imputation_model, J, seed) {
    if (verbose) {
      pct <- function(p) paste0(round(p * 100, 0), "%")
      label <- if (two_arms) paste(pct(param1), "&", pct(param2)) else pct(param1)
      cat(" \u2192 Imputing for:", method_to_impute, label, "\n")
    }

    multiply_imputed_dfs <- impute_model(dat, reason, impute[1], imputation_model, alpha = param1, J = J, seed = seed)
    if (two_arms) {
      # Drop model-fit columns added by the first imputation before passing to the second
      cols_to_drop <- c("a_control", "b_control", "a_trt", "b_trt", "a", "b", "cdf", "alpha", "iter", "impute")
      multiply_imputed_dfs <- lapply(multiply_imputed_dfs, function(d) {
        d <- d[, setdiff(names(d), cols_to_drop), drop = FALSE]
        impute_model(d, reason, impute[2], imputation_model, alpha = param2, J = 1, seed = seed)[[1]]
      })
    }

    pooled <- pool_results(multiply_imputed_dfs, cox_fit)
    pooled[[col1]] <- param1
    if (two_arms) pooled[[col2]] <- param2

    list(pooled = pooled, km_data = multiply_imputed_dfs)
  }

  #----------------------------#
  # Run all imputations
  #----------------------------#
  results <- mapply(
    run_imputation,
    param1 = param_grid$param1,
    param2 = param_grid$param2,
    MoreArgs = list(
      cox_fit = cox_fit, verbose = verbose,
      dat = dat, reason = reason, impute = impute, imputation_model = imputation_model,
      J = J, seed = seed
    ),
    SIMPLIFY = FALSE
  )

  summary_results <- dplyr::bind_rows(lapply(results, `[[`, "pooled"))
  km_data_list <- lapply(results, `[[`, "km_data")
  names(km_data_list) <- if (two_arms) {
    paste(param_grid$param1, param_grid$param2, sep = "|")
  } else {
    as.character(param_grid$param1)
  }

  #----------------------------#
  # Check tipping point
  #----------------------------#

  summary_results$tipping_point <- FALSE

  find_tip <- function(uci) {
    # returns the index within uci where the crossing occurs, or NA
    if (!any(uci < 1) || !any(uci > 1)) return(NA_integer_)
    idx <- if (deflation) tail(which(uci >= 1), 1) else head(which(uci >= 1), 1)
    if (abs(uci[idx] - 1) > 0.1)
      warning("The upper CL at tipping point was too far away from 1.0. Consider modifying the tipping_range.")
    idx
  }

  if (!two_arms) {
    tip <- find_tip(summary_results$HR_upperCI)
    if (!is.na(tip)) summary_results$tipping_point[tip] <- TRUE
  } else {
    # For each unique arm-2 parameter value (outer loop), find the tipping crossing
    # along the arm-1 axis (inner loop). Arm-1 varies within each arm-2 slice.
    for (p2 in unique(summary_results[[col2]])) {
      slice_idx <- which(summary_results[[col2]] == p2)
      tip_in_slice <- find_tip(summary_results$HR_upperCI[slice_idx])
      if (!is.na(tip_in_slice))
        summary_results$tipping_point[slice_idx[tip_in_slice]] <- TRUE
    }
  }


  if (verbose) cat("\nTipping point analysis completed successfully \u2705 \n")

  #----------------------------#
  # Return structured output
  #----------------------------#

  tipse <- list(
    original_data = dat,
    original_HR = HR,
    reason_to_impute = reason,
    arm_to_impute = impute,
    method_to_impute = method_to_impute,
    imputation_results = summary_results,
    imputation_data = km_data_list,
    seed = RNGstate
  )
  class(tipse) <- "tipse"
  return(tipse)
}
