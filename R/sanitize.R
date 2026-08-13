sanitize_dataframe <- function(df) {
  # Check required columns
  required_cols <- c("AVAL", "TRT01P", "EVENT", "CNSRRS", "MAXAVAL")
  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  # Check AVAL is numeric
  if (!is.numeric(df$AVAL)) {
    stop("Column 'AVAL' must be numeric.")
  }

  # Check TRT01P is factor
  if (!is.factor(df$TRT01P)) {
    df$TRT01P <- as.factor(df$TRT01P)
    msg <- sprintf('Column "TRT01P" is not a factor, converting to factor with default level for control arm as "%s".', levels(df$TRT01P)[1])
    warning(msg, call. = FALSE)
  }

  # Check TRT01P is factor
  if (nlevels(df$TRT01P) != 2) {
    stop("Column 'TRT01P' must be a factor containing two levels", call. = FALSE)
  }

  # Check EVENT is binary (0/1)
  if (!all(df$EVENT %in% c(0, 1))) {
    stop("Column 'EVENT' must be binary (values 0 or 1).")
  }

  # Return sanitized dataframe
  return(df)
}

# Ensure percentile_range is a numeric vector
sanitize_percentile_range <- function(percentile_range) {
  if (!is.vector(percentile_range) || !is.numeric(percentile_range)) {
    stop("`tipping_range` must be a numeric vector.")
  }
  if (any(percentile_range <= 0 | percentile_range > 100)) {
    stop("All values in `tipping_range` must be between 0 and 100.")
  }
  return(percentile_range)
}

# Ensure npts_range is a numeric vector and does not exceed available patients
sanitize_npts_range <- function(dat, reason, impute, npts_range) {
  if (!is.vector(npts_range) || !is.numeric(npts_range)) {
    stop("`tipping_range` must be a numeric vector.")
  }
  n_available <- dat %>%
    subset(CNSRRS %in% reason & TRT01P == impute) %>%
    nrow()

  if (any(npts_range > n_available)) {
    stop(paste0("All values in `tipping_range` must be <= number of available patients (", n_available, ")."))
  }
  return(npts_range)
}

sanitize_tipse <- function(tipse) {
  if (!inherits(tipse, "tipse")) {
    stop("`tipse` must be an object of class 'tipse', returned by tipping_point_model_free() or tipping_point_model_based().")
  }

  # Extract basic elements
  if (!all(c("imputation_results", "imputation_data", "original_data", "arm_to_impute", "method_to_impute") %in% names(tipse))) {
    stop("The `tipse` object is missing required components: imputation_results, imputation_data, original_data, arm_to_impute, or method_to_impute.")
  }

  method <- tipse$method_to_impute
  valid_methods <- c("percentile sampling", "landmark sampling", "hazard multiplication")
  if (!method %in% valid_methods) {
    stop("Invalid method in `tipse$method_to_impute`. Must be one of: ", paste(valid_methods, collapse = ", "))
  }
}

# Ensure hazard_deflation_range is a numeric vector
sanitize_tipping_range <- function(range) {
  if (!is.vector(range) || !is.numeric(range)) {
    stop("`tipping_range` must be a numeric vector.")
  }
  if (any(range < 0)) {
    stop("All values in `tipping_range` must be greater than 0, i.e., deflation cannot be greater than 100%.")
  }
  return(range)
}
