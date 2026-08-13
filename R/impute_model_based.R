#' Model-based imputation from parametric distributions
#'
#' @description
#' Impute data with Weibull or exponential distribution conditional on follow-up time
#'
#' @param dat data.frame containing at least 5 columns: TRT01P (treatment arm as factor), AVAL (survival time), EVENT (event indicator), CNSRRS (censoring reason) and MAXAVAL (maximum potential survival time, duration between randomization to data cut-off)
#' @param reason a string specifying the censoring reasons which require imputation. It must be one of the reasons from variable CNSRRS.
#' @param impute a string specifying the treatment arm(s) which require imputation. It must be one of the arms from variable TRT01P, the first level of TRT01P is considered as the control arm.
#' @param imputation_model a string specifying the parametric distribution used for imputation, can be "Weibull" or "exponential".
#' @param alpha hazard inflation (if treatment arm is imputed) or deflation (if control arm is imputed) rate
#' @param J numeric indicating number of imputations.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @details
#' First fit model based on the data without dropout. And then impute the survival outcome based on exponential or Weibull distribution for those who dropped out.
#'
#' @importFrom stats runif
#'
#' @return a list of data.frame from each imputation with imputed AVAL and EVENT, where original variables are kept as AVALo and EVENTo.
#' @export
#'
#' @examples
#' impute_model(
#'   dat               = codebreak200,
#'   reason            = "Early dropout",
#'   impute            = "docetaxel",
#'   imputation_model  = "weibull",
#'   alpha             = 0.7,
#'   J                 = 5,
#'   seed              = 1
#' )


impute_model <- function(dat, reason, impute, imputation_model = c("weibull", "exponential"), alpha, J, seed = NULL) {

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

  # sanitize input data
  dat <- sanitize_dataframe(dat)
  imputation_model <- match.arg(imputation_model)
  # obtain treatment arm label
  control <- levels(dat[["TRT01P"]])[1]
  trt <- levels(dat[["TRT01P"]])[2]

  # check impute in one of the arms
  if (length(impute) > 1) {
    stop("Imputation only in one arm is allowed", call. = FALSE)
  }

  if (!(impute %in% c(control, trt))) {
    stop("Argument 'impute' must be one of the arms provided in column TRT01P.")
  }

  ######## creating a copy of event and censor that will be replaced by imputed values where applicable ######
  dat <- fit_model(dat, reason = reason, impute = impute, imputation_model = imputation_model)
  dat_imp <- dat %>%
    filter(impute) %>%
    mutate(alpha = alpha)

  km_data <- list()

  dat_imp <- dat_imp %>%
    slice(rep(1:n(), times = J)) %>%
    mutate(rep = rep(1:J, each = nrow(dat_imp))) %>% # replicate J times
    mutate(U = runif(n(), min = cdf, max = 1))
  dat_imp <- dat_imp %>%
    mutate(t = b * ((-log(1 - U) - (1 - alpha) * (AVAL / b)^a) / alpha)^(1 / a)) %>%
    mutate(AVAL = ifelse(t <= MAXAVAL, t, MAXAVAL)) %>% # imputed time won't exceed data cutoff
    mutate(EVENT = (t <= MAXAVAL)) # EVENT: 1 - event and 0 - censored

  for (j in 1:J) {
    km_data[[j]] <- rbind(dat %>% filter(!impute) %>% select(-c(a, b, cdf)), dat_imp %>% filter(rep == j) %>% select(-c(a, b, cdf, alpha, rep, U, t))) %>%
      dplyr::mutate(iter = j, alpha = alpha)
  }


  return(km_data = km_data)
}
