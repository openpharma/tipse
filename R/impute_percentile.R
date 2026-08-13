#' Model-free imputation via percentile sampling
#'
#' @description
#' randomly sample from the percentile of best or worst patients (ordered by their observed times regardless of event or censoring) who do not require imputation.
#'
#' @param dat data.frame containing at least 5 columns: TRT01P (treatment arm as factor), AVAL (survival time), EVENT (event indicator), CNSRRS (censoring reason) and MAXAVAL (maximum potential survival time, duration between randomization to data cut-off)
#' @param reason a string specifying the censoring reasons which require imputation. It must be one of the reasons from variable CNSRRS.
#' @param impute a string specifying the treatment arm(s) which require imputation. It must be one of the arms from variable TRT01P, the first level of TRT01P is considered as the control arm.
#' @param percentile numeric between 1 and 100, indicating the best (or worst) percentile of subjects to sample from.
#' @param J numeric indicating number of imputations.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @details
#' We define two sets of subjects to sample from depending on the `impute` argument:
#'
#' 1. **Worst percentile of observations from treatment arm**
#'    \eqn{ \forall i \in N \mid \min\{T_i, C_i\} \leq F_{\min\{T_i, C_i\}}^{-1}(\kappa) }.
#' This set includes all indices \eqn{ i } where the minimum of \eqn{ T_i } (event time) and \eqn{ C_i } (censoring time) is **less than or equal to** the \eqn{\kappa}-th percentile of its distribution.
#'
#' 2. **Best percentile of observations control arm**
#'    \eqn{ \forall i \in N \mid \min\{T_i, C_i\} \geq F_{\min\{T_i, C_i\}}^{-1}(\kappa) }.
#' This set includes all indices \eqn{ i } where the minimum of \eqn{ T_i } and \eqn{ C_i } is **greater than or equal to** the \eqn{\kappa}-th percentile of its distribution.
#'
#' where \eqn{ F(\cdot) } denotes the cumulative distribution function (CDF) of the observed times and \eqn{F^{-1}(\kappa) } is the inverse CDF (quantile function) at percentile \eqn{\kappa}.
#'
#' @return a list of data.frame from each imputation with imputed AVAL and EVENT, where original variables are kept as AVALo and EVENTo.
#' @importFrom purrr map_dfr
#' @import dplyr
#' @export
#'
#' @examples
#' impute_percentile(
#'   dat        = codebreak200,
#'   reason     = "Early dropout",
#'   impute     = "docetaxel",
#'   percentile = 30,
#'   J          = 5,
#'   seed       = 1
#' )
impute_percentile <- function(dat, reason, impute, percentile, J, seed = NULL) {

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

  # obtain treatment arm label
  control <- levels(dat[["TRT01P"]])[1]
  trt <- levels(dat[["TRT01P"]])[2]

  # message(paste("Note:", control, "was chosen to be the control arm as the first level of variable TRT01P."))

  # check impute in one of the arms
  if (length(impute) > 1) {
    stop("Imputation only in one arm is allowed", call. = FALSE)
  }

  if (!(impute %in% c(control, trt))) {
    stop("Argument 'impute' must be one of the arms provided in column TRT01P.")
  }

  # check percentile greater than 0
  if (percentile <= 0 | percentile > 100) {
    stop("Argument 'percentile' must be greater than 0 and less or equal to 100.")
  }

  ids <- which(dat$TRT01P == impute & dat$CNSRRS %in% reason) # pts in control who discontinued
  if (length(ids) == 0) {
    stop("No subject identified meeting imputation criteria.")
  }

  ######## creating a copy of event and censor that will be replaced by imputed values where applicable ######
  dat <- dat %>% mutate(EVENTo = EVENT, AVALo = AVAL)
  dat_imp <- dat[ids, ]

  # define sampling population: best times if control arm or worst times if treatment arm
  impute_is_control <- identical(impute, control)
  dat_sample_select <- dat[-ids, ] %>%
    arrange(if (impute_is_control) {
      desc(AVAL)
    } else {
      AVAL
    }) %>%
    mutate(percentile_label = row_number() / n() * 100) %>%
    # select smallest percentile from data if parameter 'percentile' is too small
    filter(percentile_label <= max(min(percentile_label), percentile))

  if (nrow(dat_sample_select) < 1) {
    stop("There is no population to sample from, please check tipping range parameter values.")
  }

  km_data <- list()

  generate_replicate <- function(df, rep_id, n_to_impute = length(ids)) {
    df %>%
      sample_n(size = n_to_impute, replace = TRUE) %>%
      mutate(rep_sample = rep_id) # Add a replicate ID
  }

  # Generate J replicates
  dat_sample_select_rep <- purrr::map_dfr(1:J, ~ generate_replicate(dat_sample_select, .x)) %>%
    rename(AVAL_to_impute = AVAL, EVENT_to_impute = EVENT)

  dat_imp <- dat_imp %>%
    slice(rep(1:n(), times = J)) %>%
    mutate(rep = rep(1:J, each = nrow(dat_imp))) %>% # replicate J times
    bind_cols(dat_sample_select_rep %>% select(AVAL_to_impute, EVENT_to_impute, rep_sample)) %>%
    mutate(
      AVAL = if_else(AVAL_to_impute >= AVAL, AVAL_to_impute, AVAL), # use imputed time only if it's longer, since observed follow-up cannot change
      EVENT = EVENT_to_impute
    ) # use imputed event

  for (j in 1:J) {
    km_data[[j]] <- rbind(
      dat[-ids, ] %>% mutate(impute = FALSE),
      dat_imp %>% filter(rep == j) %>% select(-c(rep, AVAL_to_impute, EVENT_to_impute, rep_sample)) %>% mutate(impute = TRUE)
    ) %>%
      dplyr::mutate(iter = j, percentile = percentile)
  }

  return(km_data)
}
