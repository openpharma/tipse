#' Model-free imputation via landmark sampling
#'
#' @description
#' patients will be assigned deterministically an event time at the time of censoring  or extend the censoring time to the potential maximum follow-up of each patient.
#'
#' @param dat data.frame containing at least 5 columns: TRT01P (treatment arm as factor), AVAL (survival time), EVENT (event indicator), CNSRRS (censoring reason) and MAXAVAL (maximum potential survival time, duration between randomization to data cut-off)
#' @param npts number of patients to be imputed
#' @param reason a string specifying the censoring reasons which require imputation. It must be one of the reasons from variable CNSRRS.
#' @param impute a string specifying the treatment arm(s) which require imputation. It must be one of the arms from variable TRT01P, the first level of TRT01P is considered as the control arm.
#' @param J numeric indicating number of imputations.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @import dplyr
#' @details
#' patients will be assigned deterministically an event time at the time of censoring  or extend the censoring time to the potential maximum follow-up of each patient.

#' @return a list of data.frame from each imputation with imputed AVAL and EVENT, where original variables are kept as AVAL and EVENT.
#' @export
#'
#' @examples
#' impute_landmark(
#'   dat    = codebreak200,
#'   reason = "Early dropout",
#'   impute = "docetaxel",
#'   npts   = 5,
#'   J      = 5,
#'   seed   = 1
#' )
impute_landmark <- function(dat, reason, impute, npts, J, seed) {

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

  ids <- which(dat$TRT01P == impute & dat$CNSRRS %in% reason) # pts in control who discontinued
  if (length(ids) == 0) {
    stop("No subject identified meeting imputation criteria.")
  }

  # check ids
  if (npts < 0 | npts > length(ids)) {
    stop("Argument 'npts' must be greater than 0 and less or equal to the total number of patients meeting imputation criteria.")
  }

  ######## creating a copy of event and censor that will be replaced by imputed values where applicable ######
  dat <- dat %>% mutate(EVENTo = EVENT, AVALo = AVAL)
  dat_imp <- dat[ids, ]

  km_data <- list()

  npts_val <- npts
  impute_is_control <- identical(impute, control)
  impute_is_trt     <- identical(impute, trt)
  dat_imp <- dat_imp %>%
    slice(rep(1:n(), times = J)) %>%
    mutate(rep = rep(1:J, each = nrow(dat_imp))) %>% # replicate J times
    group_by(rep) %>%
    mutate(imp = (1:n() %in% sample(n(), size = npts_val, replace = FALSE))) %>%
    ungroup() %>%
    mutate(AVAL = ifelse(imp & impute_is_control, MAXAVAL, AVAL)) %>%
    mutate(EVENT = ifelse(imp & impute_is_trt, 1, 0))

  for (j in 1:J) km_data[[j]] <- rbind(dat[-ids, ] %>% mutate(impute = FALSE), dat_imp %>% filter(rep == j) %>% select(-c(rep, imp)) %>% mutate(impute = TRUE)) %>% mutate(iter = j, npts = npts_val)


  return(km_data = km_data)
}
