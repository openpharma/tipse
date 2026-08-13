#' Patient level data from dummy trial
#'
#' Based on re-constructed Kaplan-Meier plot from CodeBreak 200 trial (de Langen et al., 2023)
#'
#' @format
#' A data frame with 345 rows and 5 columns:
#' \describe{
#'   \item{SUBJID}{Dummy patient ID}
#'   \item{TRT01P}{Treatment arm (Sotorasib or Docetaxel)}
#'   \item{AVAL}{PFS time in days}
#'   \item{EVENT}{Indicator for PFS event}
#'   \item{CNSRRS}{Censoring reason (Early dropout or Other)}
#'   \item{MAXAVAL}{Maximum potential survival time, duration between randomization to data cut-off}
#' }
#' @source De Langen, A.J., Johnson, M.L., Mazieres, J., Dingemans, A.M.C., Mountzios, G., Pless, M., Wolf, J., Schuler, M., Lena, H., Skoulidis, F. and Yoneshima, Y., 2023. Sotorasib versus docetaxel for previously treated non-small-cell lung cancer with KRASG12C mutation: a randomised, open-label, phase 3 trial. The Lancet, 401(10378), pp.733-746.
"codebreak200"
