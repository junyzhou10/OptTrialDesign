#' @title Calculate the number of accrued samples (Not for all)
#' @description Function used for calculating number of cumulative accrued samples based on a piecewise accrual rate. Used for piecewise accrual case.
#' @param Sa Accrual duration
#' @param ra.seq Sequence of piecewise accrual rate. If length(ra.seq)>Sa, it directly assumes using ra.seq[-1] (the last one) as the accrual rate for the rest accrual period
#' @return A number of accrued samples
#' @examples cum.N(10, c(10,15,20,30,25))
#' @export

cum.N <- function(Sa, ra.seq) {
  if (length(ra.seq) > Sa) { # need to be truncated
    N = sum(ra.seq[1:floor(Sa)])+(Sa-floor(Sa))*ra.seq[ceiling(Sa)]
  } else { # need to extend
    N = sum(ra.seq) + (Sa - length(ra.seq))*ra.seq[length(ra.seq)]
  }
  return(N)
}
