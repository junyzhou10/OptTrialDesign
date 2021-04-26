#' @title Calculate the time of accrual (Not for all)
#' @description Function used for calculating the time of accrual time based on the given number of accrued samples. Used for piecewise accrual case.
#' @param N Number of total samples. If sum(ra.seq)<N, it directly assumes using ra.seq[-1] (the last one) as the accrual rate for the rest samples
#' @param ra.seq Sequence of piecewise accrual rate.
#' @param Acc Number of already accrued samples
#' @return Time of accural period

inv.N <- function(N, ra.seq, Acc){ # knowing Acc(acc of enrolled sample), return time
  if (sum(ra.seq) > N) {
    id = which(cumsum(ra.seq)<N)
    ra.seq = ra.seq[c(id, max(id)+1)]
  } else {
    r.last = ra.seq[length(ra.seq)]
    ra.seq = c(ra.seq, rep(r.last, ceiling((N-sum(ra.seq))/r.last)))
  }

  id = which(cumsum(ra.seq)>=Acc)[1]
  inv.T = id - 1 + (Acc - ifelse(id==1,0,sum(ra.seq[1:(id-1)])))/ra.seq[id]
  return(inv.T)
}
