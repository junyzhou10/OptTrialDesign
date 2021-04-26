#' @title Calculate the expected number of events (single arm, piecewise accrual rate)
#' @description A very fundamental function to calculate the expected number of events given sample size (N), study duration (S), and either one of accrual rate (ra) or accural duration (Sa)
#' @param S study duration, or equivalently, accrual time (Sa) + follow-up time (Sf)
#' @param r.seq Sequence of piecewise accrual rate. If length(ra.seq)>Sa, it directly assumes using ra.seq[-1] (the last one) as the accrual rate for the rest accrual period
#' @param Sa accrual duration, a numeric value
#' @param lambda rate parameter for event time (assuming an exponential distribution), a numeric value > 0
#' @param eta rate parameter for loss to follow-up (assuming an exponential distribution), a numeric value > 0
#' @details Only need three of S, N, ra, and Sa. Keep unit the same among these variables
#' @return Number of exepected events
#' @examples E_piece1(S = 24, Sa = 12, r.seq = c(10,15,20,30,25), lambda = log(2)/20, eta = 0)
#' @export
#' @references Kim, Kyungmann, and Anastasios A. Tsiatis. ``Study duration for clinical trials with survival response and early stopping rule." Biometrics (1990): 81-92.

E_piece1 <- function(S, Sa, r.seq, lambda, eta){
  t.seq   = unique(c(seq(0, min(S,Sa),1), min(S, Sa)))
  t.start = t.seq[1:(length(t.seq)-1)]
  t.end   = t.seq[2:length(t.seq)]
  if (length(r.seq) > length(t.start)) {
    r.seq <- r.seq[1:length(t.start)]
  } else if (length(r.seq) < length(t.start)){
    r.seq <- c(r.seq, rep(r.seq[length(r.seq)], length(t.start) - length(r.seq)))
  }

  ee = NULL
  for (i in seq(length(t.start))) {
    ee = c(ee, r.seq[i]*(H_t(lambda, eta, t.end[i], S) - H_t(lambda, eta, t.start[i], S)))
  }
  ee = sum(ee)

  return((ee)) # use round here try to become neutral, ceiling is too aggressive (thinking from solving equation, it is easier to reach required number)
}
