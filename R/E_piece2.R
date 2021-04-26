#' @title Calculate the expected number of events in two-arm studies (uniform)
#' @description A function to compute the number of events in two-arm studies under the assumption of exponential distribution for survival time and loss to follow-up time
#' @param S study duration, or equivalently, accrual time (Sa) + follow-up time (Sf)
#' @param r.seq Sequence of piecewise accrual rate. If length(ra.seq)>Sa, it directly assumes using ra.seq[-1] (the last one) as the accrual rate for the rest accrual period
#' @param Sa accrual duration, a numeric value
#' @param lambda1 rate parameter for event time (which is an exponential distribution) in arm#1, a numeric value > 0
#' @param eta1 rate parameter for censoring time (which is an exponential distribution) in arm#1, a numeric value > 0
#' @param lambda2 similar to lambda1, used for arm#2
#' @param eta2 similar to eta1, used for arm#2
#' @param alloc allocation between two arms, i.e. arm 1/(arm1 + arm2), a numeric value in (0,1)
#' @param CEILING if true, round-up the number of expected events
#' @return the number of exepected events, could be round up by specifying CEILING = TRUE
#' @examples
#' E_piece2(S = 24, Sa = 12, r.seq = c(10,15,20,30,25),
#' lambda1 = log(2)/20, lambda2 = log(2)/10, alloc = 0.66)
#' @export
#' @references Kim, Kyungmann, and Anastasios A. Tsiatis. ``Study duration for clinical trials with survival response and early stopping rule." Biometrics (1990): 81-92.

E_piece2 <- function(S = NULL, Sa = NULL, r.seq = NULL, lambda1, lambda2, eta1=0, eta2=0, alloc=0.5, CEILING = FALSE){
  ee0 <- E_piece1(S, Sa, r.seq, lambda1, eta1)
  ee1 <- E_piece1(S, Sa, r.seq, lambda2, eta2)
  ee  <- alloc * ee0 + (1-alloc) * ee1
  if ( CEILING == TRUE ) {
    ee  <- ceiling(ee)
  }
  return(ee)
}
