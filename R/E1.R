#' @title Calculate the expected number of events (single arm, uniform accrual rate)
#' @description A very fundamental function to calculate the expected number of events given sample size (N), study duration (S), and either one of accrual rate (ra) or accural duration (Sa)
#' @param S study duration, or equivalently, accrual time (Sa) + follow-up time (Sf)
#' @param N number of subjects
#' @param ra uniform accural rate, a numeric value
#' @param Sa accrual duration, a numeric value
#' @param lambda rate parameter for event time (assuming an exponential distribution), a numeric value > 0
#' @param eta rate parameter for loss to follow-up (assuming an exponential distribution), a numeric value > 0
#' @details Only need three of S, N, ra, and Sa. Keep unit the same among these variables
#' @return Number of exepected events under such trial design settings
#' @importFrom stats qnorm
#' @examples E1(S = 24, N = 120, ra = 10, lambda = log(2)/20, eta = 0)
#' @export
#' @references Kim, Kyungmann, and Anastasios A. Tsiatis. ``Study duration for clinical trials with survival response and early stopping rule." Biometrics (1990): 81-92.

E1 <- function(S = NULL, N = NULL, Sa = NULL, ra = NULL, lambda, eta){
  if ( is.null(ra) ) {
    if (Sa > S) {
      Sa = S
    }
    ee = N/Sa*lambda/(lambda + eta)*(Sa - exp(-(lambda + eta)*S)/(lambda + eta)*(exp((lambda + eta)*Sa) - 1))

  } else {
    if ( is.null(Sa) ) {
      Sa = N/ra
    }

    if (Sa > S) {
      Sa = S
    }
    ee = ra*lambda/(lambda + eta)*(Sa - exp(-(lambda + eta)*S)/(lambda + eta)*(exp((lambda + eta)*Sa) - 1))
  }
  return(ee)
}
