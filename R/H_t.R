#' @title An ancillary function for E_piece1 & E_piece2 (Not for all)
#' @description Function to calculate the integration, H(t)
#' @param lambda Event rate, assuming an exponential distribution
#' @param eta Loss to follow-up rate, assuming an exponential distribution
#' @param t current time
#' @param S study duration
#' @return A numeric value

H_t <- function(lambda, eta, t, S) {
  lambda/(lambda+eta)*(t - exp(-(lambda+eta)*(S-t))/(lambda+eta))
}
