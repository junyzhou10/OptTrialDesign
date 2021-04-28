#' @title Calculate the desired number of events at the end of study
#' @description A function to compute the number of events based on asymptotic theory of log-rank test
#' @param HR Hazard Ratio, a numeric value in (0, 1)
#' @param alpha type I error of one-sided log-rank test, a numeric value in (0, 1)
#' @param power power required to observe the significant result under current settings, a numeric value in (0, 1)
#' @param alloc allocation between two arms, i.e. arm 1/(arm1 + arm2), a numeric value in (0,1)
#' @return Number of exepected events, an integer value (rounding-up)
#' @examples E_a(HR=0.5, alpha=0.025, power=0.9, alloc=0.5)
#' @export
#' @references Schoenfeld, David. "The asymptotic properties of nonparametric tests for comparing survival distributions." Biometrika 68.1 (1981): 316-319.

E_a <- function(HR, alpha, power, alloc){
  if(HR > 1){
    HR <- 1/HR
  }
  num = (abs(qnorm(alpha)) + abs(qnorm(power)))^2
  denom = alloc*(1-alloc)*(log(HR))^2
  ceiling(num/denom*alloc) + ceiling(num/denom*(1-alloc))
}
