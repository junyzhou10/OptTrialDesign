#' @title Calculate the expected min(F(x)) via simulation
#' @description Calculate the probability of observing the median survival time from KM curve using Monte Carlo simulation of Ui given S, Sa, N, lambda, eta
#' @param S study duration, or equivalently, accrual time (Sa) + follow-up time (Sf)
#' @param Sa accrual duration, a numeric value
#' @param ra either uniform or piecewise accural rate. See details for more information
#' @param N number of subjects
#' @param lambda rate parameter for event time (assuming an exponential distribution), a numeric value > 0
#' @param eta rate parameter for loss to follow-up (assuming an exponential distribution), a numeric value > 0
#' @param nsim Number of simulations
#' @param seed Seed for reproducibility
#' @return A numeric value between (0,1), indicating the probability of observing the median survival time from KM curve. Or in other word, the last point of KM curve is under 0.5
#' @return Number of exepected events
#' @examples EF(S = 24, Sa = 12, ra = c(10,15,20,30,25), N = 320, lambda = log(2)/20, eta = 0)
#' @importFrom stats rexp runif
#' @export

EF <- function(S, Sa, ra = NULL, N, lambda, eta, nsim = 1000, seed=12345) {
  set.seed(seed)

  # Method #3: Vectorization, most efficient
  # Event time
  Ti = rexp(N*nsim, rate = lambda)
  #Enrollment time
  if (is.null(ra)) {
    Ai = runif(N*nsim, 0, Sa)
  } else {
    Pi = runif(N*nsim, 0, 1)
    Ai = sapply(Pi, function(p) inv.N(N, ra, p*N))
  }


  #Loss to follow up time
  if ( eta == 0 ) {
    #Ui is observed time
    Ui = pmin(pmax(0, S - Ai), Ti)
  } else {
    Ci = rexp(N*nsim, rate = eta)
    #Ui is observed time
    Ui = pmin(pmax(0, S - Ai), Ti, Ci)
  }

  # delta
  delta = ifelse(Ui == Ti, 1, 0)

  Ui    = matrix(Ui, nrow = N, byrow = F)
  delta = matrix(delta, nrow = N, byrow = F)
  sort.d= sapply(1:ncol(Ui), function(i) delta[,i][order(Ui[,i])])
  KMs   = apply(matrix(rep(1-1/(N - seq(1,N,1) + 1), nsim)^sort.d, nrow = N, byrow = F), 2, prod)
  KMs   = KMs[KMs!=0]
  return(mean(KMs <= 0.5))

  ## Method 2: use build in survfit function, but need to do it by loop. So it is slow
  # Ui = matrix(Ui, nrow = N, byrow = F)
  # delta = matrix(delta, nrow = N, byrow = F)
  # M = c()
  # for ( i in 1:nsim) {
  #   surv_object <- Surv(time = Ui[,i], event = delta[,i])
  #   fit1 <- survfit(surv_object ~ 1)
  #   M = c(M, summary(fit1)$surv[length(summary(fit1)$surv)])
  # }
  # return(mean(M <= 0.5))
}


