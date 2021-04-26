#' @title Find optimal clinical trial design
#' @description A function to seek the optimal trial design in terms of both financial benefits and regulatory consideration, i.e. clinical meaningful results, and data maturity requirements
#' @param alpha type I error rate (size)
#' @param POWER power (1-beta)
#' @param ra accrual rate (NULL if Sa is input)
#' @param Sa maximum enrollment period (NULL if ra is input), see details
#' @param Ea target event number to observe at the end of the study
#' @param lambda a vector of length 2, representing lambda (rate parameter of exponential distribution) for each arm (first treatment, then control)
#' @param eta a vector of length 2, representing eta (rate parameter of exponential distribution) for each arm (first treatment, then control)
#' @param alloc allocation between two arms, i.e. arm 1/(arm1 + arm2), a numeric value in (0,1)
#' @param c0 fixed cost
#' @param c1 cost per patient
#' @param c2 cost per unit time
#' @param b revenue per unit time
#' @param dL LOE(loss of exclusivity time period) - l0 (time between final analysis and market access)
#' @param Aj.Ind Indicator for clincial meaningful results. either "Median Difference" (A1) or "Median Ratio" (A2). If not used, no need to input
#' @param d0_r0 threshold value for clinical meaningful results, corresponding to the choice of Aj.Ind (see reference paper)
#' @param t0 threshold value for data maturity requirements C1. Either one input will activate the corresponding constraint
#' @param e0 threshold value for data maturity requirements C2
#' @param p0 threshold value for data maturity requirements C3 (may takes a while is specified due to simulation)
#' @param m0 threshold value for data maturity requirements C4
#' @param ... the number of simulations and seed when p0 is specified. The default is nsim = 1000, seed = 12345. Change if needed
#' @details The function requires to input one of accrual rate (ra) and accrual period (Sa). If Sa is provided,
#' a uniform accrual rate will be assumed.
#' Investigators may have a rough idea of the range of ra in general. Notice that under the assumptions, the larger ra always yield short study duration and thus a higher revenue, which is
#' more desired. Therefore, the optimal results must happen at max(ra). In reality, the accrual rate (ra) is always constrained by many limitations. So a natural choice is
#' to input max(ra) in the function, and obtain the optimal trial design from the ouput.
#' @return Optimal set that maximize the financial benefits but at meantime, yields clinical meaningful results & data maturity constraints .
#' The results are saved in a list with three elements:
#' \item{Res}{The optimal trial design. A list including optimal (S, N, ra, Sa), but also Ea, actual power (AccPower), data maturity indicators (id.c1-4), total sale (Sale), total cost (Cost), total revenue (Rev), probability of observing clinical meaningful results (PA)}
#' \item{All.dat}{A table of all potential trial designs}
#' \item{Valid.set}{A table of all valid trial designs, from which to find the optimal one}
#' @examples
#' findOpt(ra = NULL, Sa = 15, Ea=88, lambda=c(log(2)/20, log(2)/10), eta=c(log(2)/120, log(2)/100),
#' b=18, c0=10, c1=0.6, c2=0.1, dL=180, Aj.Ind = "Median Ratio", d0_r0 = 1.4, t0 = 6)
#'
#' # input ra, and more data constraints
#' findOpt(ra = c(10,20,30,40,36), Ea=88, lambda=c(log(2)/20, log(2)/10), eta=c(log(2)/120, log(2)/100),
#' b=18, c0=10, c1=0.6, c2=0.1, dL=180, Aj.Ind = "Median Ratio", d0_r0 = 1.4, t0 = 10, e0 = 0.6, p0 = 0.7, m0 = 6)
#' @importFrom stats pnorm rexp runif uniroot
#'
#' @export
#' @references (Nektar Therapeutics) An Optimal Design Strategy for Phase III Clinical Trials with Time-To-Event Endpoint
findOpt <- function(alpha=0.05, POWER=0.9, ra = NULL, Sa = NULL, Ea, lambda, eta, alloc = 0.5, c0, c1, c2, b, dL, d0_r0, Aj.Ind = "NotUsed", t0 = NA, e0 = NA, p0 = NA, m0 = NA, ...){
  PA = NULL
  Sa0 = Sa
  b0  = b # just some backup

  ## prepare ra.seq
  if (length(ra) == 3) { # implies (r0, r1, r_max)
    r0    = ra[1]
    r1    = ra[2]
    r_max = ra[3]
    ra0 = c(seq(r0, r_max, r1), r_max)
  } else {
    ra0 = ra
  }

  ## objective: max revenue / min cost
  if (!is.null(b) & !is.null(dL)) {
    obj.task = "Rev"
  } else {
    obj.task = "Cost"
  }

  if ( is.null(ra) + is.null(Sa) != 1 ) {
    warning('WARNING: Exactly one of ra or Sa should input!')
    return(NULL)
  }

  if ( length(lambda) == 2 ) {
    lambda1 = lambda[1]; lambda2 = lambda[2]
    if ( length(eta) == 2) {
      eta1 = eta[1]; eta2 = eta[2]
    } else {
      eta1 = eta2 = eta
    }


  } else if ( length(lambda) == 1 ) {
    lambda1 = lambda2 = lambda
    if (length(eta) == 2) {
      eta1 = eta[1]; eta2 = eta[2]
    } else {
      eta1 = eta2 = eta
    }
  } else {
    print('WARNING: length of lambda/eta should be at most 2!')
    return(NULL)
  }

  #########  For a line search, the upper/lower bound of N is critical;
  #### Lower bound: When S -> Inf, N value is determined
  #### Upper bound: 1. natural bound by S=Sa; 2. Or if data maturity constraints are activated, then min(N_natural, N_constraint)
  N.min <- ceiling(Ea/(alloc*lambda1 / (lambda1 + eta1) + (1-alloc)* lambda2/( lambda2 +  eta2)))
  if (is.null(ra)) {
    N.max = solve2given2(Ea = Ea, S = Sa0, Sa = Sa0, lambda1 = lambda1, lambda2 = lambda2, eta1 = eta1, eta2 = eta2)$N
  } else if (length(ra)== 1) {
    obj = function(x) {
      E2(S=x, Sa=x, ra = ra, lambda1 = lambda1, lambda2 = lambda2, eta1 = eta1, eta2 = eta2, alloc=alloc) - Ea
    }
    e.Sa = uniroot(obj, interval = c(0.1, 50), extendInt = 'upX')$root
    N.max = floor(e.Sa * ra)
  } else if (length(ra) > 1) {
    obj = function(x) {
      E_piece2(S=x, Sa=x, r.seq = ra, lambda1 = lambda1, lambda2 = lambda2, eta1 = eta1, eta2 = eta2, alloc=alloc) - Ea
    }
    e.Sa = uniroot(obj, interval = c(0.1, 50), extendInt = 'upX')$root
    N.max = floor(cum.N(e.Sa, ra))
  }
  N.seq = seq(N.min, N.max, 1) # available sets

  ## Previous N.max is globally true. But then we need to check data maturity to see if there is any further restriction
  id.c1 <- id.c2 <- id.c3 <- id.c4 <- 0 # some indicators; 1: activated ; 0: not activated


  # direct output all potential combinations; then we will further check c1,c3,c4
  out = sapply(N.seq, function(n) {
    new = solve2given2(Ea = Ea, N = n, Sa = Sa0, ra = ra0, lambda1 = lambda1, lambda2 = lambda2, eta1 = eta1, eta2 = eta2)
    N = new$N
    S = new$S
    Sa = new$Sa
    ra = new$ra

    if (length(ra0)>1) {
      E_1 = E_piece1(S=S, Sa=Sa0, r.seq=ra0, lambda=lambda1, eta=eta1)*alloc
      E_2 = E_piece1(S=S, Sa=Sa0, r.seq=ra0, lambda=lambda2, eta=eta2)*(1-alloc)
      ra  = NA
    } else {
      E_1 = E1(S=S, N, Sa=Sa0, ra=ra0, lambda=lambda1, eta=eta1)*alloc
      E_2 = E1(S=S, N, Sa=Sa0, ra=ra0, lambda=lambda2, eta=eta2)*(1-alloc)
    }

    if ( Aj.Ind == 'Median Difference') {
      d0 = d0_r0
      PA = 1 - pnorm( (d0 - log(2)/lambda1 + log(2)/lambda2)/(log(2)*sqrt(1/(lambda1^2*E_1)+1/(lambda2^2*E_2))) )
    } else if (Aj.Ind == 'Median Ratio') {
      r0 = d0_r0
      PA = 1 - pnorm( log(r0*lambda1/lambda2)/sqrt(1/E_1 + 1/E_2) )
    } else {
      PA = 1
    }

    rev = ifelse(obj.task == "Cost", - (c0 + c1*N + c2*S), PA*POWER*b*(dL-S) - (c0 + c1*N + c2*S))
    cost = c0 + c1*N + c2*S
    sale = rev + cost
    return(c(N, S, Sa, ra, PA, rev, cost, sale))
  })
  out0 <- out <- t(out)
  colnames(out0) <- colnames(out) <- c("N", "S", "Sa", "ra", "PA", "revenue", "costs", "sales")

  # Check for C2
  if ( !is.na(e0) ) {
    id.c2 = 1
    out = out[out[,1] <= Ea/e0, ]
  }

  # check C1
  if ( !is.na(t0) ) {
    id.c1 = 1
    out = out[(out[,2]-out[,3])>=t0, ]
    if ( dim(out)[1] < 1) {
      print('ERROR: Improper value of t0. Please change to a smaller one!')
      return(NULL)
    }
  }

  # Check for C4
  if ( !is.na(m0) ) {
    id.c4 = 1
    p.R = sapply(seq(dim(out)[1]), function(i){
      ignore = ifelse(is.na(out[i,4]), rr<-ra0, rr<-out[i,4])
      pr = 1 - min(1, max(0, cum.N(out[i,2]-m0, rr)/out[i,1] ) ) * (alloc*exp(-eta1 * m0) + (1-alloc)*exp(-eta2*m0))
      return(pr)
    })
    out = out[p.R<=0.5, ]
  }

  # Check for C3
  if ( !is.na(p0) ) {
    id.c3 = 1
    for (i in seq(dim(out)[1],1)) {
      ss = out[i, 2]; sa = out[i, 3]; nn = out[i, 1];
      ignore = ifelse(is.na(out[i,4]), rr<-ra0, rr<-out[i,4])
      minF = EF(S=ss, Sa=sa, ra = round(rr*alloc), N=round(nn*alloc), lambda = lambda1, eta = eta1, ...)*EF(S=ss, Sa=sa, ra = rr-round(rr*alloc), N=nn-round(nn*alloc), lambda = lambda2, eta = eta2, ...)
      if (minF > p0) {
        break
      }
    }
    if (i==1) {
      warning("Input p0 is not obtainable under such setting. Please try a smaller one! \n")
      break
    }
    out = out[1:i, ]
  }

  ## winner set:
  if (obj.task == "Rev") {
    opt.set = out[which.max(out[, 6]), ]
  } else {
    opt.set = out[which.min(out[, 7]), ]
  }

  AccPower = pnorm( sqrt(alloc*(1-alloc) * (log(lambda1/lambda2))^2 * Ea) - abs(qnorm(alpha)) )
  ignore   = ifelse(is.na(opt.set[4]), rr<-ra0, rr<-opt.set[4])
  Res      = list(AccPower = AccPower, N = opt.set[1], S = opt.set[2], Sa = opt.set[3], ra = rr, # notice ra could be NA as output which indicates piecewise
                  id.c1=id.c1, id.c2=id.c2, id.c3=id.c3, id.c4=id.c4, Sale = opt.set[8], Cost = opt.set[7], Rev = opt.set[6], PA =opt.set[5], N.min = N.min, N.max = N.max)
  return(list(Res = Res, All.dat = out0, Valid.set = out))
}

