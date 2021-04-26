#' @title Given two of N, S, ra, Sa to figure out the rest two
#' @description A fundamental function to yield a study design to achieve the desired number of events (Ea) given any two of study duration (S), total sample size (N), accural rate (ra), and accrual period (Sa)
#' @param Ea Desired number of event to observe at the end of study
#' @param S study duration, or equivalently, accrual time (Sa) + follow-up time (Sf)
#' @param N number of subjects
#' @param ra either uniform or piecewise accural rate. See details for more information
#' @param Sa accrual duration, a numeric value
#' @param lambda1 rate parameter for event time (which is an exponential distribution) in arm#1, a numeric value > 0
#' @param eta1 rate parameter for censoring time (which is an exponential distribution) in arm#1, a numeric value > 0
#' @param lambda2 similar to lambda1, used for arm#2
#' @param eta2 similar to eta1, used for arm#2#'
#' @param alloc allocation between two arms, i.e. arm 1/(arm1 + arm2), a numeric value in (0,1)
#' @details If ra is not provided, the function assumes a uniform accrual rate, as there is no information for other type of accrual.
#' If piecewise accrual rate, please notice that it is appropriate to give (ra, Sa) or (ra, N) which determines the accrual pattern.
#' If providing (ra, S), the solution is found as Sa equals to S
#' @return A list of (N, S, ra, Sa), which is a
#' @import stats
#' @export
#' @examples
#' # providing (S, Sa)
#' solve2given2(Ea = 88, S = 24, Sa = 12, lambda1 = log(2)/20, lambda2 = log(2)/10)
#'
#' # providing (N, ra)
#' solve2given2(Ea = 88, N = 320, ra = c(10,15,20,30,25), lambda1 = log(2)/20, lambda2 = log(2)/10)
#' @references Kim, Kyungmann, and Anastasios A. Tsiatis. ``Study duration for clinical trials with survival response and early stopping rule." Biometrics (1990): 81-92.


solve2given2 <- function(Ea, S=NULL, N=NULL, Sa=NULL, ra=NULL, alloc=0.5, lambda1, lambda2, eta1=0, eta2=0){
  CheckInputs = sum(c(is.null(S), is.null(N), is.null(Sa), is.null(ra)))
  if (CheckInputs != 2) {
    cat('ERROR: Please input exactly two of S, N, Sa, and ra!\n')
    return(NULL)
  }


  if ( is.null(ra) ){ # then a uniform ra is assumed, so it is more staightforward
    if (is.null(S)) { # use N, Sa
      ra = N / Sa
      obj <- function(s) E2(S = s, N, Sa, ra, lambda1, lambda2, eta1, eta2, alloc) - Ea
      S = uniroot(obj, interval = c(0,100), extendInt = 'upX')$root
    } else if (is.null(N)) { # use S, Sa
      obj <- function(s) E2(S, N, Sa, ra = s, lambda1, lambda2, eta1, eta2, alloc) - Ea
      ra = uniroot(obj, interval = c(0,100), extendInt = 'upX')$root
      N = ra * min(Sa, S)
    } else if (is.null(Sa)) { # use N, S
      obj <- function(s) E2(S, N, Sa, ra = s, lambda1, lambda2, eta1, eta2, alloc) - Ea
      ra = uniroot(obj, interval = c(0,100), extendInt = 'upX')$root
      Sa = min(N / ra, S)
    }
  } else { # ra is given
    # ra allows to be
    # 1. uniform: r_u;
    # 2. linear with max: (r_max, r0, r1);
    # 3. piecewise: (r_1, r_2, r_3,...r_k) after r_k, assume all equals to r_k
    if (length(ra) == 1) {# indicates r_u, then must have S, Sa (N is not needed)
      if (is.null(S)) { # ra + N/Sa
        Sa = ifelse(is.null(Sa), N/ra, Sa)
        obj <- function(s) E2(S = s, N, Sa, ra, lambda1, lambda2, eta1, eta2, alloc) - Ea
        S = uniroot(obj, interval = c(0,100), extendInt = 'upX')$root
        N = ifelse(is.null(N), ra*Sa, N)
      } else { # ra + S
        obj <- function(s) E2(S, N, Sa = s, ra, lambda1, lambda2, eta1, eta2, alloc) - Ea
        Sa = uniroot(obj, interval = c(0,S), extendInt = 'no')$root
        N = ra * Sa
      }
    } else { # both linear and piecewise could be solved in piecewise fashion
      if (length(ra) == 3) { # implies (r0, r1, r_max)
        r0    = ra[1]
        r1    = ra[2]
        r_max = ra[3]
        r.seq = c(seq(r0, r_max, r1), r_max)
      } else {
        r.seq = ra
      }

      if (!is.null(S)) {
        warning("When it is piecewise accrual rate, the solution is determined as Sa will be set to S, and all N, S is known.")
        Sa = S
        N = cum.N(Sa, r.seq)
      } else if (!is.null(N)) { # N, ra -> Sa
        Sa = inv.N(N, r.seq, N)
        obj <- function(s) E_piece2(S = s, Sa, r.seq, lambda1, lambda2, eta1, eta2, alloc) - Ea
        S = uniroot(obj, interval = c(0.1, 100), extendInt = 'upX')$root
      } else if (!is.null(Sa)) { # Sa, ra -> N
        obj <- function(s) E_piece2(S = s, Sa, r.seq, lambda1, lambda2, eta1, eta2, alloc) - Ea
        S = uniroot(obj, interval = c(0.1, 100), extendInt = 'upX')$root
        if (Sa >= length(r.seq)) {
          N = sum(r.seq) + (Sa-length(r.seq))*r.seq[length(r.seq)]
        } else if (Sa < length(r.seq)) {
          warning("Accrual more people than Sa, ra will be truncated!")
          cut.id = floor(Sa)
          N = sum(r.seq[1:cut.id]) + (Sa - cut.id)*r.seq[cut.id+1]
        }
      }
    }
  }

  return(list(S = round(S, 3), N = ceiling(N), Sa = round(Sa, 3), ra = round(ra, 3)))
}
