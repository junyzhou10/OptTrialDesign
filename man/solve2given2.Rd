% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solve2given2.R
\name{solve2given2}
\alias{solve2given2}
\title{Given two of N, S, ra, Sa to figure out the rest two}
\usage{
solve2given2(
  Ea,
  S = NULL,
  N = NULL,
  Sa = NULL,
  ra = NULL,
  alloc = 0.5,
  lambda1,
  lambda2,
  eta1 = 0,
  eta2 = 0
)
}
\arguments{
\item{Ea}{Desired number of event to observe at the end of study}

\item{S}{study duration, or equivalently, accrual time (Sa) + follow-up time (Sf)}

\item{N}{number of subjects}

\item{Sa}{accrual duration, a numeric value}

\item{ra}{either uniform or piecewise accural rate. See details for more information}

\item{alloc}{allocation between two arms, i.e. arm 1/(arm1 + arm2), a numeric value in (0,1)}

\item{lambda1}{rate parameter for event time (which is an exponential distribution) in arm#1, a numeric value > 0}

\item{lambda2}{similar to lambda1, used for arm#2}

\item{eta1}{rate parameter for censoring time (which is an exponential distribution) in arm#1, a numeric value > 0}

\item{eta2}{similar to eta1, used for arm#2#'}
}
\value{
A list of (N, S, ra, Sa), which is a
}
\description{
A fundamental function to yield a study design to achieve the desired number of events (Ea) given any two of study duration (S), total sample size (N), accural rate (ra), and accrual period (Sa)
}
\details{
If ra is not provided, the function assumes a uniform accrual rate, as there is no information for other type of accrual.
If piecewise accrual rate, please notice that it is appropriate to give (ra, Sa) or (ra, N) which determines the accrual pattern.
If providing (ra, S), the solution is found as Sa equals to S
}
\examples{
# providing (S, Sa)
solve2given2(Ea = 88, S = 24, Sa = 12, lambda1 = log(2)/20, lambda2 = log(2)/10)

# providing (N, ra)
solve2given2(Ea = 88, N = 320, ra = c(10,15,20,30,25), lambda1 = log(2)/20, lambda2 = log(2)/10)
}
\references{
Kim, Kyungmann, and Anastasios A. Tsiatis. ``Study duration for clinical trials with survival response and early stopping rule." Biometrics (1990): 81-92.
}
