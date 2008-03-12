\name{nls.lm.control}
\alias{nls.lm.control}
\encoding{UTF-8}
\title{Control various aspects of the Levenberg-Marquardt algorithm}
\description{
  Allow the user to set some characteristics 
  Levenberg-Marquardt nonlinear least squares algorithm implemented
  in \code{nls.lm}.
}
\usage{
nls.lm.control(ftol = sqrt(.Machine$double.eps),
ptol = sqrt(.Machine$double.eps), gtol = 0, diag = list(), epsfcn = 0,
factor = 100, maxfev = integer(), maxiter = 50, nprint = 0)
}
\arguments{
  \item{ftol}{non-negative numeric. Termination occurs when
      both the actual and predicted relative reductions in the sum of
      squares are at most \code{ftol}. Therefore, \code{ftol} measures
      the relative error desired in the sum of squares.}
    \item{ptol}{non-negative numeric. Termination occurs when
      the relative error between two consecutive iterates is at most
      \code{ptol}. Therefore, \code{ptol} measures the relative error
      desired in the approximate solution.}
    \item{gtol}{non-negative numeric. Termination occurs when
      the cosine of the angle between result of \code{fn} evaluation
      \eqn{fvec} and any column of the Jacobian is at most \code{gtol}
      in absolute value. Therefore, \code{gtol} measures the
      orthogonality desired between the function vector and the
      columns of the Jacobian.}
    \item{diag}{a list or numeric vector containing positive
      entries that serve as multiplicative scale factors for the
      parameters. Length of \code{diag} should be equal to that of
      \code{par}. If not, user-provided \code{diag} is ignored and
      \code{diag} is internally set.}
    \item{epsfcn}{(used if \code{jac} is not provided) is a
      numeric used in determining a suitable step for the
      forward-difference approximation. This approximation assumes
      that the relative errors in the functions are of the order of
      \code{epsfcn}. If \code{epsfcn} is less than the machine
      precision, it is assumed that the relative errors in the
      functions are of the order of the machine precision.}
    \item{factor}{positive numeric, used in determining the
      initial step bound.  This bound is set to the product of
      \code{factor} and the \eqn{|\code{diag}*\code{par}|} if nonzero,
      or else to \code{factor} itself. In most cases \code{factor}
      should lie in the interval (0.1,100). 100 is a generally
      recommended value.}
    \item{maxfev}{integer; termination occurs
      when the number of calls to \code{fn} has reached \code{maxfev}.
      Note that \code{nls.lm} sets the value of \code{maxfev} to 
      \code{100*(length(par) + 1)} if 
      \code{maxfev = integer()}, where \code{par} is the list or
      vector of parameters to be optimized.  }
    \item{maxiter}{positive integer. Termination occurs
      when the number of iterations reaches \code{maxiter}.}
    \item{nprint}{is an integer, that enables controlled
      printing of iterates if it is positive. In this case, estimates
      of \code{par} are printed at the beginning of the first
      iteration and every \code{nprint} iterations thereafter and
      immediately prior to return. If \code{nprint} is not positive,
      no tracing information on the progress of the optimization is
      produced.}
  }
  \value{
    A \code{list} with exactly nine components:
    \item{ftol}{}
    \item{ptol}{}
    \item{gtol}{}
    \item{diag}{}
    \item{epsfcn}{}
    \item{factor}{}
    \item{maxfev}{}
    \item{nprint}{}
    with meanings as explained under \sQuote{Arguments}.
}
\references{
  J.J. Moré, "The Levenberg-Marquardt algorithm: implementation and
  theory," in \emph{Lecture Notes in Mathematics}
  \bold{630}: Numerical Analysis, G.A. Watson (Ed.),
  Springer-Verlag: Berlin, 1978, pp. 105-116.
}

\seealso{\code{\link{nls.lm}}}
\examples{
nls.lm.control(maxiter = 4)
}
\keyword{nonlinear}
\keyword{optimize}
\keyword{regression}

