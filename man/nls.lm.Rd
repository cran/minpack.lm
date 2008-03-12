\name{nls.lm}
\alias{nls.lm}
\encoding{UTF-8}
\title{Addresses NLS problems with the Levenberg-Marquardt algorithm}
\description{
  The purpose of \code{nls.lm} is to minimize the sum square of the
  vector returned by the function \code{fn}, by a modification of the
  Levenberg-Marquardt algorithm. The user may also provide a 
  function \code{jac} which calculates the Jacobian.
}
\usage{
nls.lm(par, fn, jac = NULL, control = nls.lm.control(), \dots)
}
\arguments{
 \item{par}{A list or numeric vector of starting estimates. If
   \code{par} is a list, then each element must be of length 1. }
 \item{fn}{A function that returns a vector of residuals, the sum square
   of which is to be minimized.  The first argument of \code{fn} must be
   \code{par}. }
 \item{jac}{A function to return the Jacobian for the \code{fn} function.}
 \item{control}{
    An optional list of control settings.  See \code{\link{nls.lm.control}} for
    the names of the settable control values and their effect.
}
 \item{\dots}{Further arguments to be passed to \code{fn} and \code{jac}.}
}
\details{
  Both functions \code{fn} and \code{jac} (if provided) must return
  numeric vectors. Length of the vector returned by \code{fn} must
  not be lower than the length of \code{par}. The vector returned by
  \code{jac} must have length equal to
  \eqn{length(\code{fn}(\code{par}, \dots))\cdot length(\code{par})}{%
       length(\code{fn}(\code{par}, \dots)) * length(\code{par})}.

  The \code{control} argument is a list;  see \code{\link{nls.lm.control}} for
    details.
 
    \bold{Successful completion.}\cr
    \cr
    The accuracy of \code{nls.lm} is controlled by the convergence
    parameters \code{ftol}, \code{ptol}, and \code{gtol}. These
    parameters are used in tests which make three types of comparisons
    between the approximation \eqn{par} and a solution
    \eqn{par_0}{par0}. \code{nls.lm} terminates when any of the tests
    is satisfied. If any of the convergence parameters is less than
    the machine precision, then \code{nls.lm} only attempts to satisfy
    the test defined by the machine precision. Further progress is not
    usually possible.\cr
    The tests assume that \code{fn} as well as \code{jac} are
    reasonably well behaved.  If this condition is not satisfied, then
    \code{nls.lm} may incorrectly indicate convergence. The validity
    of the answer can be checked, for example, by rerunning
    \code{nls.lm} with tighter tolerances.\cr
    \cr
    \emph{First convergence test.}\cr
    If \eqn{|z|} denotes the Euclidean norm of a vector \eqn{z}, then
    this test attempts to guarantee that
    \deqn{|fvec| < (1 + \code{ftol})\,|fvec_0|,}{%
              |fvec| < (1 + \code{ftol})|fvec0|,}
    where \eqn{fvec_0}{fvec0} denotes the result of \code{fn} function
    evaluated at \eqn{par_0}{par0}. If this condition is satisfied
    with \code{ftol} \eqn{\simeq 10^{-k}}{~ 10^(-k)}, then the final
    residual norm \eqn{|fvec|} has \eqn{k} significant decimal digits
    and \code{info} is set to 1 (or to 3 if the second test is also
    satisfied). Unless high precision solutions are required, the
    recommended value for \code{ftol} is the square root of the machine
    precision.\cr
    \cr
    \emph{Second convergence test.}\cr
    If \eqn{D} is the diagonal matrix whose entries are defined by the
    array \code{diag}, then this test attempt to guarantee that
        \deqn{|D\,(par - par_0)| < \code{ptol}\,|D\,par_0|,}{%
          |D*(par - par0)| < \code{ptol}|D*par0|,}
	If this condition is satisfied with \code{ptol} \eqn{\simeq
	  10^{-k}}{~ 10^(-k)}, then the larger components of
	\eqn{(D\,par)}{D*par} have \eqn{k} significant decimal digits and
	\code{info} is set to 2 (or to 3 if the first test is also
	satisfied). There is a danger that the smaller components of
	\eqn{(D\,par)}{D*par} may have large relative errors, but if
	\code{diag} is internally set, then the accuracy of the components
	of \eqn{par} is usually related to their sensitivity. Unless high
	precision solutions are required, the recommended value for
	\code{ptol} is the square root of the machine precision.\cr
	\cr
	\emph{Third convergence test.}\cr
    This test is satisfied when the cosine of the angle between the
    result of \code{fn} evaluation \eqn{fvec} and any column of the
    Jacobian at \eqn{par} is at most \code{gtol} in absolute value.
    There is no clear relationship between this test and the accuracy
    of \code{nls.lm}, and furthermore, the test is equally well
    satisfied at other critical points, namely maximizers and saddle
    points.  Therefore, termination caused by this test (\code{info} =
    4) should be examined carefully. The recommended value for
    \code{gtol} is zero.\cr
    \cr
  \bold{Unsuccessful completion.}\cr
    \cr
    Unsuccessful termination of \code{nls.lm} can be due to improper
    input parameters, arithmetic interrupts, an excessive number of
    function evaluations, or an excessive number of iterations. \cr
    \cr
    \emph{Improper input parameters.}\cr
    \code{info} is set to 0 if \eqn{length(\code{par}) = 0}, or
    \eqn{length(fvec) < length(\code{par})}, or \code{ftol} \eqn{< 0},
    or \code{ptol} \eqn{< 0}, or \code{gtol} \eqn{< 0}, or \code{maxfev}
    \eqn{\leq 0}{<= 0}, or \code{factor} \eqn{\leq 0}{<= 0}.\cr
    \cr
    \emph{Arithmetic interrupts.}\cr
    If these interrupts occur in the \code{fn} function during an
    early stage of the computation, they may be caused by an
    unacceptable choice of \eqn{par} by \code{nls.lm}. In this case,
    it may be possible to remedy the situation by rerunning
    \code{nls.lm} with a smaller value of \code{factor}.\cr
    \cr
    \emph{Excessive number of function evaluations.}\cr
    A reasonable value for \code{maxfev} is \eqn{100\cdot
    (length(\code{par}) + 1)}{100*(length(\code{par}) + 1)}. If the
    number of calls to \code{fn} reaches \code{maxfev}, then this
    indicates that the routine is converging very slowly as measured
    by the progress of \eqn{fvec} and \code{info} is set to 5. In this
    case, it may be helpful to force \code{diag} to be internally set.

    \emph{Excessive number of function iterations.}\cr
    The allowed number of iterations defaults to 50, can be increased if
    desired. \cr

    The list returned by \code{nls.lm} has methods 
    for the generic functions \code{\link{coef}},
    \code{\link{deviance}}, \code{\link{df.residual}},
    \code{\link{print}}, \code{\link{residuals}}, \code{\link{summary}},
    and \code{\link{vcov}}.

  }
\value{
  A list with components:
  \item{par}{The best set of parameters found.}
  \item{hessian}{A symmetric matrix giving an estimate of the Hessian
    at the solution found.}
  \item{fvec}{The result of the last \code{fn} evaluation; that is, the
  residuals. }
  \item{info}{\code{info} is an integer code indicating
    the reason for termination.
    \describe{
      \item{0}{Improper input parameters.}
    \item{1}{Both actual and predicted relative reductions in the
      sum of squares are at most \code{ftol}.}
    \item{2}{Relative error between two consecutive iterates is
      at most \code{ptol}.}
    \item{3}{Conditions for \code{info} = 1 and \code{info} = 2 both hold.}
    \item{4}{The cosine of the angle between \code{fvec} and any column
      of the Jacobian is at most \code{gtol} in absolute value.}
    \item{5}{Number of calls to \code{fn} has reached \code{maxfev}.}
    \item{6}{\code{ftol} is too small. No further reduction in the sum
      of squares is possible.}
    \item{7}{\code{ptol} is too small. No further improvement in the
      approximate solution \code{par} is possible.}
    \item{8}{\code{gtol} is too small. \code{fvec} is orthogonal to the
      columns of the Jacobian to machine precision.}
    \item{9}{The number of iterations has reached \code{maxiter}.}
  }}
  \item{message}{character string indicating reason for termination}.
  \item{diag}{The result list of \code{diag}. See \bold{Details}.}
  \item{niter}{The number of iterations completed before termination.}
  \item{deviance}{The sum of the squared residual vector.}
}
\references{
  J.J. Moré, "The Levenberg-Marquardt algorithm: implementation and
  theory," in \emph{Lecture Notes in Mathematics}
  \bold{630}: Numerical Analysis, G.A. Watson (Ed.),
  Springer-Verlag: Berlin, 1978, pp. 105-116.
}
\note{
  The public domain FORTRAN sources of MINPACK package by J.J. Moré,
  implementing the Levenberg-Marquardt algorithm were downloaded from
  \url{http://ftp.netlib.org/minpack}, and left largely unchanged. 
  The contents of this manual page are largely extracted from
  the comments of MINPACK sources.
}

\seealso{\code{\link{optim}}, \code{\link{nls}}, \code{\link{nls.lm.control}}}

\examples{
f <- function(T, tau, N0, a, f0) {
    expr <- expression(N0*exp(-T/tau)*(1 + a*cos(f0*T)))
    eval(expr)
}
j <- function(T, tau, N0, a, f0) {
    expr <- expression(N0*exp(-T/tau)*(1 + a*cos(f0*T)))
    c(eval(D(expr, "tau")),
      eval(D(expr, "N0" )),
      eval(D(expr, "a"  )),
      eval(D(expr, "f0" )))
}

T <- seq(0, 8, len=501)
p <- c("tau" = 2.2, "N0" = 1000, "a" = 0.25, "f0" = 8)

N <- do.call("f", c(list(T = T), as.list(p)))
N <- rnorm(length(N), mean=N, sd=sqrt(N))

plot(T, N, bg = "black", pch = 21, cex = 0.5)

fcn     <- function(p, T, N, N.Err, fcall, jcall)
    (N - do.call("fcall", c(list(T = T), as.list(p))))/N.Err
fcn.jac <- function(p, T, N, N.Err, fcall, jcall) {
    N.Err <- rep(N.Err, length(p))
    -do.call("jcall", c(list(T = T), as.list(p)))/N.Err
}

guess <- c("tau" = 2.2, "N0" = 1500, "a" = 0.25, "f0" = 10)

out <- nls.lm(par = guess, fn = fcn, #jac = fcn.jac,
              fcall = f, jcall = j,
              T = T, N = N, N.Err = sqrt(N), 
              control = nls.lm.control(nprint=1, diag = numeric()))

N1 <- do.call("f", c(list(T = T), out$par))  # N1 == N - sqrt(N)*out$fvec

lines(T, N1, col="blue", lwd=2)

summary(out) 
}
\keyword{nonlinear}
\keyword{optimize}
\keyword{regression}

