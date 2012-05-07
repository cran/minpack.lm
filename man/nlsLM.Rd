\name{nlsLM}
\alias{nlsLM}
\encoding{UTF-8}
\title{Standard 'nls' framework that uses 'nls.lm' for fitting}
\description{
  \code{nlsLM} is a modified version of \code{\link{nls}} that uses \code{nls.lm} for fitting. 
  Since an object of class 'nls' is returned, all generic functions such as \code{\link{anova}}, 
  \code{\link{coef}}, \code{\link{confint}}, \code{\link{deviance}}, \code{\link{df.residual}}, 
  \code{\link{fitted}}, \code{\link{formula}}, \code{\link{logLik}}, \code{\link{predict}},
  \code{\link{print}}, \code{\link{profile}}, \code{\link{residuals}}, \code{\link{summary}},
  \code{\link{update}}, \code{\link{vcov}} and \code{\link{weights}} are applicable.  
}

\usage{
nlsLM(formula, data = parent.frame(), start, jac = NULL, 
      algorithm = "LM", control = nls.lm.control(), 
      lower = NULL, upper = NULL, trace = FALSE, subset, 
      weights, na.action, model = FALSE, \dots)
}

\arguments{
   \item{formula}{a nonlinear model \code{\link{formula}} including variables and parameters. Will be coerced to a formula if necessary.}
  \item{data}{an optional data frame in which to evaluate the variables in\code{formula} and \code{weights}.  
  Can also be a list or an environment, but not a matrix.}
  \item{start}{a named list or named numeric vector of starting estimates.}
  \item{jac}{A function to return the Jacobian.}
  \item{algorithm}{only method \code{"LM"} (Levenberg-Marquardt) is implemented.}
  \item{control}{an optional list of control settings.  See \code{\link{nls.lm.control}} for 
  the names of the settable control values and their effect.}
  \item{lower}{A numeric vector of lower bounds on each parameter. If not given, the default lower bound for each parameter is set to \code{-Inf}.}
  \item{upper}{A numeric vector of upper bounds on each parameter. If not given, the default upper bound for each parameter is set to \code{Inf}.}
  \item{trace}{logical value indicating if a trace of the iteration progress should be printed. Default is \code{FALSE}. If \code{TRUE}, the residual (weighted) sum-of-squares and the parameter values are printed at the conclusion of each iteration.}
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
  \item{weights}{an optional numeric vector of (fixed) weights.  When present, the objective function is weighted least squares.}
  \item{na.action}{a function which indicates what should happen when the data contain \code{NA}s. The default is set by the \code{na.action} setting of \code{\link{options}}, and is \code{\link{na.fail}} if that is unset. The \sQuote{factory-fresh} default is \code{\link{na.omit}}. Value \code{\link{na.exclude}} can be useful.}
  \item{model}{logical. If true, the model frame is returned as part of the object. Default is \code{FALSE}.}
  \item{\dots}{Additional optional arguments. None are used at present.}
}

\details{
The standard \code{\link{nls}} function was modified in several ways to incorporate the Levenberg-Marquardt type \code{\link{nls.lm}} fitting algorithm. The \code{formula} is transformed into a function that returns a vector of (weighted) residuals whose sum square is minimized by \code{\link{nls.lm}}. The optimized parameters are then transferred to \code{stats:::nlsModel} in order to obtain an object of class 'nlsModel'. The internal C function \code{C_nls_iter} and \code{stats:::nls_port_fit} were removed to avoid subsequent "Gauss-Newton", "port" or "plinear"" types of optimization of \code{nlsModel}. Several other small modifications were made in order to make all generic functions work on the output.
}

\value{
  A list of
  \item{m}{an \code{nlsModel} object incorporating the model.}
  \item{data}{the expression that was passed to \code{nls} as the data argument.  The actual data values are present in the environment of the \code{m} component.}
  \item{call}{the matched call.}
  \item{convInfo}{a list with convergence information.}
  \item{control}{the control \code{list} used, see the \code{control} argument.}
  \item{na.action}{the \code{"na.action"} attribute (if any) of the model frame.}
  \item{dataClasses}{the \code{"dataClasses"} attribute (if any) of the \code{"terms"} attribute of the model frame.}
  \item{model}{if \code{model = TRUE}, the model frame.}
  \item{weights}{if \code{weights} is supplied, the weights.} 
}

\references{
Bates, D. M. and Watts, D. G. (1988)
\emph{Nonlinear Regression Analysis and Its Applications}, Wiley

Bates, D. M. and Chambers, J. M. (1992)
\emph{Nonlinear models.}
Chapter 10 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.

J.J. MorÃ©, "The Levenberg-Marquardt algorithm: implementation and theory," 
in \emph{Lecture Notes in Mathematics}
\bold{630}: Numerical Analysis, G.A. Watson (Ed.), Springer-Verlag: Berlin, 1978, pp. 105-116.
}

\author{Andrej-Nikolai Spiess and Katharine M. Mullen}

\seealso{\code{\link{nls.lm}}, \code{\link{nls}}, \code{\link{nls.lm.control}}, \code{\link{optim}}}

\examples{
### Examples from 'nls' doc ###
DNase1 <- subset(DNase, Run == 1)
## using a selfStart model
fm1DNase1 <- nlsLM(density ~ SSlogis(log(conc), Asym, xmid, scal), DNase1)
## using logistic formula
fm2DNase1 <- nlsLM(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
                 data = DNase1, 
                 start = list(Asym = 3, xmid = 0, scal = 1))

## all generics are applicable
coef(fm1DNase1)
confint(fm1DNase1)
deviance(fm1DNase1)
df.residual(fm1DNase1)
fitted(fm1DNase1)
formula(fm1DNase1)
logLik(fm1DNase1)
predict(fm1DNase1)
print(fm1DNase1)
profile(fm1DNase1)
residuals(fm1DNase1)
summary(fm1DNase1)
update(fm1DNase1)
vcov(fm1DNase1)
weights(fm1DNase1)

## weighted nonlinear regression using 
## inverse squared variance of the response
## gives same results as original 'nls' function
Treated <- Puromycin[Puromycin$state == "treated", ]
var.Treated <- tapply(Treated$rate, Treated$conc, var)
var.Treated <- rep(var.Treated, each = 2)
Pur.wt1 <- nls(rate ~ (Vm * conc)/(K + conc), data = Treated, 
               start = list(Vm = 200, K = 0.1), weights = 1/var.Treated^2)
Pur.wt2 <- nlsLM(rate ~ (Vm * conc)/(K + conc), data = Treated, 
               start = list(Vm = 200, K = 0.1), weights = 1/var.Treated^2)
all.equal(coef(Pur.wt1), coef(Pur.wt2))

## 'nlsLM' can fit zero-noise data
## in contrast to 'nls'
x <- 1:10
y <- 2*x + 3                           
\dontrun{
nls(y ~ a + b * x, start = list(a = 0.12345, b = 0.54321))
}
nlsLM(y ~ a + b * x, start = list(a = 0.12345, b = 0.54321))

### Examples from 'nls.lm' doc
## values over which to simulate data 
x <- seq(0,5, length = 100)
## model based on a list of parameters 
getPred <- function(parS, xx) parS$a * exp(xx * parS$b) + parS$c 
## parameter values used to simulate data
pp <- list(a = 9,b = -1, c = 6) 
## simulated data with noise  
simDNoisy <- getPred(pp, x) + rnorm(length(x), sd = .1)
## make model
mod <- nlsLM(simDNoisy ~ a * exp(b * x) + c, 
             start = c(a = 3, b = -0.001, c = 1), 
             trace = TRUE)     
## plot data
plot(x, simDNoisy, main = "data")
## plot fitted values
lines(x, fitted(mod), col = 2, lwd = 2)

## create declining cosine
## with noise
TT <- seq(0, 8, length = 501)
tau <- 2.2
N0 <- 1000
a <- 0.25
f0 <- 8
Ndet <- N0 * exp(-TT/tau) * (1 + a * cos(f0 * TT))
N <- Ndet +  rnorm(length(Ndet), mean = Ndet, sd = .01 * max(Ndet))
## make model
mod <- nlsLM(N ~ N0 * exp(-TT/tau) * (1 + a * cos(f0 * TT)), 
             start = c(tau = 2.2, N0 = 1500, a = 0.25, f0 = 10), 
             trace = TRUE)  

## plot data
plot(TT, N, main = "data")
## plot fitted values
lines(TT, fitted(mod), col = 2, lwd = 2)
}
\keyword{nonlinear}
\keyword{optimize}
\keyword{regression}
