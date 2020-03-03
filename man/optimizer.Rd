\name{optimizer}
\alias{optimizer}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Minimize a function}
\description{This function serves as a wrapper around \code{\link{optimize}}, \code{\link{optim}}, and \code{ctmm}'s partial-Newton optimization routine, with standardized arguments and return values. It finds the optimal parameters that minimize a function, whether it be a cost, loss, risk, or negative log-likelihood function.}
\usage{
optimizer(par,fn,...,method="pNewton",lower=-Inf,upper=Inf,period=FALSE,reset=identity,control=list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{Initial parameter guess.}
  \item{fn}{Function to be minimized with first argument \code{par} and optional argument \code{zero} (see 'Details' below).}
  \item{...}{Optional arguments fed to \code{fn}.}
  \item{method}{Optimization algorithm (see 'Details' below).}
  \item{lower}{Lower bound for parameters.}
  \item{upper}{Upper bound for parameters.}
  \item{period}{Period of circular parameters if not \code{FALSE}.}
  \item{reset}{Optional function to re-center parameters, if symmetry permits, to prevent numerical underflow.}
  \item{control}{Argument list for the optimization routine (see 'Details' below).}
}

\details{Only \code{method='pNewton'} will work in both one dimension and multiple dimensions. Any other \code{method} argument will be ignored in one dimension, in favor of \code{\link{optimize}} with a backup evaluation of \code{\link{nlm}} (under a log-link) for cases where \code{\link{optimize}} is known to fail. In multiple dimensions, methods other than \code{pNewton} include those detailed in \code{\link{optim}}.

\code{method='pNewton'} is \code{ctmm}'s partial-Newton optimizer, which is a quasi-Newton method that is more accurate than BFGS-based methods when the gradient of \code{fn} must be calculated numerically. In short, while BFGS-based methods provide a single rank-1 update to the Hessian matrix per iteration, the partial-Newton algorithm provides \code{length(par)+1} rank-1 updates to the Hessian matrix per iteration, at the same computational cost. Furthermore, \code{length(par)} of those updates have better numerical precision than the BFGS update, meaning that they can be used at smaller step sizes to obtain better numerical precision. The \code{pNewton} optimizer also supports several features not found in other \code{R} optimizers: the \code{zero} argument, the \code{period} argument, and parallelization.

The \code{zero} argument is an optional argument in \code{fn} supported by \code{method='pNewton'}. Briefly, if you rewrite a negative log-likelihood of the form \eqn{fn = \sum_{i=1}^n fn_i} as \eqn{fn = \sum_{i=1}^n ( fn_i - zero/n ) + zero}, where \code{zero} is the current estimate of the minimum value of \code{fn}, then the sum becomes approximately "zeroed" and so the variance in numerical errors caused by the difference in magnitude between \code{fn} and \code{fn_i} is mitigated. In practice, without the \code{zero} argument, log-likelihood functions grow in magnitude with increasing data and then require increasing numerical precision to resolve the same differences in log-likelihood. But absolute differences in log-likelihoods (on the order of 1) are always important, even though most optimization routines more naturally consider relative differences as being important.

The \code{period} argument informs \code{method='pNewton'} if parameters is circular, such as with angles, and what their periods are.

The \code{control} list can take the folowing arguments, with defaults shown:
\describe{
\item{\code{precision=1/2}}{Fraction of machine numerical precision to target in the maximized likelihood value. The optimal \code{par} will have half this precision. On most computers, \code{precision=1} is approximately 16 decimal digits of precision for the objective function and 8 for the optimal \code{par}.}
\item{\code{maxit=.Machine$integer.max}}{Maximum number of iterations allowed for optimization.}
\item{\code{parscale=pmin(abs(par),abs(par-lower),abs(upper-par))}}{The natural scale of the parameters such that variations in \code{par} on the order of \code{parscale} produce variations in \code{fn} on the order of one.}
\item{\code{trace=FALSE}}{Return step-by-step progress on optimization.}
\item{\code{cores=1}}{Perform \code{cores} evaluations of \code{fn} in parallel, if running in UNIX. \code{cores<=0} will use all available cores, save \code{abs(cores)}. This feature is only supported by \code{method='pNewton'} and is only useful if \code{fn} is slow to evaluate, \code{length(par)>1}, and the total number of parallel evaluations required does not trigger fork-bomb detection by the OS.}
}


}

\value{Returns a list with components \code{par} for the optimal parameters, \code{value} for the minimum value of \code{fn}, and possibly other components depending on the optimization routine employed.}

%\references{}

\author{ C. H. Fleming. }

\note{\code{method='pNewton'} is very stringent about achieving its \code{precision} target and assumes that \code{fn} has small enough numerical errors (permitting the use of argument \code{zero}) to achieve that \code{precision} target. If the numerical errors in \code{fn} are too large, then the optimizer can fail to converge. \code{\link{ctmm.fit}} standardizes its input data before optimization, and back-transforms afterwards, as one method to minimize numerical errors in \code{fn}.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{optim}}, \code{\link{optimize}}, \code{\link{nlm}} }

%\examples{}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
