\name{Log}
\alias{Log}
\alias{Exp}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Log transformation of parameter estimates and their uncertainties}
\description{Methods for log transforming individual parameter estimates and their uncertainty estimates for use in meta-analytic regression, and then back-transforming mean-log parameter estimates back to mean parameter estimates.
}
\usage{
Log(x,variable="area",debias=TRUE,...)

Exp(est,VAR.est=0,VAR=0,VAR.VAR=0,variable="area",debias=TRUE,level=0.95,units=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A list of \code{UD} objects, \code{UD} \code{summary} objects, or \code{speed} objects.}
  \item{variable}{Can be \code{"area"}, \code{"diffusion"}, \code{"speed"}, \code{"tau position"}, or \code{"tau velocity"}.}
  \item{debias}{Apply \eqn{\log\chi^2} and \eqn{\log\chi} bias corrections if \code{TRUE}.}
  \item{...}{Further arguments passed.}
  \item{est}{Point estimate of the mean log-parameter.}
  \item{VAR.est}{Uncertainty in the mean log-parameter estimate (square standard error).}
  \item{VAR}{Variance in the log-parameters.}
  \item{VAR.VAR}{Uncertainty in the log-paramter variance estimate (square standard error).}
  \item{level}{Confidence level for parameter estimates.}
  \item{units}{Convert result to natural units.}
}

%\details{ }

\value{\code{Log} returns a list with two slots, \code{log} and \code{VAR.log}, corresponding to the point estimates and variance estimates of the logged variables.

\code{Exp} returns a confidence intervals for the back-transformed mean parameter estimate.}

%\section{Warnings}{}

%\references{}


\author{ C. H. Fleming. }

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{meta}}, \code{\link{mean}}. }

\examples{\donttest{
# load package and data
library(ctmm)
data(buffalo)

# fit movement models
FITS <- AKDES <- list()
for(i in 1:length(buffalo))
{
  GUESS <- ctmm.guess(buffalo[[i]],interactive=FALSE)
  # use ctmm.select unless you are certain that the selected model is OUF
  FITS[[i]] <- ctmm.fit(buffalo[[i]],GUESS)
}

# calculate AKDES on a consistent grid
AKDES <- akde(buffalo,FITS)

# extract 95\% areas
AREAS <- lapply(AKDES,summary)

# log transform for further meta-analysis
LOG <- Log(AREAS)

LOG
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
