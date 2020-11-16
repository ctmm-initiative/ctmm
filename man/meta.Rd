\name{meta}
\alias{meta}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Meta-analysis of movement-model parameters}
\description{These functions estimate population-level parameters from individual movement models and related estimates, including AKDE home-range areas, while taking into account estimation uncertainty.
}
\usage{
meta(x,level=0.95,level.UD=0.95,method="MLE",IC="AICc",boot=FALSE,error=0.01,debias=TRUE,
     units=TRUE,plot=TRUE,sort=FALSE,mean=TRUE,col="black",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A list of \code{ctmm} movement-model objects, \code{UD} objects, or \code{UD} \code{summary} output, constituting a sampled population, or a list of such lists, each constituting a sampled sub-population.}
  \item{level}{Confidence level for parameter estimates.}
  \item{level.UD}{Coverage level for home-range estimates. E.g., 50\% core home range.}
%  \item{level.pop}{Coverage level for the population of individuals.}
  \item{method}{Statistical estimator used---either maximum likelihood estimation (\code{"MLE"}) or `best linear unbiased estimator' (\code{"BLUE"}).}
  \item{IC}{Information criterion to determine whether or not population variation can be estimated. Can be \code{"AICc"}, \code{AIC}, or \code{"BIC"}.}
  \item{boot}{Perform a parametric bootstrap for confidence intervals and first-order bias correction if \code{debias=TRUE}.}
  \item{error}{Relative error tolerance for parametric bootstrap.}
  \item{debias}{Apply Bessel's inverse-Gaussian correction if \code{method="MLE"}, REML if \code{method="BLUE"}, and an additional first-order correction if \code{boot=TRUE}.}
  \item{units}{Convert result to natural units.}
  \item{plot}{Generate a meta-analysis forest plot.}
  \item{sort}{Sort individuals by their point estimates in forest plot.}
  \item{mean}{Include population mean estimate in forest plot.}
  \item{col}{Color(s) for individual labels and error bars.}
  \item{...}{Further arguments passed to \code{plot}.}
}

\details{
So-far only the meta-analysis of home-range areas is implemented. More details will be provided in an upcomming manuscript.

For both estimator \code{methods}, the same underlying model is assumed.
}

\value{If \code{x} constitutes a sampled population, then \code{meta} returns a table with rows corresponding to population coverage levels and columns corresponding to confidence levels of the estimates.
If \code{x} constitutes a list of sampled sub-populations, then \code{meta} returns confidence intervals on the sub-population estimate ratios.}

%\section{Warnings}{}

%\references{}


\author{ C. H. Fleming. }

\note{The AICc formula is approximated via the Gaussian relation.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde}}, \code{\link{ctmm.fit}}. }

%\examples{\donttest{ }}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
