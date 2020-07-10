\name{meta}
\alias{meta}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Meta-analysis of movement-model parameters}
\description{These functions estimate population-level parameters from individual movement models and related estimates, including AKDE home-range areas, while taking into account estimation uncertainty.
}
\usage{
meta(x,level=0.95,level.UD=0.95,level.pop=0.95,robust=FALSE,IC="AIC",boot=FALSE,
     error=0.01,debias=TRUE,units=TRUE,plot=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A list of \code{ctmm} movement-model objects, \code{UD} objects, or \code{UD} \code{summary} output.}
  \item{level}{Confidence level for parameter estimates.}
  \item{level.UD}{Coverage level for home-range estimates. E.g., 50\% core home range.}
  \item{level.pop}{Coverage level for the population of individuals.}
  \item{robust}{Estimate the population median rather than mean.}
  \item{IC}{Information criterion to determine whether or not population variation can be estimated. Can be \code{"AIC"}, \code{"BIC"}, \code{"LOOCV"}, or none (\code{NA}).}
  \item{boot}{Perform a parametric bootstrap to remove first-order bias and refine confidence intervals.}
  \item{error}{Relative error tolerance for parametric bootstrap.}
  \item{debias}{Use bootstrap results to remove first-order bias from point estimates.}
  \item{units}{Convert result to natural units.}
  \item{plot}{Generate a meta-analysis forest plot.}
  \item{...}{Further arguments passed to \code{plot}.}
}

\details{
So-far only the meta-analysis of home-range areas is implemented. More details will be provided in an upcomming manuscript.
}

\value{A table with rows corresponding to population coverage levels and columns corresponding to confidence levels of the estimates.}

%\section{Warnings}{}

%\references{}


\author{ C. H. Fleming. }

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde}}, \code{\link{ctmm.fit}}. }

%\examples{\donttest{ }}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
