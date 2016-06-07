\name{summary.UD}
\alias{summary.UD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize a range distribution
}
\description{ This function returns a list of biologically interesting parameters in human readable format, as derived from an autocorrelated kernel density estimate.}
\usage{
\S3method{summary}{UD}(object,level.UD=0.95,level=0.95,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An \code{akde} autocorrelated kernel-density estimate from the output of \code{akde}. }
  \item{level.UD}{ Confidence level for the home-range area. E.g., the 50\% core area. }
  \item{level}{ Confidence level for the above area estimate. E.g., the 95\% confidence interval of the 50\% core area.}
  \item{...}{Unused options.}
}

%\details{}

\value{
A \code{matrix} with low, maximum likelihood, and high estimates for the following parameters:
\describe{
\item{\code{area}}{The home-range area with fraction of inclusion \code{level.UD}. E.g., the 50\% core home range is estimated with \code{level.UD=0.50}, and 95\% confidence intervals are placed on that area estimate with \code{level=0.95}.}
}
}

\references{
C. H. Fleming and W. F. Fagan and T. Mueller and K. A. Olson and P. Leimgruber and J. M. Calabrese (2015).
Rigorous home-range estimation with movement data: A new autocorrelated kernel-density estimator.
\href{http://www.esajournals.org/doi/abs/10.1890/14-2010.1}{Ecology, 96(5), 1182-1188}.
}

\author{ C. H. Fleming. }

\note{ Prior to \code{ctmm} v0.3.1, AKDEs included only errors due to autocorrelation uncertainty, which are insignificant in cases such as IID data.
Starting in v0.3.1, \code{akde} calculated an effective sample size \code{DOF.H} and used this to estimate area uncertainty under a chi-square approxmation.
Starting in v0.3.2, this method was improved to use \code{DOF.area} in the Gaussian reference function approximation.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde}}. }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)

# Extract movement data for a single animal
cilla <- buffalo[[1]]

# Find the best OU movement model
# also see help(variogram.fit)
GUESS <- ctmm(tau=6*24*60^2)
FIT <- ctmm.fit(cilla,GUESS)

# Estimate and summarize the AKDE
UD <- akde(cilla,FIT)
summary(UD)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
