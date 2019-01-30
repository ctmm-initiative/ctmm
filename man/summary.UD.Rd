\name{summary.UD}
\alias{summary.UD}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarize a range distribution}
\description{ This function returns a list of biologically interesting parameters in human readable format, as derived from an autocorrelated kernel density estimate.}
\usage{
\S3method{summary}{UD}(object,level=0.95,level.UD=0.95,units=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An \code{akde} autocorrelated kernel-density estimate from the output of \code{akde}. }
  \item{level}{ Confidence level for the above area estimate. E.g., the 95\% confidence interval of the 50\% core area.}
  \item{level.UD}{ Coverage level for the home-range area. E.g., the 50\% core area. }
  \item{units}{Convert result to natural units.}
  \item{...}{Unused options.}
}

%\details{}

\value{
A \code{matrix} with low, maximum likelihood, and high estimates for the following parameters:
\describe{
\item{\code{area}}{The home-range area with fraction of inclusion \code{level.UD}. E.g., the 50\% core home range is estimated with \code{level.UD=0.50}, and 95\% confidence intervals are placed on that area estimate with \code{level=0.95}.

This kernel density estimate differs from the Gaussian estimate of \code{\link{summary.ctmm}}. The Gaussian estimate has more statistical efficiency, but is less related to space use for non-Gaussian processes.}
}
}

\references{
C. H. Fleming, J. M. Calabrese.
A new kernel-density estimator for accurate home-range and species-range area estimation.
\href{http://onlinelibrary.wiley.com/wol1/doi/10.1111/2041-210X.12673/abstract}{Methods in Ecology and Evolution, 8:5, 571-579 (2016)}.
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
Cilla <- buffalo$Cilla

# Fit a movement model
GUESS <- ctmm.guess(Cilla,interactive=FALSE)
FIT <- ctmm.fit(Cilla,GUESS)

# Estimate and summarize the AKDE
UD <- akde(Cilla,FIT)
summary(UD)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
