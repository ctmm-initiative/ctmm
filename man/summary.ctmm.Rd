\name{summary.ctmm}
\alias{summary.ctmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize a continuous-time movement model
}
\description{ This function returns a list of biologically interesting parameters in human readable format, as derived from a continuous-time movement model.}
\usage{
\S3method{summary}{ctmm}(object,level=0.95,level.UD=0.95,units=TRUE,IC="AICc",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{ctmm} movement-model object from the output of \code{ctmm.fit}. }
  \item{level}{ Confidence level for parameter estimates. }
  \item{level.UD}{ Confidence level for the Gaussian home-range area. }
  \item{units}{Convert result to natural units.}
  \item{IC}{Information criteria for sorting lists of \code{ctmm} objects.}
  \item{...}{Unused options.}
}

%\details{}

\value{
If summary is called with a single \code{ctmm} object output from \code{\link{ctmm.fit}}, then a list is returned with effective sample size array \code{DOF} and parameter estimate table \code{CI}, with low, maximum likelihood, and high estimates for the following possible parameters:
\describe{
\item{\code{tau}}{The autocorrelation timescales.}
\item{\code{area}}{The Gaussian home-range area, where the point estimate has a significance level of \code{level.UD}. I.e., the core home range is where the animal is located 50\% of the time with \code{level.UD=0.50}.
This point estimate itself is subject to uncertainty, and is given confidence intervals derived from \code{level}.}
\item{\code{speed}}{The Gaussian root-mean-square (RMS) velocity, which is a convenient measure of average speed.}
}
If summary is called on a list of \code{ctmm} objects output from \code{\link{ctmm.select}}, then a table is returned with the model names and AIC differences, where "IID" denotes the uncorrelated bi-variate Gaussian model, "OU" denotes the continuous-position Ornstein-Uhlenbeck model, and "OUF" denotes the continuous-velocity Ornstein-Uhlenbeck-F model.
}

%\references{}

\author{ C. H. Fleming. }

\note{ Confidence intervals on the autocorrelation timescales assume they are sufficiently greater than zero and less than infinity.

In \code{ctmm} v0.3.4 the speed estimate was fixed to be the RMS velocity and not \eqn{1/\sqrt{2}} times the RMS velocity.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{ctmm.fit}}, \code{\link{ctmm.select}}. }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)

# Extract movement data for a single animal
cilla <- buffalo[[1]]

# Find the best OU movement model
# also see help(variogram.fit)
GUESS <- ctmm(tau=60*60*24*10)
FIT <- ctmm.fit(cilla,GUESS)

# Tell us something interpretable
summary(FIT)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
