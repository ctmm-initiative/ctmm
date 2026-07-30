\name{ci}
\alias{ci.beta}
\alias{ci.chisq}
\alias{ci.invgauss}
\alias{ci.lognorm}
\alias{ci.norm}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Coverage Interval Functions}
\description{These functions calculate coverage intervals with various distributions.}
\usage{
# beta sampling distribution
ci.beta(est,var,se=sqrt(var),level=0.95,...)

# chi-squared sampling distribution
ci.chisq(est,var,se=sqrt(var),level=0.95,...)

# inverse-Gaussian sampling distribution
ci.invgauss(est,var,se=sqrt(var),level=0.95,...)

# log-normal sampling distribution
ci.lognorm(est,var,se=sqrt(var),level=0.95,...)

# normal sampling distribution
ci.norm(est,var,se=sqrt(var),level=0.95,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{est}{The point estimate.}
  \item{var}{The variance of the point estimate.}
  \item{sd}{The standard error of the point estimate.}
  \item{level}{The coverage probability of the confidence interval.}
  \item{alpha}{The }
  \item{...}{Optional arguments.}
}

\details{Parameter estimation and variable prediction often yield point estimates and variance estimates,
which are conventionally used to calculate confidence intervals and prediction intervals under the assumption of a normal distribution.
While this assumption is asymptotically optimal, it may be very biased and not respect boundaries for finite sample sizes.
These functions all calculate asymptotically optimal coverage intervals based on alternative distributions with a specified mean and variance.
}

\value{Each function returns a coverage interval array with the lower limit, point estimate, and upper limit.}

%\references{}

\author{ C. H. Fleming. }

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{}

\examples{
est <- 0.5
var <- 0.07

# support on (0,1)
ci.beta(est,var)

# support on (0,Inf)
ci.chisq(est,var)
ci.invgauss(est,var)
ci.lognorm(est,var)

#support on (-Inf,Inf)
ci.norm(est,var)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
