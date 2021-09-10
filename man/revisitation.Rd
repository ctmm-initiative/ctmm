\name{revisitation}
\alias{revisitation}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate an revisitation distribution estimate}
\description{
This function estimates the distribution of revisitations from \code{telemetry} data and a continuous-time movement model.
}
\usage{
revisitation(data,UD,error=0.001,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ 2D timeseries telemetry data represented as a \code{telemetry} object or list of objects. }
  \item{UD}{A \code{UD} object from the output of \code{\link{akde}}.}
  \item{error}{Target probability error.}
  \item{...}{Arguments passed to \code{\link{akde}}.}
}

%\details{}

\value{Returns a \code{UD} object.}

%\references{}

\author{ C. H. Fleming. }

%\note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde}}, \code{\link{occurrence}} }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)
DATA <- buffalo$Cilla

# calculate fit guess object
GUESS <- ctmm.guess(DATA,interactive=FALSE)
# in general, you should be running ctmm.select here instead of ctmm.fit
FIT <- ctmm.fit(DATA,GUESS)

# Compute akde object
UD <- akde(DATA,FIT)

# compute revisitation distribution
RD <- revisitation(DATA,UD)

# Plot data with revisitation distribution
plot(DATA,RD)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
