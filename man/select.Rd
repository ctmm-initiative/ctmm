\name{select}
\alias{lasso}
\alias{marquee}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spatial selection methods for telemetry objects.}
\description{Methods to segment or subset telemety objects based on polygon lasso and rectangular marquee selectors.}
\usage{
lasso(object,...)

marquee(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{telemetry} object or list of such objects.}
  \item{...}{Additional arguments passed to \code{plot}.}
}

\details{Allows the user to subset telemetry data into two groups (interior and exterior), based on a hand-drawn polygon lasso or rectangular marquee.}

\value{Returns a named list telemetry objects, twice the length of the input \code{object}, where the first half are the interior subsets and the second half are the exterior subsets.}

%\references{}

\author{C. H. Fleming.}

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{plot.telemetry}}  }

%% This example is interactive!
\examples{\donttest{
# Load package and data
library(ctmm)
data(wolf)

# Extract wolf Luna
Luna <- wolf$Luna

# Select resident data
Luna.sub <- lasso(Luna)

# You can now work with the resident and dispersive data separately
names(Luna.sub)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
