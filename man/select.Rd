\name{select}
\alias{lasso}
\alias{marquee}
\alias{cleave}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spatial selection methods for telemetry objects.}
\description{Methods to segment or subset telemety objects based on polygon lasso, rectangular marquee, and time slider selectors.}
\usage{
lasso(object,...)

marquee(object,...)

cleave(object,fraction=0.5,name="CLEFT",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{telemetry} object or list of such objects.}
  \item{fraction}{Initial split index, as fraction of total (will change to time).}
  \item{name}{Name of list to store cleft telemetry objects to.}
  \item{...}{Additional arguments passed to \code{plot}.}
}

\details{\code{lasso} and \code{marquee} allow the user to subset telemetry data into two groups (interior and exterior), based on a hand-drawn polygon lasso or rectangular marquee. \code{cleave} allows the user to split the data into two halves at a particular time selected via slider.}

\value{\code{lasso} and \code{marquee} return a named list telemetry objects, twice the length of the input \code{object}, where the first half are the interior subsets and the second half are the exterior subsets. \code{cleave} stores a similar list of telemetry objects to \code{name} on button press.}

%\references{}

\author{C. H. Fleming.}

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{plot.telemetry}}  }

%% This example is interactive!
\examples{
# This example is interactive
if(interactive())
{
  # Load package and data
  library(ctmm)
  data(wolf)

  # Extract wolf Luna
  Luna <- wolf$Luna

  # Select resident data
  Luna.sub <- lasso(Luna)

  # You can now work with the resident and dispersive data separately
  names(Luna.sub)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
