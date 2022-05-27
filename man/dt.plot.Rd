\name{dt.plot}
\alias{dt.plot}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions for diagnosing sampling schedules}
\description{
Produces a log-scale plot of the sorted sampling intervals for inspection.
}
\usage{
dt.plot(data,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{telemetry} object.}
  \item{...}{Additional options passed to \code{plot}.}
}


\details{Horizontal lines are included at common sampling intervals (e.g., 1-hour) and dimmed horizontal lines are included at common subdivisions (e.g., 30-minutes).}

%\value{}

%\references{}

\author{C. H. Fleming.}

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{as.telemetry}}. }

\examples{
# Load package and data
library(ctmm)
data(gazelle)

# Plot the sampling intervals
dt.plot(gazelle)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
