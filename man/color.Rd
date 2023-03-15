\name{color}
\alias{annotate}
\alias{color}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Color telemetry objects by time}
\description{
These functions facilitate the coloring of tracks by annotating tracking data with time/location specific information and computing color arguments for \code{plot}.
}
\usage{
annotate(object,by="all",cores=1,...)

color(object,by="time",col.fn=NULL,alpha=1,dt=NULL,cores=1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{telemetry} object or list of objects. \code{color} can also take \code{ctmm} and \code{UD} objects.}
  \item{by}{What to \code{annotate} or \code{color} times by. Options include \code{"individual"}, \code{"time"}, \code{"sun"}, \code{"moon"}, \code{"season"}, and \code{"tropic"} (see Details below). \code{ctmm} and \code{UD} objects can only be colored by \code{"individual"}.}
  \item{col.fn}{Optional coloring function that can take a [0,1] interval and alpha channel argument.}
  \item{alpha}{Base alpha channel value.}
  \item{dt}{Sampling interval specification for making oversampled times more transparent. If \code{NULL}, the median will be used. Disabled if zero.}
  \item{cores}{Number of annotations or overlap calculations to peform in parallel. \code{cores=0} will use all cores, while \code{cores<0} will reserve \code{abs(cores)}.}
  \item{...}{Additional arguments.}
}

\details{
Annotated \code{telemetry} objects are required for \code{color} \code{by} arguments \code{"sun"}, \code{"moon"}, \code{"season"}, or \code{"tropic"}.

\code{by="time"} colors tracking data with a gradient that increases in time.
\code{by="sun"} colors according to the sine of the sun's altitude, which is proportional to solar flux during daylight hours.
\code{by="moon"} colors according to the illuminated fraction of the moon.
\code{by="season"} colors according to the length of the day, and therefore corresponds to the local season.
\code{by="tropic"} currently colors according to the calender day, but will eventually be upgraded to tropical-year cycle. The default \code{col.fn} argument runs from blue to red with increasing time, sunlight, moonlight, or day length.

\code{by="individual"} assigns colors to minimize the maximum combined spatial and color overlap.
Finding the best color assignment is an \eqn{NP}-hard problem that is here approximated in \eqn{O(N^3)} time with a custom greedy algorithm.

Other named columns in the \code{telemetry} object can also be used with \code{color}, by specifying the column name with \code{by}.
}


\value{\code{annotate} returns an annotated telemetry object with extra columns to facilitate coloring.
\code{color} returns a valid \code{col} argument for \code{{plot.telemetry}}.
}

%\references{}

\author{ C. H. Fleming. }

%\note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{plot.telemetry}} }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)

# assign distinct colors to buffalo
COL <- color(buffalo,by='individual')
# Notice the separation into RGB and CMY for maximum contrast
plot(buffalo,col=COL)

# annotate buffalo with sunlight data and compute colors
buffalo <- annotate(buffalo,cores=2) # CRAN policy limits to 2 cores
COL <- color(buffalo,by='sun')

# use North-preserving projection and plot
projection(buffalo) <- median(buffalo)
plot(buffalo,col=COL)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
