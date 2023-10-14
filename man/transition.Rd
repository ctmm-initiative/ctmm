\name{transition}
\alias{transition}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Make a sequence of transition slide figures}
\description{This function generates a time-ordered sequence of transition slide images from a single tracking dataset.
}
\usage{
transition(data,n=3,filename="transition",height=2160,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{telemetry} object.}
  \item{n}{The desired number of slides to create.}
  \item{filename}{The base filename of the generated figures.}
  \item{height}{.}
  \item{...}{Additional arguments passed to \code{\link{plot.telemetry}}.}
}

\details{
\code{transition} partitions the tracking data into \code{n} equal-time segments, which are plotted with color on the \eqn{n^{\mathrm{th}}} slides. These are intended to be used in transition slides between \code{n} indexed presentation topics.
}

\value{Generates \code{n+1} PNG files as a side effect.}

%\references{}

\author{ C. H. Fleming. }

\note{ Currently, there is a black border that needs to be removed, such as with the LaTeX code:
\code{\\includegraphics[height=0.9\\textheight,trim={3 3 3 3},clip]{transition-1.6.png}} }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{plot.telemetry}} }

%\examples{\donttest{}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
