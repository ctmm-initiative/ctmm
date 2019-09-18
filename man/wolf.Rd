\name{wolf}
\alias{wolf}
\encoding{UTF-8}
\docType{data}

\title{Maned wolf GPS dataset from The Maned Wolf Conservation Program.}

\description{\code{x-y} projected GPS data on 8 Maned wolves.
Please contact Rogerio Cunha de Paula (rogercunha@gmail.com) if you want to publish with these data.}

\usage{data("wolf")}

\format{ A list of 8 \code{telemetry} objects.}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

%\references{}

\seealso{ \code{\link{as.telemetry}}, \code{\link{plot.telemetry}}, \code{\link{buffalo}}, \code{\link{coati}}, \code{\link{gazelle}}, \code{\link{pelican}}, \code{\link{turtle}}. }

\examples{
# Load package and data
library(ctmm)
data("wolf")

# Plot a wolf's locations
plot(wolf[[8]])
}

\keyword{datasets}
