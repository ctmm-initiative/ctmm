\name{pelican}
\alias{pelican}
\encoding{UTF-8}
\docType{data}

\title{Brown Pelican GPS and ARGOS data.}

\description{GPS and ARGOS data on a single brown pelican (Pelecanus occidentalis).
Please contact Autumn-Lynn Harrison (HarrisonAL@si.edu) if you want to publish with these data.

Funding for Brown Pelican tracking was provided by the Friends of the National Zoo Conservation Research Grant and ConocoPhillips Global Signature Program. Field support provided by D. Brinker.}

\usage{data("pelican")}

\format{ A list of 2 \code{telemetry} objects.}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{}

\seealso{ \code{\link{as.telemetry}}, \code{\link{plot.telemetry}}, \code{\link{buffalo}}, \code{\link{coati}}, \code{\link{gazelle}}, \code{\link{jaguar}}, \code{\link{turtle}}, \code{\link{wolf}}. }

\examples{
# Load package and data
library(ctmm)
data("pelican")
names(pelican)

# Plot all sampled locations
plot(pelican,col=c('blue','red'))
}

\keyword{datasets}
