\name{jaguar}
\alias{jaguar}
\encoding{UTF-8}
\docType{data}

\title{Jaguar data from the Jaguar movement database.}

\description{\code{x-y} projected GPS data on 4 jaguar.
Please cite Morato et al (2018) when publishing with these data.}

\usage{data("jaguar")}

\format{ A list of 4 \code{telemetry} objects.}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}

\references{
R. G. Morato et al,
``Jaguar movement database: a GPS‐based movement dataset of an apex predator in the Neotropic'',
\href{https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2379}{Ecology, 99:7, 1691-1691 (2018)}.
}

\seealso{ \code{\link{as.telemetry}}, \code{\link{plot.telemetry}}, \code{\link{buffalo}}, \code{\link{coati}}, \code{\link{gazelle}}, \code{\link{pelican}}, \code{\link{turtle}}, \code{\link{wolf}}. }

\examples{
# Load package and data
library(ctmm)
data("jaguar")

# Plot all jaguar locations
plot(jaguar,col=rainbow(length(jaguar)))
}

\keyword{datasets}
