\name{gazelle}
\alias{gazelle}
\encoding{UTF-8}
\docType{data}

\title{Mongolian gazelle GPS dataset from the Mongolia's Eastern Steppe.}

\description{\code{x-y} projected GPS data on 36 Mongolian gazelle.}

\usage{data("gazelle")}

\format{ A list of 36 \code{telemetry} objects.}

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
C. H. Fleming, J. M. Calabrese, T. Mueller, K.A. Olson, P. Leimgruber, and W. F. Fagan. Data from: From fine-scale foraging to home ranges: A semi-variance approach to identifying movement modes across spatiotemporal scales. Dryad Digital Repository (2014) \doi{10.5061/dryad.45157}.
}

\seealso{ \code{\link{as.telemetry}}, \code{\link{plot.telemetry}}, \code{\link{buffalo}}, \code{\link{coati}}, \code{\link{jaguar}}, \code{\link{pelican}}, \code{\link{turtle}}, \code{\link{wolf}}. }

\examples{
# Load package and data
library(ctmm)
data("gazelle")

# Plot a gazelle's locations
plot(gazelle[[18]])
}

\keyword{datasets}
