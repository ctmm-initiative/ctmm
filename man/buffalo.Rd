\name{buffalo}
\alias{buffalo}
\docType{data}

\title{African buffalo GPS dataset from Kruger National Park, South Africa.}

\description{GPS data on six African buffalo, provided for demonstration purposes only! Contact data owner Paul Cross if you wish to do more with the data.}

\usage{data("buffalo")}

\format{ A list of 6 \code{telemetry} objects.}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{
Cross, P., Heisey, D., Bowers, J., Hay, C., Wolhuter, J., Buss, P., Hofmeyr, M., Michel, A., Bengis, R., Bird, T. et al (2009).
Disease, predation and demography: assessing the impacts of bovine tuberculosis on African buffalo by monitoring at individual and population levels.
Journal of Applied Ecology, 46, 467-475. 
}

\note{In \code{ctmm} v0.3.2 the erroneous location fix \code{606} was removed from \code{buffalo[[4]]} "Pepper".}

\seealso{ \code{\link{as.telemetry}}, \code{\link{gazelle}}, \code{\link{plot.telemetry}}. }

\examples{
# Load package and data
library(ctmm)
data("buffalo")

# Extract movement data for a single animal
cilla <- buffalo[[1]]

# Plot all sampled locations
plot(cilla)
}

\keyword{datasets}
