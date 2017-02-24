\name{variogram.fit}
\alias{ctmm.guess}
\alias{variogram.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Visually fit a movement model to a variogram
}
\description{This function plots a \code{variogram} object overlayed with a continuous-time movement model guesstimated from the variogram's shape. Sliders are given to adjust the parameter guesstimates and the result can be saved to a global variable. The intention of this function is to facilitate good starting guesses for \code{ctmm.fit}, starting with a prototype hypothesis argument \code{CTMM}, which can contain features such as \code{isotropic}, \code{range}, \code{circle}, etc..}

\usage{
ctmm.guess(data,CTMM=ctmm(),variogram=NULL,name="GUESS",interactive=TRUE)

variogram.fit(variogram,CTMM=ctmm(),name="GUESS",fraction=0.5,interactive=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{telemetry} object.}
  \item{CTMM}{Optional model prototype or initial guesstimate of the model parameters, in \code{ctmm} object format.}
  \item{name}{Name of the global variable to store the guesstimate in.}
  \item{interactive}{Boolean denoting whether to render the initial guess with interactive sliders or store the result silently.}
  \item{variogram}{ A \code{variogram} object from the output of \code{variogram}. }
  \item{fraction}{Initial fraction of the variogram to render.}
  \item{...}{Optional parameters passed to \code{plot.variogram}.}
}


\details{
By default, \code{sigma} is the asymptote of the variogram and \code{tau} is an array of autocorrelation timescales. The position timescale is roughly the time lag it takes of the variogram to reach 63\% of its asymptote. The velocity autocorrelation timescale visually corresponds to width of the concave bowl shape at the beginning of the variogram. If \code{CTMM=ctmm(range=FALSE)}, \code{sigma} is the asymptotic slope of the variogram and only the velocity timescale is finite.

By default, parameter values are estimated from the shape of the variogram. If this fails, the \code{CTMM} option can provide alternative initial guesstimates.
}

%\value{}

%\references{}

\author{ C. H. Fleming. }

\note{ If the \code{manipulate} package is unavailable, then \code{interactive} is set to \code{FALSE}. }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{ctmm.fit}}, \code{\link{plot.variogram}}, \code{\link{variogram}}. }

\examples{
#Load package and data
library(ctmm)
data(buffalo)

#Extract movement data for a single animal
cilla <- buffalo[[1]]

#Calculate variogram
SVF <- variogram(cilla)

# generate a visual fit of the variogram (requires RStudio)
variogram.fit(SVF)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
