\name{as.telemetry}
\alias{as.telemetry}
\alias{as.telemetry.character}
\alias{as.telemetry.data.frame}
\alias{as.telemetry.Move}
\alias{summary.telemetry}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Import, coerce, and summarize MoveBank data
}
\description{Functions to import MoveBank csv files, \code{data.frame}, and \code{Move} objects, coerce them into \code{telemetry} objects, and summarize them.}

\usage{
as.telemetry(object,timeformat="",timezone="UTC",projection=NULL,UERE=NULL,...)

\S3method{as.telemetry}{character}(object,timeformat="",timezone="UTC",projection=NULL,UERE=NULL,...)

\S3method{as.telemetry}{data.frame}(object,timeformat="",timezone="UTC",projection=NULL,UERE=NULL,...)

\S3method{as.telemetry}{Move}(object,timeformat="",timezone="UTC",projection=NULL,UERE=NULL,...)

\S3method{summary}{telemetry}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A MoveBank CSV filename, MoveBank \code{data.frame} object, or \code{Move} object to coarce, or a \code{telemetry} object to summarize.}
  \item{timeformat}{Format argument for \code{\link{strptime}}.}
  \item{timezone}{Timezone argument for \code{\link{strptime}}.}
  \item{projection}{Optional PROJ.4 projection argument to be fed to \code{\link[rgdal]{project}}.}
  \item{UERE}{"User equivalent range errors" needed to translate DOP values into location uncertainty.}
  \item{...}{Optional arguments to be fed to \code{\link{read.csv}}.}
}

\details{
If no projection argument is specified, a two-point equidistant projection is calculated that should be good for most range resident and migratory species.
Global migrations that are not along one geodesic (locally straight line) will probably suffer distortion.

\code{as.telemetry} assumes \href{http://rstudio.com}{Movebank naming conventions}.
Sufficient MoveBank columns include \code{individual.local.identifier} (or \code{tag.local.identifier}), \code{timestamp}, \code{location.long} and \code{location.lat}, while the optional Movebank columns are \code{eobs.horizontal.accuracy.estimate} and \code{GPS.HDOP}.
In the future, we will be supporting more GPS and ARGOS columns, including speed, heading, further error classifications, etc.

The \code{GPS.HDOP} value requires the device's \code{UERE} value to be translated into telemetry errors. The \code{UERE} represents the RMS error given ideal satellite conditions. Therefore, if your device is purported to have an accuracy of 2-10 meters, then the \code{UERE} is likely around 2 meters. Consult your device manual and manufacturer for a specific value or see \code{\link{uere}} for UERE calibration.
}

\value{\code{as.telemetry} returns a single \code{telemetry} object or list of \code{telemetry} objects if multiple animals are identified.}

%\references{}

\author{C. H. Fleming, B. Kranstauber, G. Péron, and K. Safi.}

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{plot.telemetry}}, \code{\link{SpatialPoints.telemetry}}, \code{\link{uere}}.  }

%\examples{
%# Load package and data
%library(ctmm)
%data(buffalo)

%# Extract movement data for a single animal
%billy <- buffalo[[1]]

%# Plot the data
%plot.telemetry(billy)
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
