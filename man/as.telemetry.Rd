\name{as.telemetry}
\alias{as.telemetry}
\alias{as.telemetry.character}
\alias{as.telemetry.data.frame}
\alias{as.telemetry.Move}
\alias{summary.telemetry}
\alias{head}
\alias{head.telemetry}
\alias{tail}
\alias{tail.telemetry}
\alias{tbind}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Import, coerce, summarize, and combine MoveBank data}
\description{Functions to import MoveBank csv files, \code{data.frame}, and \code{Move} objects, coerce them into \code{telemetry} objects, summarize them, and combine data from multiple tracking devices.}

\usage{
as.telemetry(object,timeformat="auto",timezone="UTC",projection=NULL,datum="WGS84",
             dt.hot=NA,timeout=Inf,na.rm="row",mark.rm=FALSE,keep=FALSE,drop=TRUE,...)

\S3method{as.telemetry}{character}(object,timeformat="auto",timezone="UTC",projection=NULL,datum="WGS84",
             dt.hot=NA,timeout=Inf,na.rm="row",mark.rm=FALSE,keep=FALSE,drop=TRUE,...)

\S3method{as.telemetry}{data.frame}(object,timeformat="auto",timezone="UTC",projection=NULL,datum="WGS84",
             dt.hot=NA,timeout=Inf,na.rm="row",mark.rm=FALSE,keep=FALSE,drop=TRUE,...)

\S3method{as.telemetry}{Move}(object,timeformat="auto",timezone="UTC",projection=NULL,datum="WGS84",
             dt.hot=NA,timeout=Inf,na.rm="row",mark.rm=FALSE,keep=FALSE,drop=TRUE,...)

\S3method{summary}{telemetry}(object,...)

head(x,...)

\S3method{head}{telemetry}(x,n=6L,...)

tail(x,...)

\S3method{tail}{telemetry}(x,n=6L,...)

tbind(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A MoveBank CSV filename, MoveBank \code{data.frame} object, or \code{Move} object to coerce, or a \code{telemetry} object to summarize.}
  \item{timeformat}{Format argument for \code{\link{strptime}}, corresponding to the input data. Alternatively \code{timeformat="auto"} will attempt to infer the timestamp format with \code{\link[parsedate]{parse_date}}.}
  \item{timezone}{Timezone argument for \code{\link{strptime}}, corresponding to the input data.}
  \item{projection}{Optional PROJ projection argument for the \bold{output} telemetry object.}
  \item{datum}{Optional argument to specify the \bold{input} longitude-latitude or UTM datum. The default is WGS84.}
  \item{dt.hot}{Time-interval threshold at which GPS location fixes can be considered as ``hot'' and location estimate precisions may be smaller (regardless of DOP value) for assigning \code{"hot"} and \code{"cold"} location classes.}
  \item{timeout}{GPS location fix timeout value (seconds) for assigning a \code{"timed-out"} location class.}
  %\item{class.2D}{Assign location classes \code{"2D"} and \code{"3D"} to fixes with \code{3} and \code{>3} satellites.}
  \item{na.rm}{If some values are \code{NA} in the data frame, are the rows (times) deleted or are the columns (data types) deleted.}
  \item{mark.rm}{Delete Movebank manually marked outliers. Also see \code{\link{outlie}}.}
  \item{keep}{Retain additonal columns after coercion. \code{keep=TRUE} retains all columns, while individual columns to retain can be specified by name.}
  \item{drop}{Only return a \code{telemetry} object for one individual if \code{TRUE}. Always return a \code{list} of telemetry objects if \code{FALSE}.}
  \item{...}{\code{telemetry} objects or a list of such objects, for \code{tbind()}. Optional arguments to be fed to \code{\link{fread}} or \code{\link{read.csv}}, in the case of compressed files, for \code{as.telemetry()}.}
  \item{x}{\code{telemetry} object.}
  \item{n}{Number of rows to return, if positive, or number of rows to omit, if negative.}
}

\details{For data that have not been corralled throuh Movebank, timestamps either need to be provided in a POSIX format (see the output of \code{Sys.time()}) or supplied with a \code{timeformat} argument for interpretation (see \code{\link{strptime}}). Alternatively, you can try your luck with \code{timeformat="auto"}, and \code{\link[parsedate]{parse_date}} will attempt to infer the format.

If no projection argument is specified, a two-point equidistant projection is calculated that should be good for most range resident and migratory species.
Global migrations that are not along one geodesic (locally straight line) will probably suffer distortion.

\code{as.telemetry()} assumes \href{https://www.movebank.org/node/2381}{Movebank naming conventions}.
Sufficient MoveBank columns include \code{individual.local.identifier} (or \code{tag.local.identifier}), \code{timestamp}, \code{location.long} and \code{location.lat}, while the optional Movebank columns include (e-obs) \code{eobs.horizontal.accuracy.estimate}, (Telonics) \code{GPS.Horizontal.Error}, \code{GPS.HDOP}, (Argos) \code{Argos.orientation}, \code{Argos.semi.minor} and \code{Argos.semi.major} or \code{Argos.location.class}, etc..
To have all columns detected and not overwrite eachother's information, \bold{it is best to have only one tracking device model per file imported}.
Multiple deployments on a single individual can be merged afterwards, using \code{tbind()}.

%The \code{GPS.HDOP} value requires the device's UERE value to be translated into telemetry errors. The UERE represents the RMS error given ideal (\eqn{HDOP=1}) satellite conditions. Therefore, if your device is purported to have an accuracy of 10 meters, then the UERE is likely around 10 meters. Consult your device manual and manufacturer for a specific value or see \code{\link{uere}} for UERE calibration.
}

\value{\code{as.telemetry} returns a single \code{telemetry} object or list of \code{telemetry} objects if multiple animals are identified.

\code{as.telemetry} will always report the smallest sampling interval, as a message, and the number repeating timestamps, as a warning. Tiny sampling intervals (and repeating timestamps) can sometimes result from misformated timestamps or an incorrect \code{timeformat} argument. However, even if geniune, they can necessitate data cleaning (\link{outlie}) or location-error modeling (\code{vignette('error')}).
}

%\references{}

\author{C. H. Fleming, X. Dong, B. Kranstauber, G. Péron, and K. Safi.}

\note{
Prior to v1.1.1, \code{datum} was required to be a full PROJ string, but starting with v1.1.1 \code{datum} is just taken to be the PROJ datum argument.
}

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
