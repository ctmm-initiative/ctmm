\name{export}
\alias{export}
\alias{as.sf}
\alias{raster,UD-method}
\alias{SpatialPoints.telemetry}
\alias{SpatialPointsDataFrame.telemetry}
\alias{SpatialPolygonsDataFrame.telemetry}
\alias{SpatialPolygonsDataFrame.UD}
\alias{writeRaster,UD,character-method}
\alias{writeShapefile}
\alias{writeShapefile.telemetry}
\alias{writeShapefile.UD}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Export ctmm data formats}
\description{Functions to export \code{ctmm} data formats into common \code{sp}, \code{sf}, \code{raster}, and ESRI formats.}

\usage{
as.sf(x,error=FALSE,...)

SpatialPoints.telemetry(object,...)

SpatialPointsDataFrame.telemetry(object,...)

SpatialPolygonsDataFrame.telemetry(object,level.UD=0.95,...)

SpatialPolygonsDataFrame.UD(object,level.UD=0.95,level=0.95,...)

\S4method{raster}{UD}(x,DF="CDF",...)

writeShapefile(object,folder,file=NULL,...)

\S3method{writeShapefile}{telemetry}(object,folder,file=NULL,error=TRUE,level.UD=0.95,...)

\S3method{writeShapefile}{UD}(object,folder,file=NULL,level.UD=0.95,level=0.95,...)

\S4method{writeRaster}{UD,character}(x,filename,format,DF="CDF",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ \code{telemetry} or \code{UD} object. }
  \item{error}{Export \code{telemetry} location error circles/ellipses as polygons if \code{TRUE}.}
  \item{object}{ \code{telemetry} or \code{UD} object. }
  \item{level.UD}{ Coverage level of the UD area. I.e., the 50\% core home range would be given by \code{level.UD=0.50}. }
  \item{level}{ Confidence level for the magnitude of the above area. I.e., the 95\% CI of the core home range area. }
  \item{DF}{Rasterize the probability density function \code{"PDF"}, probability mass function \code{"PMF"}, or cumulative distribution function \code{"CDF"}.}
  \item{folder}{Character name of folder for shapefile.}
  \item{file}{Character name of files for shapefile.}
  \item{filename}{Character name of file for raster file.}
  \item{format}{Character format, if not inferred from \code{filename} extension (see \code{\link[raster]{writeRaster}}).}
  \item{...}{ Optional arguments passed to \code{\link[raster]{writeRaster}}, \code{\link[rgdal]{writeOGR}}, etc..}
}

\details{\code{as.sf} exports \code{ctmm} objects to the \code{sf} format. Arguments to \code{ctmm} \code{Spatial}* export functions can also be used, such as \code{level.UD} and \code{level}.

\code{Spatial}* functions export \code{ctmm} objects to \code{sp} formats.

\code{raster} exports \code{UD} object point-estimates distribution functions (\code{DF}) to \code{raster} objects.
\code{DF="PDF"} gives the average probability density per cell,
\code{DF="PMF"} gives the total probability per cell,
and \code{DF="CDF"} gives the cumulative probability.

\code{writeShapefile} writes a shapefile to disk, with UD polygons corresponding to the low-CI, point-estimate, and high-CI home-range area estimates.

\code{writeRaster} writes a raster file to disk, with pixel values corresponding to the distribution function \code{DF}.}

\value{\code{as.sf} returns an \code{sf} object for the input points or polygons, with individual identity and other information retained.

\code{SpatialPoints.telemetry} returns a single \code{spatialPoints} object for the \code{x}-\code{y} locations, without individual identity and other information retained.

\code{SpatialPointsDataFrame.telemetry} returns a \code{SpatialPointsDataFrame} with the individual identities and other data recorded in the data frame retained.

\code{SpatialPolygonsDataFrame.telemetry} returns a \code{SpatialPolygonsDataFrame} that encodes the location estimate's error circles/ellipses.

\code{SpatialPolygonsDataFrame.UD} returns a \code{SpatialPolygonsDataFrame} of the low-CI, point-estimate, and high-CI home-range area estimates, in the appropriate order for plotting.

\code{raster} returns a \code{raster} of the point-estimate distribution function \code{DF}, given a \code{UD} object.
}

%\references{}

\author{C. H. Fleming and K. Safi.}

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde}}, \code{\link{as.telemetry}}, \code{\link{occurrence}}. }

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
