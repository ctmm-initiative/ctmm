\name{export}
\alias{export}
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
\description{Functions to export \code{ctmm} data formats into common \code{sp}, \code{raster}, and ESRI formats.}

\usage{
\S4method{raster}{UD}(x,DF="CDF",...)

SpatialPoints.telemetry(object,...)

SpatialPointsDataFrame.telemetry(object,...)

SpatialPolygonsDataFrame.telemetry(object,level.UD=0.95,...)

SpatialPolygonsDataFrame.UD(object,level.UD=0.95,level=0.95,...)

\S4method{writeRaster}{UD,character}(x,filename,format,DF="CDF",...)

writeShapefile(object,folder,file=NULL,...)

\S3method{writeShapefile}{telemetry}(object,folder,file=NULL,error=TRUE,level.UD=0.95,...)

\S3method{writeShapefile}{UD}(object,folder,file=NULL,level.UD=0.95,level=0.95,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ \code{UD} object. }
  \item{DF}{Rasterize the probability density function \code{"PDF"}, probability mass function \code{"PMF"}, or cumulative distribution function \code{"CDF"}.}
  \item{object}{ \code{telemetry} or \code{UD} object. }
  %\item{UERE}{Device-dependent telemetry error in meters. Only necessary for uncalibrated data.}
  \item{level.UD}{ Confidence level of the UD area. I.e., the 50\% core home range would be given by \code{level.UD=0.50}. }
  \item{level}{ Confidence level for the magnitude of the above area. I.e., the 95\% CI of the core home range area. }
  \item{filename}{Character name of file for raster file.}
  \item{format}{Character format, if not inferred from \code{filename} extension (see \code{\link[raster]{writeRaster}}.)}
  \item{folder}{Character name of folder for shapefile.}
  \item{file}{Character name of files for shapefile.}
  \item{error}{Export error circles/ellipses as polygons if \code{TRUE}.}
  \item{...}{ Optional arguments passed to \code{\link[raster]{writeRaster}}, \code{\link[rgdal]{writeOGR}}, etc..}
}

\details{\code{writeRaster} writes a raster file to disk, with pixel values corresponding to the density function.

\code{writeShapefile} writes a shapefile to disk, with UD polygons corresponding to the low, ML, and high home-range area estimates.}

\value{\code{raster.UD} returns a \code{raster} of the maximum likelihood (ML) distribution function \code{DF}.
\code{DF=PDF} gives the average probability density per cell,
\code{DF=PMF} gives the total probability per cell,
and \code{DF=CDF} gives the cumulative probability.

\code{SpatialPoints.telemetry} returns a single \code{spatialPoints} object for the \code{x}-\code{y} locations, with individual identity information lost.

\code{SpatialPointsDataFrame.telemetry} returns a \code{SpatialPointsDataFrame} with the individual identities recorded in the data frame.

\code{SpatialPolygonsDataFrame.telemetry} returns a \code{SpatialPolygonsDataFrame} that encodes the location estimate's error circles/ellipses.

\code{SpatialPolygonsDataFrame.UD} returns a \code{SpatialPolygonsDataFrame} of the low, ML, and high home-range area estimates, in the appropriate order for plotting.
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
