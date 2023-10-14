\name{plot.telemetry}
\alias{plot}
\alias{plot.telemetry}
\alias{zoom,list-method}
\alias{zoom,telemetry-method}
\alias{zoom,UD-method}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting methods for telemetry objects}
\description{
Produces simple plots of \code{telemetry} objects, possibly overlayed with a Gaussian \code{ctmm} movement model or a \code{UD} utilization distribution.
}
\usage{
plot(x,y,...)

\S3method{plot}{telemetry}(x,CTMM=NULL,UD=NULL,col.bg="white",cex=NULL,col="red",lwd=1,pch=1,type='p',
     error=TRUE,transparency.error=0.25,velocity=FALSE,DF="CDF",col.UD="blue",
     col.grid="white",labels=NULL,convex=FALSE,level=0.95,level.UD=0.95,col.level="black",
     lwd.level=1,SP=NULL,border.SP=TRUE,col.SP=NA,R=NULL,col.R="green",legend=FALSE,
     fraction=1,xlim=NULL,ylim=NULL,ext=NULL,units=TRUE,add=FALSE,...)

\S4method{zoom}{list}(x,...)

\S4method{zoom}{telemetry}(x,fraction=1,...)

\S4method{zoom}{UD}(x,fraction=1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{telemetry}, \code{ctmm} or \code{UD} object, or list of such objects.}
  \item{y}{Unused option.}
  \item{CTMM}{Optional Gaussian \code{ctmm} movement model from the output of \code{ctmm.fit} or list of such objects.}
  \item{UD}{Optional \code{UD} object such as from the output of \code{akde} or list of such objects.}
  \item{col.bg}{Background color}
  \item{cex}{Relative size of plotting symbols. Only used when \code{error=FALSE}, because \code{error=TRUE} uses the location-error radius instead of \code{cex}.}
  \item{col}{Color option for telemetry data. Can be an array or list of arrays.}
  \item{lwd}{Line widths of \code{telemetry} points.}
  \item{pch}{Plotting symbol. Can be an array or list of arrays.}
  \item{type}{How plot points are connected. Can be an array.}
  \item{error}{Plot error circles/ellipses if present in the data. \code{error=2} will fill in the circles and \code{error=3} will plot densities instead. \code{error=FALSE} will disable this feature.}
  \item{transparency.error}{Transparency scaling for erroneous locations when \code{error=1:2}. \code{trans=0} disables transparancy. Should be no greater than \code{1}.}
  \item{velocity}{Plot velocity vectors if present in the data.}
  \item{DF}{Plot the maximum likelihood probability density function \code{"PDF"} or cumulative distribution function \code{"CDF"}.}
  \item{col.UD}{Color option for the density function. Can be an array.}
  \item{col.grid}{Color option for the maximum likelihood \code{akde} bandwidth grid. \code{col.grid=NA} will disable the plotting of the bandwidth grid.}
  \item{labels}{Labels for UD contours. Can be an array or list of arrays.}
  \item{convex}{Plot convex coverage-area contours if \code{TRUE}. By default, the highest density region (HDR) contours are plotted.}
  \item{level}{ Confidence levels placed on the contour estimates themselves. I.e., the above 50\% core home-range area can be estimated with 95\% confidence via \code{level=0.95}. \code{level=NA} will disable the plotting of confidence intervals.}
  \item{level.UD}{Coverage level of Gaussian \code{ctmm} model or \code{UD} estimate contours to be displayed. I.e., \code{level.UD=0.50} can yield the 50\% core home range within the rendered contours.}
  \item{col.level}{Color option for home-range contours. Can be an array.}
  \item{lwd.level}{Line widths of \code{UD} contours.}
  \item{SP}{\code{SpatialPolygonsDataFrame} object for plotting a shapefile base layer.}
  \item{border.SP}{Color option for shapefile polygon boundaries.}
  \item{col.SP}{Color option for shapefile polygon regions.}
  \item{R}{Background raster, such as habitat \code{\link{suitability}}.}
  \item{col.R}{Color option for background raster.}
  \item{legend}{Plot a color legend for background raster.}
  \item{fraction}{Quantile fraction of the data, Gaussian \code{ctmm}, or \code{UD} range to plot, whichever is larger.}
  \item{xlim}{The \code{x} limits \code{c(x1, x2)} of the plot (in SI units).}
  \item{ylim}{The \code{y} limits \code{c(y1, y2)} of the plot (in SI units).}
  \item{ext}{Plot extent alternative to \code{xlim} and \code{ylim} (see \code{\link{extent}}).}
  \item{units}{Convert axes to natural units.}
  \item{add}{Setting to \code{TRUE} will disable the unit conversions and base layer plot, so that \code{plot.telemetry} can be overlayed atop other outputs more easily.}
  \item{\dots}{ Additional options passed to \code{plot}.}
}


\details{
Confidence intervals placed on the \code{ctmm} Gaussian home-range contour estimates only represent uncertainty in the area's magnitude and not uncertainty in the mean location, eccentricity, or orientation angle. For \code{akde} \code{UD} estimates, the provided contours also only represent uncertainty in the magnitude of the area. With \code{akde} estimates, it is also important to note the scale of the bandwidth and, by default, grid cells are plotted with \code{akde} contours such that their length and width matches that of a bandwidth kernels' standard deviation in each direction. Therefore, this grid provides a visual approximation of the kernel-density estimate's ``resolution''. Grid lines can be disabled with the argument \code{col.grid=NA}.}

\value{Returns a plot of \eqn{x} vs. \eqn{y}, and, if specified, Gaussian \code{ctmm} distribution or \code{UD}.
\code{akde} \code{UD} plots also come with a standard resolution grid.
\code{zoom} includes a zoom slider to manipulate \code{fraction}.}

%\references{}

\author{C. H. Fleming.}

\note{If \code{xlim} or \code{ylim} are provided, then the smaller or absent range will be expanded to ensure \code{asp=1}.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde}}, \code{\link{ctmm.fit}}, \code{\link{plot}}, \code{\link{SpatialPoints.telemetry}}.  }

\examples{
# Load package and data
library(ctmm)
data(buffalo)

# Plot the data
plot(buffalo,col=rainbow(length(buffalo)))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
