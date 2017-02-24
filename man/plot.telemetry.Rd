\name{plot.telemetry}
\alias{plot}
\alias{plot.telemetry}
\alias{zoom,list-method}
\alias{zoom,telemetry-method}
\alias{zoom,UD-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting methods for telemetry objects.
}
\description{
Produces simple plots of \code{telemetry} objects, possibly overlayed with a Gaussian \code{ctmm} movement model or a \code{UD} utilization distribution.
}
\usage{
plot(x,y,...)

\S3method{plot}{telemetry}(x, CTMM=NULL, UD=NULL, level.UD=0.95, level=0.95, DF="CDF", col="red",
col.level="black", col.DF="blue", col.grid="grey", pch=1, fraction=1, add=FALSE,
xlim=NULL, ylim=NULL, cex=1, lwd=1, ...)

\S4method{zoom}{list}(x,...)

\S4method{zoom}{telemetry}(x,fraction=1,...)

\S4method{zoom}{UD}(x,fraction=1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{telemetry} or \code{UD} object.}
  \item{y}{Unused option.}
  \item{CTMM}{Optional Gaussian \code{ctmm} movement model from the output of \code{ctmm.fit} or list of such objects.}
  \item{UD}{Optional \code{UD} object such as from the output of \code{akde} or list of such objects.}
  \item{level.UD}{Confidence level of Gaussian \code{ctmm} model or \code{UD} estimate contours to be displayed. I.e., \code{level.UD=0.50} can yield the 50\% core home range within the rendered contours.}
  \item{level}{ Confidence levels placed on the contour estimates themselves. I.e., the above 50\% core home-range area can be estimated with 95\% confidence via \code{level=0.95}. }
  \item{DF}{Plot the maximum likelihood probability density function \code{"PDF"} or cumulative distribution function \code{"CDF"}.}
  \item{col}{Color option for telemetry data. Can be an array or list of arrays.}
  \item{col.level}{Color option for home-range contours. Can be an array.}
  \item{col.DF}{Color option for the density function. Can be an array.}
  \item{col.grid}{Color option for the maximum likelihood \code{akde} bandwidth grid.}
  \item{pch}{Plotting symbol. Can be an array or list of arrays.}
  \item{fraction}{Fraction of the data, Gaussian \code{ctmm}, or \code{UD} range to plot, whichever is larger.}
  \item{add}{Setting to \code{TRUE} will disable the unit conversions and base layer plot, so that \code{plot.telemetry} can be overlayed atop other outputs more easily.}
  \item{xlim}{The \code{x} limits \code{c(x1, x2)} of the plot.}
  \item{ylim}{The \code{y} limits \code{c(y1, y2)} of the plot.}
  \item{cex}{Relative size of plotting symbols. Only used when errors are missing.}
  \item{lwd}{Line widths of \code{UD} contours.}
  \item{\dots}{ Additional options passed to \code{plot}.}
}


\details{
Confidence intervals placed on the \code{ctmm} Gaussian home-range contour estimates only represent uncertainty in the area's magnitude and not uncertainty in the mean location, eccentricity, or orientation angle. For \code{akde} \code{UD} estimates, the provided contours also only represent uncertainty in the magnitude of the area. With \code{akde} estimates, it is also important to note the scale of the bandwidth and, by default, grid cells are plotted with \code{akde} contours such that their length and width matches that of a bandwidth kernels' standard deviation in each direction. Therefore, this grid provides a visual approximation of the kernel-density estimate's ``resolution''.}

\value{Returns a plot of \eqn{x} vs. \eqn{y}, and, if specified, Gaussian \code{ctmm} distribution or \code{UD}.
\code{akde} \code{UD} plots also come with a standard resolution grid.
\code{zoom} includes a zoom slider to manipulate \code{fraction}.}

\references{
C. H. Fleming, J. M. Calabrese, T. Mueller, K.A. Olson, P. Leimgruber, and W. F. Fagan. (2014). From fine-scale foraging to home ranges: A semi-variance approach to identifying movement modes across spatiotemporal scales. \href{http://www.jstor.org/discover/10.1086/675504}{The American Naturalist, 183(5), E154-E167.}

C. H. Fleming and J. M. Calabrese and T. Mueller and K. A. Olson and P. Leimgruber and W. F. Fagan (2014).
Non-Markovian maximum likelihood estimation of autocorrelated movement processes
\href{http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12176/abstract}{Methods in Ecology and Evolution, 5(5) 462-472}.

C. H. Fleming and W. F. Fagan and T. Mueller and K. A. Olson and P. Leimgruber and J. M. Calabrese (2015).
Rigorous home-range estimation with movement data: A new autocorrelated kernel-density estimator.
\href{http://www.esajournals.org/doi/abs/10.1890/14-2010.1}{Ecology, DOI:10.1890/14-2010.1}.
}

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
