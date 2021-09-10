\name{akde}
\alias{akde}
\alias{akde.telemetry}
\alias{akde.list}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate an autocorrelated kernel density estimate}
\description{
This function calculates autocorrelated kernel density home-range estimates from \code{telemetry} data and a corresponding continuous-time movement model.
}
\usage{
akde(data,CTMM,VMM=NULL,variable="utilization",debias=TRUE,weights=FALSE,smooth=TRUE,
     error=0.001,res=10,grid=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ 2D timeseries telemetry data represented as a \code{telemetry} object or list of objects. }
  \item{CTMM}{A \code{ctmm} movement model from the output of \code{ctmm.fit} or list of objects.}
  \item{VMM}{An optional vertical \code{ctmm} object for 3D home-range calculation.}
  \item{variable}{Not yet supported.}
  \item{debias}{Debias the distribution for area estimation (AKDEc).}
  \item{smooth}{"Smooth" out errors from the data.}
  \item{weights}{Optimally weight the data to account for temporal sampling bias (See \code{\link{bandwidth}} for details).}
  \item{error}{Target probability error.}
  \item{res}{Number of grid points along each axis, relative to the bandwidth.}
  \item{grid}{Optional grid specification via \code{raster}, \code{UD}, or list of arguments (See `Details' below).}
  \item{...}{Arguments passed to all instances of \code{\link{bandwidth}}.}
}

\details{For weighted AKDE, please note additional \code{...} arguments passed to \code{\link{bandwidth}}, which can have a large impact on computation time in certain cases.

When feeding in lists of \code{telemetry} and \code{ctmm} objects, all UDs will be calculated on the same grid. These UDs can be averaged with the \code{\link{mean.UD}} command.

If a \code{UD} or \code{raster} object is supplied in the \code{grid} argument, then the estimate will be calculated on the same grid. Alternatively, a list of grid arguments can be supplied, with any of the following components:
\describe{
\item{\code{r}}{A list with vectors \code{x} and \code{y} that define the grid-cell midpoints.}
\item{\code{dr}}{A vector setting the \code{x} and \code{y} cell widths in meters. Equivalent to \code{\link{res}} for \code{raster} objects.}
\item{\code{extent}}{The \eqn{x}-\eqn{y} extent of the grid cells, formatted as from the output of \code{\link{extent}}.}
\item{\code{align.to.origin}}{Logical value indicating that cell midpoint locations are aligned to be an integer number of \code{dr} steps from the projection origin.}
}
}

\value{Returns a \code{UD} object: a list with the sampled grid line locations \code{r$x} and \code{r$y}, the extent of each grid cell \code{dr}, the probability density and cumulative distribution functions evaluated on the sampled grid locations \code{PDF} & \code{CDF}, the optimal bandwidth matrix \code{H}, and the effective sample size of the data in \code{DOF.H}.
}

\references{
C. H. Fleming, W. F. Fagan, T. Mueller, K. A. Olson, P. Leimgruber, J. M. Calabrese,
``Rigorous home-range estimation with movement data: A new autocorrelated kernel-density estimator'',
Ecology, 96:5, 1182-1188 (2015) \doi{10.1890/14-2010.1}.

C. H. Fleming, J. M. Calabrese,
``A new kernel-density estimator for accurate home-range and species-range area estimation'',
Methods in Ecology and Evolution, 8:5, 571-579 (2017) \doi{10.1111/2041-210X.12673}.

C. H. Fleming, D. Sheldon, W. F. Fagan, P. Leimgruber, T. Mueller, D. Nandintsetseg, M. J. Noonan, K. A. Olson, E. Setyawan, A. Sianipar, J. M. Calabrese,
``Correcting for missing and irregular data in home-range estimation'',
Ecological Applications, 28:4, 1003-1010 (2018) \doi{10.1002/eap.1704}.

%D. Marcotte. (1996). Fast variogram computation with FFT. \emph{Computers and Geosciences} 22(10), 1175-1186.
}

\author{ C. H. Fleming and K. Winner. }

\note{ In the case of coarse grids, the value of \code{PDF} in a grid cell corresponds to the average probability density over the entire rectangular cell.

The \code{PDF} estimate is not re-normalized to 1, and may fall short of this by the target numerical \code{error}. If inspecting quantiles that are very far from the data, the quantiles may hit the grid boundary or become erratic, making it necessary to reduce the numerical \code{error} target. However, default arguments should be able to render any quantiles of reasonable accuracy.

Prior to \code{ctmm} v0.3.2, the default AKDE method was the autocorrelated Gaussian reference function bandwidth.
Starting in v0.3.2, the default AKDE method is the autocorrelated Gaussian reference function bandwidth with debiased area.

Prior to \code{ctmm} v0.3.1, AKDEs included only errors due to autocorrelation uncertainty, which are insignificant in cases such as IID data.
Starting in v0.3.1, \code{akde} calculated an effective sample size \code{DOF.H} and used this to estimate area uncertainty under a Gaussian reference function approxmation.
In v0.3.2, this method was further improved to use \code{DOF.area} from the Gaussian reference function approximation.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{bandwidth}}, \code{\link{mean.UD}}, \code{\link{raster,UD-method}}, \code{\link{revisitation}} }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)
DATA <- buffalo$Cilla

# calculate fit guess object
GUESS <- ctmm.guess(DATA,interactive=FALSE)
# in general, you should be running ctmm.select here instead of ctmm.fit
FIT <- ctmm.fit(DATA,GUESS)

# Compute akde object
UD <- akde(DATA,FIT)

# Plot data with AKDE
plot(DATA,UD=UD)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
