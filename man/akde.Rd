\name{akde}
\alias{akde}
\alias{akde.telemetry}
\alias{akde.list}
\alias{mean.UD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate an autocorrelated kernel density estimate}
\description{
This function calculates autocorrelated kernel density home-range estimates from \code{telemetry} data and a corresponding continuous-time movement model.
}
\usage{
akde(data,CTMM,VMM=NULL,debias=TRUE,smooth=TRUE,error=0.001,res=10,grid=NULL,...)

\S3method{akde}{telemetry}(data,CTMM,VMM=NULL,debias=TRUE,smooth=TRUE,error=0.001,res=10,grid=NULL,...)

\S3method{akde}{list}(data,CTMM,VMM=NULL,debias=TRUE,smooth=TRUE,error=0.001,res=10,grid=NULL,...)

\S3method{mean}{UD}(x,...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ 2D timeseries telemetry data represented as a \code{telemetry} object or list of objects. }
  \item{CTMM}{A \code{ctmm} movement model from the output of \code{ctmm.fit} or list of objects.}
  \item{VMM}{An optional vertical \code{ctmm} object for 3D home-range calculation.}
  \item{debias}{Debias the distribution for area estimation (AKDEc).}
  \item{smooth}{"Smooth" out errors from the data.}
  \item{error}{Target probability error.}
  \item{res}{Number of grid points along each axis, relative to the bandwidth.}
  \item{grid}{Optional grid specification with columns labeled \code{x} and \code{y}. Not yet supported.}
  \item{...}{Arguments passed to all instances of \code{\link{bandwidth}}, such as \code{weights}.}
  \item{x}{A list of \code{UD}s calculated on the same grid.}
}

%\details{}

\value{Returns a \code{UD} object: a list with the sampled grid line locations \code{r$x} and \code{r$y}, the extent of each grid cell \code{dr}, the probability density and cumulative distribution functions evaluated on the sampled grid locations \code{PDF} & \code{CDF}, the optimal bandwidth matrix \code{H}, and the effective sample size of the data in \code{DOF.H}.

For weighted AKDE, please note additional \code{...} arguments passed to \code{\link{bandwidth}} and the \code{weights=TRUE} argument, specifically.

When feeding in lists of \code{telemetry} and \code{ctmm} objects, all UDs will be calculated on the same grid. These UDs can be averaged with the \code{mean} command, however this is not an optimal way to calculate population ranges.
}

\references{
C. H. Fleming, W. F. Fagan, T. Mueller, K. A. Olson, P. Leimgruber, J. M. Calabrese.
Rigorous home-range estimation with movement data: A new autocorrelated kernel-density estimator.
\href{http://www.esajournals.org/doi/abs/10.1890/14-2010.1}{Ecology, 96:5, 1182-1188 (2015)}.

C. H. Fleming, J. M. Calabrese.
A new kernel-density estimator for accurate home-range and species-range area estimation.
\href{http://onlinelibrary.wiley.com/wol1/doi/10.1111/2041-210X.12673/abstract}{Methods in Ecology and Evolution, 8:5, 571-579 (2016)}.

%D. Marcotte. (1996). Fast variogram computation with FFT. \emph{Computers and Geosciences} 22(10), 1175-1186.
}

\author{ C. H. Fleming and K. Winner. }

\note{ In the case of coarse grids, the value of \code{PDF} in a grid cell corresponds to the average probability density over the entire rectangular cell.

Prior to \code{ctmm} v0.3.2, the default AKDE method was the autocorrelated Gaussian reference function bandwidth.
Starting in v0.3.2, the default AKDE method is the autocorrelated Gaussian reference function bandwidth with debiased area.

Prior to \code{ctmm} v0.3.1, AKDEs included only errors due to autocorrelation uncertainty, which are insignificant in cases such as IID data.
Starting in v0.3.1, \code{akde} calculated an effective sample size \code{DOF.H} and used this to estimate area uncertainty under a Gaussian reference function approxmation.
In v0.3.2, this method was further improved to use \code{DOF.area} from the Gaussian reference function approximation.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{bandwidth}}, \code{\link{raster,UD-method}} }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)
Cilla <- buffalo$Cilla

GUESS <- ctmm.guess(Cilla,interactive=FALSE)
FIT <- ctmm.fit(Cilla,GUESS)

# Compute akde object
UD <- akde(Cilla,FIT)

# Plot data with AKDE contours
plot(Cilla,UD=UD)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
