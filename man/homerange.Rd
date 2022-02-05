\name{homerange}
\alias{homerange}
\alias{agde}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate a range distribution estimate}
\description{
This function estimates the range distribution from \code{telemetry} data and a continuous-time movement model.
}
\usage{
homerange(data,CTMM,method="AKDE",...)

agde(CTMM,R=list(),variable="utilization",error=0.001,res=100,grid=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ 2D timeseries telemetry data represented as a \code{telemetry} object. }
  \item{CTMM}{A \code{ctmm} movement model from the output of \code{ctmm.fit}.}
  \item{method}{Which range distribution method to use. Can be "AKDE" or "AGDE".}
  \item{...}{Arguments passed to the method call or \code{\link{bandwidth}}.}
  \item{R}{A named list of raster covariates if \code{CTMM} contains an RSF model}
  \item{variable}{Not yet supported.}
  \item{error}{Target probability error.}
  \item{res}{Number of grid points along each axis, relative to the location covariance.}
  \item{grid}{Optional grid specification via \code{raster}, \code{UD}, or list of arguments (See `Details' in \code{\link{akde}}).}
}

\details{\code{agde} calculates autocorrelated Gaussian and RSF home-range areas.

Please consult \code{\link{akde}} for further details on \code{method="AKDE"}.}

\value{Returns a \code{UD} object.}

%\references{
%C. H. Fleming and W. F. Fagan and T. Mueller and K. A. Olson and P. Leimgruber and J. M. Calabrese (2015).
%Rigorous home-range estimation with movement data: A new autocorrelated kernel-density estimator.
%\href{http://www.esajournals.org/doi/abs/10.1890/14-2010.1}{Ecology, 96(5), 1182-1188}.

%D. Marcotte. (1996). Fast variogram computation with FFT. \emph{Computers and Geosciences} 22(10), 1175-1186.
%}

\author{ C. H. Fleming. }

%\note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde}}, \code{\link{raster,UD-method}} }

%\examples{\donttest{
%# Load package and data
%library(ctmm)
%data(buffalo)
%cilla <- buffalo[[1]]

%# Fit a continuous-velocity model with tau ~ c(10 days, 1 hour)
%# see help(variogram.fit)
%m2 <- ctmm(tau=c(10*24*60^2,60^2))
%M2 <- ctmm.fit(cilla,m2)

%# Compute akde object
%KD2 <- akde(cilla,M2)

%# Plot data with AKDE contours
%plot(cilla,AKDE=KD2)
%}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
