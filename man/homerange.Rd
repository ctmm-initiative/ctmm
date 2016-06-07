\name{homerange}
\alias{homerange}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate a range distribution estimate
}
\description{
This function estimates the range distribution from \code{telemetry} data and a continuous-time movement model.
}
\usage{
homerange(data,CTMM,method="AKDE",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ 2D timeseries telemetry data represented as a \code{telemetry} object. }
  \item{CTMM}{A \code{ctmm} movement model from the output of \code{ctmm.fit}.}
  \item{method}{Which range distribution method to use. For now only "AKDE" is supported.}
  \item{...}{Arguments passed to the method call.}
}

\details{Please consult \code{\link{akde}} for further details.}

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
