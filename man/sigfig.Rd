\name{format}
\alias{sigfig}
\alias{dimfig}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Scientific formatting of numbers}
\description{
Functions for concisely representing dimensionful quantities and uncertain quantities.
}
\usage{
dimfig(data,dimension,thresh=1,...)

sigfig(est,VAR=NULL,SD=NULL,level=0.95,digits=2,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A numerical vector of dimensionful quantities represented in SI units.}
  \item{dimension}{One of \code{"length"}, \code{"area"}, \code{"time"}, \code{"frequency"}, \code{"speed"}, \code{"diffusion"}, or \code{"mass"}.}
  \item{thresh}{Threshold quantity for switching between units. E.g., 100 cm is represented as 1 m only if \code{thresh>=1}.}
  \item{est}{Can be either confidence-interval estimates with rows {(lower-limit,point-estimate,upper-limit)} or point estimates (with \code{VAR} or \code{SD} also specified).}
  \item{VAR}{Variance in the sampling distribution of \code{x}.}
  \item{SD}{Standard deviation in the sampling distribution of \code{x}.}
  \item{level}{Confidence level for designating the numerical precision of the significant digits.}
  \item{digits}{Number of significant digits to retain.}
  \item{...}{Not currently used.}
}

\details{\code{dimfig} chooses the set of units that provides the most concise representation for \code{data}, and \code{sigfig} concisely represents statistical estimates with a fixed number of significant digits.}


\value{\code{dimfig} returns a \code{list} with slots for the converted \code{data} and the name of the most concise units. \code{sigfig} returns a character string that is formated with the specified number of significant digits.}

%\references{
%C. H. Fleming and W. F. Fagan and T. Mueller and K. A. Olson and P. Leimgruber and J. M. Calabrese (2015).
%Rigorous home-range estimation with movement data: A new autocorrelated kernel-density estimator.
%\href{http://www.esajournals.org/doi/abs/10.1890/14-2010.1}{Ecology, 96(5), 1182-1188}.

%D. Marcotte. (1996). Fast variogram computation with FFT. \emph{Computers and Geosciences} 22(10), 1175-1186.
%}

\author{ C. H. Fleming. }

% \note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{\%#\%}} }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)
DATA <- buffalo$Cilla

GUESS <- ctmm.guess(DATA,interactive=FALSE)
# in general, you want to run ctmm.select instead
FIT <- ctmm.fit(DATA,GUESS)

# raw summary (SI units)
summary(FIT,units=FALSE)

# default summary (concise units)
summary(FIT,units=TRUE)

# text-formatted summary
sigfig( summary(FIT)$CI )
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
