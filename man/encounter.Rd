\name{encounter}
\alias{encounter}
\alias{cde}
\encoding{UTF-8}
%\alias{overlap.ctmm}
%\alias{overlap.telemetry}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Encounter statistics}

\description{Functions to calculate encounter probabilities and the conditional location distribution of where encounters take place (conditional on said encounters taking place), as described in Noonan et al (2021).}

\usage{encounter(object,debias=FALSE,level=0.95,normalize=FALSE,self=TRUE,...)

cde(object,include=NULL,exclude=NULL,debias=FALSE,...) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{list} of aligned \code{UD} objects.}
  \item{debias}{Approximate bias corrections [IN DEVELOPMENT].}
  \item{level}{Confidence level for relative encounter rates.}
  \item{normalize}{Normalize relative encounter rates by the average uncorrelated self-encounter rate.}
  \item{self}{Fix the self-interaction rate appropriately.}
  \item{include}{A matrix of interactions to include in the calculation (see Details below).}
  \item{exclude}{A matrix of interactions to exclude in the calculation (see Details below).}
  \item{...}{Additional arguments for future use.}
}

\details{
Encounter probabilities are standardized to 1 meter, and must be multiplied by the square encounter radius (in meters), to obtain other values. If \code{normalize=FALSE}, the relative encounter rates have units of \eqn{1/m^2} and tend to be very small numbers for very large home-range areas. If \code{normalize=TRUE}, the relative encounter rates are normalized by the average uncorrelated self-encounter rate, which is an arbitrary value that provides a convenient scaling.

The \code{include} argument is a matrix that indicates which interactions are considered in the calculation.
By default, \code{include = 1 - diag(length(object))}, which implies that all interactions are considered aside from self-interactions. Alternatively, \code{exclude = 1 - include} can be specified, and is by-default \code{exclude = diag(length(object))}, which implies that only self-encounters are excluded.
}

\value{\code{encounter} produces an array of standardized encounter probabilities with CIs, while \code{cde} produces a single \code{UD} object.}

%% ~put references to the literature/web site here ~
\references{
M. J. Noonan, R. Martinez-Garcia, G. H. Davis, M. C. Crofoot, R. Kays, B. T. Hirsch, D. Caillaud, E. Payne, A. Sih, D. L. Sinn, O. Spiegel, W. F. Fagan, C. H. Fleming, J. M. Calabrese, ``Estimating encounter location distributions from animal tracking data'', Methods in Ecology and Evolution (2021) \doi{10.1111/2041-210X.13597}.
}

\author{C. H. Fleming}

\note{Prior to v1.2.0, \code{encounter()} calculated the CDE and \code{rates()} calculated relative encounter probabilities.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde}}, \code{\link{overlap}} }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)

# fit models for first two buffalo
GUESS <- lapply(buffalo[1:2], function(b) ctmm.guess(b,interactive=FALSE) )
# in general, you should use ctmm.select here
FITS <- lapply(1:2, function(i) ctmm.fit(buffalo[[i]],GUESS[[i]]) )
names(FITS) <- names(buffalo[1:2])

# create aligned UDs
UDS <- akde(buffalo[1:2],FITS)

# calculate 100-meter encounter probabilities
P <- encounter(UDS)
P$CI * 100^2

# calculate CDE
CDE <- cde(UDS)

# plot data and encounter distribution
plot(buffalo[1:2],col=c('red','blue'),UD=CDE,col.DF='purple',col.level='purple',col.grid=NA)
}}
