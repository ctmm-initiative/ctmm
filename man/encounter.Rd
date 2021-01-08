\name{encounter}
\alias{encounter}
\encoding{UTF-8}
%\alias{overlap.ctmm}
%\alias{overlap.telemetry}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Calculate the conditional location distribution of ecounters}

\description{Functions to calculate the location distribution of where encounters take place, conditional on said encounters taking place, as described in Noonan et al (2020).}

\usage{ encounter(object,include=NULL,exclude=NULL,debias=FALSE,...) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{list} of aligned \code{UD} objects.}
  \item{include}{A matrix of interactions to include in the calculation (see Details below).}
  \item{exclude}{A matrix of interactions to exclude in the calculation (see Details below).}
  \item{debias}{Approximate GRF bias correction (in development).}
  \item{...}{Additional arguments for future use.}
}

\details{
The \code{include} argument is a matrix that indicates which interactions are considered in the calculation.
By default, \code{include = 1 - diag(length(object))}, which implies that all interactions are considered aside from self-interactions. Alternatively, \code{exclude = 1 - include} can be specified, and is by-default \code{exclude = diag(length(object))}, which implies that only self-encounters are excluded.
}

\value{A \code{UD} object.}

%% ~put references to the literature/web site here ~
\references{
M. J. Noonan, R. Martinez-Garcia, G. H. Davis, M. C. Crofoot, R. Kays, B. T. Hirsch, D. Caillaud, E. Payne, A. Sih, D. L. Sinn, O. Spiegel, W. F. Fagan, C. H. Fleming, J. M. Calabrese, ``Estimating encounter location distributions from animal tracking data'', bioRxiv 2020.08.24.261628 (2020) \doi{10.1101/2020.08.24.261628}.
}

\author{C. H. Fleming}

%\note{}

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

# calculate CDE
CDE <- encounter(UDS)

# plot data and encounter distribution
plot(buffalo[1:2],col=c('red','blue'),UD=CDE,col.DF='purple',col.level='purple',col.grid=NA)
}}
