\name{overlap}
\alias{overlap}
%\alias{overlap.ctmm}
%\alias{overlap.telemetry}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Calculate the overlap between two stationary distributions}

\description{This function calculates a useful measure of similarity between distributions known as the \emph{Bhattacharyya coefficient} in statistics and simply the \emph{fidelity} or \emph{overlap} in quantum and statistical mechanics. It is roughly speaking the ratio of the intersection area to the average individual area. When applied to \code{ctmm} objects, this function returns the overlap of the two Gaussian distributions. When applied to \code{telemetry} or (aligned) \code{UD} objects with corresponding movement models, this function returns the overlap of their (autocorrelated) kernel density estimates.}

\usage{ overlap(object,CTMM=NULL,level=0.95,debias=TRUE,...) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{list} of \code{ctmm} fit, \code{telemetry}, or aligned \code{UD} objects to compare.}
  \item{CTMM}{A \code{list} of \code{ctmm} fit objects corresponding to \code{object} list.}
  \item{level}{The confidence level desired for the output.}
  \item{debias}{Approximate debiasing of the overlap.}
  \item{...}{Additional arguments relevant for \code{\link{akde}}, such as \code{res} and \code{weights}.}
}

%\details{}

\value{A table of confidence intervals on the overlap estimate. A value of \code{1} implies that the two distributions are identical, while a value of \code{0} implies that the two distributions share no area in common. Corresponding \code{ctmm} objects are necessary to provide confidence intervals and debiasing for the point esitmates.}

%% ~put references to the literature/web site here ~
\references{
K. Winner, M. J. Noonan, C. H. Fleming, K. Olson, T. Mueller, D. Sheldon, J. M. Calabrese.
Statistical inference for home range overlap.
\href{https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13027}{Methods in Ecology and Evolution, DOI:10.1111/2041-210X.13027 (2018)}.
}

\author{C. H. Fleming and K. Winner}

\note{Uncertainties in \code{CTMM1} and \code{CTMM2} are propagated into the overlap estimate under the approximation that the Bhattacharyya distance is a chi-square random variable. Debiasing makes further approximations noted in Winner & Noonan et al (2018).}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde}}, \code{\link{ctmm.fit}} }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)

# fit models for first two buffalo
GUESS <- lapply(buffalo[1:2], function(b) ctmm.guess(b,interactive=FALSE) )
FITS <- lapply(1:2, function(i) ctmm.fit(buffalo[[i]],GUESS[[i]]) )
names(FITS) <- names(buffalo[1:2])

# Gaussian overlap between these two buffalo
overlap(FITS)

# AKDE overlap between these two buffalo
overlap(buffalo[1:2],FITS)
}}
