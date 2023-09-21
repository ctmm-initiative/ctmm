\name{overlap}
\alias{overlap}
\encoding{UTF-8}
%\alias{overlap.ctmm}
%\alias{overlap.telemetry}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Calculate the overlap between two stationary distributions}

\description{This function calculates a useful measure of similarity between distributions known as the \emph{Bhattacharyya coefficient} in statistics and simply the \emph{fidelity} or \emph{overlap} in quantum and statistical mechanics. It is roughly speaking the ratio of the intersection area to the average individual area, but it is a direct comparison between the density functions and does not require an arbitrary quantile to be specified. When applied to \code{ctmm} objects, this function returns the overlap of the two Gaussian distributions. When applied to aligned \code{UD} objects with corresponding movement models, this function returns the overlap of their (autocorrelated) kernel density estimates.}

\usage{ overlap(object,method="Bhattacharyya",level=0.95,debias=TRUE,...) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{list} of \code{ctmm} fit or aligned \code{UD} objects to compare.}
  \item{method}{Can be \code{"Bhattacharyya"} or \code{"Encounter"} (see Details below).}
%  \item{CTMM}{A \code{list} of \code{ctmm} fit objects corresponding to \code{object} \code{telemetry} list.}
  \item{level}{The confidence level desired for the output.}
  \item{debias}{Approximate debiasing of the overlap.}
  \item{...}{Not currently used.}
}

\details{The default \code{method="Bhattacharyya"} estimates the standard overlap measure \eqn{\int\int \sqrt{p(x,y) \, q(x,y)} \, dx \, dy} between the distributions \eqn{p(x,y)} and \eqn{q(x,y)},
while \code{method="encounter"} estimates the non-standard measure \eqn{\frac{\int\int p(x,y) \, q(x,y) \, dx \, dy}{\sqrt{\int\int p(x',y')^2 \, dx' dy' \int\int q(x'',y'')^2 \, dx'' dy''}}},
which has a numerator proportional to the uncorrelated encounter probability.
Both measures lie between 0 and 1, where 0 indicates no shared support and 1 indicates identical distributions.
}

\value{An object with slots \code{DOF}, containing the effective sample sizes, and \code{CI} containing a table of confidence intervals on the overlap estimates. A value of \code{1} implies that the two distributions are identical, while a value of \code{0} implies that the two distributions share no area in common.}

%% ~put references to the literature/web site here ~
\references{
K. Winner, M. J. Noonan, C. H. Fleming, K. Olson, T. Mueller, D. Sheldon, J. M. Calabrese.
``Statistical inference for home range overlap'',
Methods in Ecology and Evolution, 9:7, 1679-1691 (2018) \doi{10.1111/2041-210X.13027}.
}

\author{C. H. Fleming and K. Winner}

\note{In \code{ctmm} v0.5.2, direct support for \code{telemetry} objects was dropped and the \code{CTMM} argument was depreciated for \code{UD} objects, simplifying usage.

Uncertainties in the model fits are propagated into the overlap estimate under the approximation that the Bhattacharyya distance is a chi-square random variable. Debiasing makes further approximations noted in Winner & Noonan et al (2018).}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde}}, \code{\link{ctmm.fit}}, \code{\link{distance}}, \code{\link{encounter}} }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)

# fit models for first two buffalo
GUESS <- lapply(buffalo[1:2], function(b) ctmm.guess(b,interactive=FALSE) )
# using ctmm.fit here for speed, but you should almost always use ctmm.select
FITS <- lapply(1:2, function(i) ctmm.fit(buffalo[[i]],GUESS[[i]]) )
names(FITS) <- names(buffalo[1:2])

# Gaussian overlap between these two buffalo
overlap(FITS)

# AKDE overlap between these two buffalo
# create aligned UDs
UDS <- akde(buffalo[1:2],FITS)
# evaluate overlap
overlap(UDS)
}}
