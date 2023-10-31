\name{mean.ctmm}
\alias{mean.ctmm}
\alias{mean.UD}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Average movement models and autocorrelated kernel density estimates}
\description{
These functions calculate population averages of continuous-time movement models and utilization distributions.}
\usage{
\S3method{mean}{ctmm}(x,formula,weights=NULL,sample=TRUE,debias=TRUE,IC="AIC",trace=TRUE,...)

\S3method{mean}{UD}(x,weights=NULL,sample=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A list of \code{ctmm} objects calculated in the same projection or \code{UD} objects calculated on the compatible grids.}
  \item{formula}{A \code{formula} object with the predictors of a functional response for RSF models [IN DEVELOPMENT].}
  \item{weights}{A vector of numeric weights with the same length as \code{x}, specifying the relative frequency of each distribution in \code{x}.}
  \item{sample}{\code{x} represents a sample of a larger population if \code{TRUE}, or the entire statistical population if \code{FALSE}.}
  \item{debias}{Include \eqn{\log-\chi^2} and REML bias corrections.}
  \item{IC}{Model selection criterion for the anisotropy of the distribution of mean locations and covariance matrices.}
  \item{trace}{Report location and autocovariance model selection results.}
  \item{...}{Additional arguments for future use.}
}

\details{
When applied to a list of \code{ctmm} objects, \code{mean} calculates an average movement model with populaton variability estimates.
The population model is taken to be multivariate normal and log-normal.
The population mean location represents an arithmetic mean, while the population mean home-range areas, RMS speeds, and diffusion rates represent geometric means.
Location-error estimates are not correctly averaged yet.

When applied to a list of \code{UD} objects, \code{mean} calculates a weighted average of autocorrelated kernel density home-range estimates from \code{akde}. The point estimates are correct, but the confidence-interval calculation is not yet complete.

By default, uniform weights are used (\code{weights=rep(1,length(x))}). This can be sensible for averaging over individuals. For averaging over periods of time, users should consider weighting by the proportion of time spent in each distribution. For example, if an animal spends 4 months in its winter range, \code{x[[1]]}, and 7 months in its summer range, \code{x[[2]]}, then the annual range (sans migration corridor) would be calculated with \code{weights=c(4,7)}.

All UDs need to be calculated on the same grid (see \code{\link{overlap}} for an example).
}

\value{
When applied to a list of \code{ctmm} objects, \code{mean} returns a \code{ctmm} object with additional population variability parameter estimates.

When applied to a list of \code{UD} objects, \code{mean} returns a \code{UD} object: a list with the sampled grid line locations \code{r$x} and \code{r$y}, the extent of each grid cell \code{dr}, the probability density and cumulative distribution functions evaluated on the sampled grid locations \code{PDF} & \code{CDF}, the optimal bandwidth matrix \code{H}, and the effective sample size of the data in \code{DOF.H}.
}

%\references{}

\author{ C. H. Fleming}

%\note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde}}, \code{\link{ctmm.select}} }

%\examples{}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
