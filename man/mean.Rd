\name{mean.UD}
\alias{mean.UD}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Average autocorrelated kernel density estimates}
\description{
This function calculates a weighted average of autocorrelated kernel density home-range estimates from \code{akde}. The point estimates are correct, but the confidence-interval calculation is not yet complete.
}
\usage{
\S3method{mean}{UD}(x,weights=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A list of \code{UD}s calculated on the same grid.}
  \item{weights}{A vector of numeric weights with the same length as \code{x}, specifying the relative frequency of each distribution in \code{x}.}
  \item{...}{Additional arguments for future use.}
}

\details{
By default, uniform weights are used (\code{weights=rep(1,length(x))}). This can be sensible for averaging over individuals. For averaging over periods of time, users should consider weighting by the proportion of time spent in each distribution. For example, if an animal spends 4 months in its winter range, \code{x[[1]]}, and 7 months in its summer range, \code{x[[2]]}, then the annual range (sans migration corridor) would be calculated with \code{weights=c(4,7)}.

All UDs need to be calculated on the same grid (see \code{\link{overlap}} for an example).
}

\value{Returns a \code{UD} object: a list with the sampled grid line locations \code{r$x} and \code{r$y}, the extent of each grid cell \code{dr}, the probability density and cumulative distribution functions evaluated on the sampled grid locations \code{PDF} & \code{CDF}, the optimal bandwidth matrix \code{H}, and the effective sample size of the data in \code{DOF.H}.
}

%\references{}

\author{ C. H. Fleming}

\note{ The confidence intervals are not yet calculated correctly.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde}}}

%\examples{}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
