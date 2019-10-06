\name{mean.variogram}
\alias{mean.variogram}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute a number-weighted average of variogram objects}
\description{This function takes a list of variogram objects and calculates its number-weighted average variogram.}

\usage{
\S3method{mean}{variogram}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A \code{variogram} object or list of such objects to be averaged.}
  \item{...}{Additional variograms if specified individually.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Returns a \code{variogram} object which is a dataframe containing the lag, the semi-variance estimate at that lag, and the approximate degrees of freedom associated with the semi-variance estimate.
}
\references{
C. H. Fleming, J. M. Calabrese, T. Mueller, K.A. Olson, P. Leimgruber, W. F. Fagan,
``From fine-scale foraging to home ranges: A semi-variance approach to identifying movement modes across spatiotemporal scales'',
\href{http://www.jstor.org/discover/10.1086/675504}{The American Naturalist, 183:5, E154-E167 (2014)}.
}
\author{
J. M. Calabrese and C. H. Fleming
}
\note{Variogram averaging should only be used when there is a degree of similarity across individual variograms.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{plot.variogram}}, \code{\link{variogram}}. }

\examples{\donttest{
# Load package and data
library(ctmm)
data(buffalo)

# Calculate a list of variograms for all similar individuals in the dataset
# the 4th buffalo has a different sampling rate
SVFS <- lapply( buffalo[-4] , variogram )
# alternatively, we could variogram all at coarsest scale with variogram option dt

# Calculate the average variogram
SVF <- mean(SVFS)

# Plot the mean variogram
plot(SVF)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
