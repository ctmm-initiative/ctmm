\name{plot.variogram}
\alias{plot.variogram}
\alias{zoom,variogram-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting methods for variogram objects.}
\description{
Produces simple plots of \code{varigram} objects (semi-variance vs. time lag) and model semi-variance functions, with approximate confidence intervals around the semi-variance estimates.
}
\usage{
\S3method{plot}{variogram}(x,CTMM=NULL,level=0.95,fraction=0.5,col="black",col.CTMM="red",...)

\S4method{zoom}{variogram}(x,fraction=0.5,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A \code{variogram} object calculated using \code{\link{variogram}}.}
  \item{CTMM}{A \code{ctmm} movement model object in the same format as the output of \code{ctmm.fit} or \code{variogram.fit}.}
  \item{level}{Confidence level of confidence bands (95\% default CIs). Can be an array.}
  \item{fraction}{
The proportion of the variogram object, \code{variogram}, that will be plotted. By convention, half is shown. The tail end is generally garbage.
}
  \item{col}{Color for the empirical variogram. Can be an array.}
  \item{col.CTMM}{Color for the model. Can be an array.}
  \item{\dots}{
Additional \code{plot} function parameters.
}
}

%\details{ }

\value{Returns a plot of semi-variance vs. time lag, with the empirical variogram in blue and the \code{ctmm} semi-variance function in red if specified. \code{zoom} includes a log-scale zoom slider to manipulate \code{fraction}.}

\references{
C. H. Fleming, J. M. Calabrese, T. Mueller, K.A. Olson, P. Leimgruber, and W. F. Fagan. (2014). From fine-scale foraging to home ranges: A semi-variance approach to identifying movement modes across spatiotemporal scales. \href{http://www.jstor.org/discover/10.1086/675504}{The American Naturalist, 183(5), E154-E167.}
}
\author{J. M. Calabrese and C. H. Fleming}

\note{The errors of the empirical variogram are correlated. Smooth trends are not necessarily significant.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{ctmm.fit}}, \code{\link{plot}}, \code{\link{variogram}}, \code{\link{variogram.fit}}. }

\examples{
# Load package and data
library(ctmm)
data(buffalo)

# Extract movement data for a single animal
cilla <- buffalo[[1]]

# Calculate variogram
SVF <- variogram(cilla)

# Plot the variogram
plot(SVF)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
