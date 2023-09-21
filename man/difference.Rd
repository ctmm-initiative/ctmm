\name{difference}
\alias{difference}
\alias{midpoint}
\alias{proximity}
\alias{distances}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate the proximity of two individuals}
\description{
Given a pair of \code{telemetry} objects and \code{ctmm} movement models, predict their location differences or midpoints at shared times and estimate their distances.
}
\usage{
difference(data,CTMM,t=NULL,...)

midpoint(data,CTMM,t=NULL,complete=FALSE,...)

distances(data,CTMM,t=NULL,level=0.95,...)

proximity(data,CTMM,t=NULL,GUESS=ctmm(error=TRUE),debias=TRUE,level=0.95,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{list} of two \code{telemetry} objects.}
  \item{CTMM}{A \code{list} of two \code{ctmm} movement-model objects.}
  \item{t}{An optional vector of times or range of times over which to predict the location differences.}
  \item{complete}{Additionally calculate timestamps and geographic coordinates.}
  \item{level}{Confidence level for the distance/proximity estimate.}
  \item{GUESS}{An optional \code{ctmm} object to specify the candidate model parameters of the location differences.}
  \item{debias}{Include inverse-\eqn{\chi^2} bias corrections.}
  \item{...}{Options passed to \code{\link{ctmm.select}}.}
}
\details{The \code{difference} function predicts the location difference vectors, \eqn{(x_A-x_B,y_A-y_B)}, for a pair of individuals, \eqn{\{A,B\}}, at overlapping times. The \code{midpoint} function predicts the location midpoints, \eqn{(x_A+x_B,y_A+y_B)/2}, for a pair of individuals. The \code{distances} function further estimates the instantaneous distances between individuals. The \code{proximity} function fits an autocorrelation model to the output of \code{difference}, and then compares the mean-square distance between the individuals to what you would expect if the two individuals were moving independently.
}

\value{\code{difference} and \code{midpoint} output \code{telemetry} objects of the location differences and midpoints with prediction covariances. \code{distances} outputs a \code{data.frame} of distance estimates with confidence intervals. \code{proximity} outputs a ratio estimate with confidence intervals, where values \code{<1} indiciate that the two individuals are closer on average than expected for independent movement, \code{1} is consistent with independent movement, and values \code{>1} indicate that the individuals are farther from each other on average than expected for independent movement. Therefore, if the CIs contain 1, then the distance is insignificant with a p-value threshold of \code{1-level} (two-sided) or half that for a one-sided test.}

%\references{}

\author{ C. H. Fleming. }

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{ctmm.select}}, \code{\link{predict.ctmm}} }

\examples{\donttest{
#Load package
library(ctmm)

# load buffalo data
data(buffalo)

# select two buffalo that overlap in space and time
DATA <- buffalo[c(1,3)]
# plot the two buffalo
plot(DATA,col=c('red','blue'))

FITS <- list()
for(i in 1:2)
{
  GUESS <- ctmm.guess(DATA[[i]],interactive=FALSE)
  # in general, you want to use ctmm.select
  FITS[[i]] <- ctmm.fit(DATA[[i]],GUESS)
}

# calculate difference vectors
DIFF <- difference(DATA,FITS)
# plot the difference vectors with prediction-error ellipses
plot(DIFF)

# calculate the proximity statistic
# disabling location error for speed
proximity(DATA,FITS,GUESS=ctmm(error=FALSE))
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
