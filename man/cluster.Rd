\name{cluster}
\alias{cluster}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Clustering of movement-model parameters}
\description{These functions cluster and classify individual movement models and related estimates, including AKDE home-range areas, while taking into account estimation uncertainty.
}
\usage{
cluster(x,level=0.95,level.UD=0.95,debias=TRUE,IC="BIC",units=TRUE,plot=TRUE,sort=FALSE,
        ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A list of \code{ctmm} movement-model objects, \code{UD} objects, or \code{UD} \code{summary} output, constituting a sampled population, or a list of such lists, each constituting a sampled sub-population.}
  \item{level}{Confidence level for parameter estimates.}
  \item{level.UD}{Coverage level for home-range estimates. E.g., 50\% core home range.}
  \item{debias}{Apply Bessel's inverse-Gaussian correction and various other bias corrections.}
  \item{IC}{Information criterion to determine whether or not population variation can be estimated. Can be \code{"AICc"}, \code{AIC}, or \code{"BIC"}.}
  \item{units}{Convert result to natural units.}
  \item{plot}{Generate a meta-analysis forest plot with two means.}
  \item{sort}{Sort individuals by their point estimates in forest plot.}
  \item{...}{Further arguments passed to \code{plot}.}
}

\details{
So-far only the clustering of home-range areas is implemented. More details will be provided in an upcomming manuscript.
}

\value{A \code{list} with elements \code{P} and \code{CI},
where \code{P} is an array of individual membership probabilities for sub-population 1,
and \code{CI} is a table with rows corresponding to the sub-population means, coefficients of variation, and membership probabilities, and the ratio of sub-population means.
}

%\section{Warnings}{}

%\references{}


\author{ C. H. Fleming. }

\note{The AICc formula is approximated via the Gaussian relation.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde}}, \code{\link{ctmm.fit}}, \code{\link{meta}}. }

\examples{\donttest{
# load package and data
library(ctmm)
data(buffalo)

# fit movement models
FITS <- AKDES <- list()
for(i in 1:length(buffalo))
{
  GUESS <- ctmm.guess(buffalo[[i]],interactive=FALSE)
  # use ctmm.select unless you are certain that the selected model is OUF
  FITS[[i]] <- ctmm.fit(buffalo[[i]],GUESS,trace=2)
}

# calculate AKDES on a consistent grid
AKDES <- akde(buffalo,FITS,trace=1)

# color to be spatially distinct
COL <- color(AKDES,by='individual')

# plot AKDEs
plot(AKDES,col.DF=COL,col.level=COL,col.grid=NA,level=NA)

# cluster-analysis of buffalo
cluster(AKDES,sort=TRUE)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
