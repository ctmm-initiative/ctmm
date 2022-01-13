\name{meta}
\alias{meta}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Meta-analysis of movement-model parameters}
\description{These functions estimate population-level mean parameters from individual movement models and related estimates, including AKDE home-range areas, while taking into account estimation uncertainty.
}
\usage{
meta(x,variable="area",level=0.95,level.UD=0.95,method="MLE",IC="AICc",boot=FALSE,
     error=0.01,debias=TRUE,verbose=FALSE,units=TRUE,plot=TRUE,sort=FALSE,mean=TRUE,
     col="black",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A list of \code{ctmm} movement-model objects, \code{UD} objects, or \code{UD} \code{summary} output, constituting a sampled population, or a list of such lists, each constituting a sampled sub-population.}
  \item{variable}{Variable of interest. Can be \code{"area"}, \code{"diffusion"}, \code{"speed"}, \code{"tau position"}, or \code{"tau velocity"}.}
  \item{level}{Confidence level for parameter estimates.}
  \item{level.UD}{Coverage level for home-range estimates. E.g., 50\% core home range.}
%  \item{level.pop}{Coverage level for the population of individuals.}
  \item{method}{Statistical estimator used---either maximum likelihood estimation based (\code{"MLE"}) or approximate `best linear unbiased estimator' (\code{"BLUE"})---for comparison purposes.}
  \item{IC}{Information criterion to determine whether or not population variation can be estimated. Can be \code{"AICc"}, \code{AIC}, or \code{"BIC"}.}
  \item{boot}{Perform a parametric bootstrap for confidence intervals and first-order bias correction if \code{debias=TRUE}.}
  \item{error}{Relative error tolerance for parametric bootstrap.}
  \item{debias}{Apply Bessel's inverse-Gaussian correction and various other bias corrections if \code{method="MLE"}, REML if \code{method="BLUE"}, and an additional first-order correction if \code{boot=TRUE}.}
  \item{verbose}{Return a list of both cross-population and sub-population analyses if \code{TRUE} and \code{x} is a list of sub-population lists.}
  \item{units}{Convert result to natural units.}
  \item{plot}{Generate a meta-analysis forest plot.}
  \item{sort}{Sort individuals by their point estimates in forest plot.}
  \item{mean}{Include population mean estimate in forest plot.}
  \item{col}{Color(s) for individual labels and error bars.}
  \item{...}{Further arguments passed to \code{plot}.}
}

\details{
\code{meta} employs a custom \eqn{\chi^2-IG} hierarchical model to calculate relatively unbiased population mean estimates of positive scale parameters,
including home-range areas, diffusion rates, mean speeds, and autocorrelation timescales.
Population coefficient of variation (CoV) estimates are also provided.
More details will be provided in an upcomming manuscript.
}

\value{If \code{x} constitutes a sampled population, then \code{meta} returns a table with rows corresponding to the population mean and coefficient of variation.
If \code{x} constitutes a list of sampled sub-populations, then \code{meta} returns confidence intervals on the sub-population mean ratios.}

%\section{Warnings}{}

%\references{}


\author{ C. H. Fleming. }

\note{The AICc formula is approximated via the Gaussian relation.

Confidence intervals depicted in the forest plot are \eqn{\chi^2} and may differ from the output of \code{summary()} in the case of mean speed and timescale parameters.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde}}, \code{\link{cluster}}, \code{\link{ctmm.fit}}. }

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

# meta-analysis of buffalo home-range areas
meta(AKDES,col=c(COL,'black'),sort=TRUE)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
