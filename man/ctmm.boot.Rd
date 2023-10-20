\name{ctmm.boot}
\alias{ctmm.boot}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parametric bootstrap continuous-time movement models}
\description{This function allows the point estimates and confidence intervals of an initial estimated movement model to be improved by parametric boostrap, as described in Fleming et al (2019).
}
\usage{
ctmm.boot(data,CTMM,method=CTMM$method,AICc=FALSE,iterate=FALSE,robust=FALSE,error=0.01,
          clamp=0.001,cores=1,trace=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Timeseries data represented as a \code{telemetry} object. }
  \item{CTMM}{ A \code{ctmm} movement-model object from the output of \code{ctmm.fit} containing the initial parameter estimates.}
  \item{method}{Fitting method to use: \code{"ML"}, \code{"HREML"}, \code{"pREML"}, \code{"pHREML"}, or \code{"REML"}. See \code{\link{ctmm.fit}} for descriptions.}
  \item{AICc}{Run dual set of simulations to approximate AICc values via Kullback–Leibler divergence. Otherwise, only the AIC is updated.}
%  \item{multiplicative}{Removes multiplicative bias (rather than additive) bias from positive scale parameters like area.}
  \item{iterate}{Iteratively solve for the parameters such that the average estimate (of \code{method}) is that of the data, whereas with \code{iterate=FALSE} only the first-order correction is calculated from the initial estimate.}
  \item{robust}{Uses robust estimates of the average and covariation for debiasing. Useful when parameters are near boundaries.}
  \item{error}{Relative standard error target for bootstrap ensemble estimates and nonlinear iterations.}
  \item{clamp}{Fix the COV/CoV estimate to the initial COV/CoV estimate, until \code{error} \eqn{\ll} \code{clamp}.}
  \item{cores}{Number of simulations to run in parallel. \code{cores=NULL} will use all cores, while \code{cores<0} will reserve \code{abs(cores)}.}
  \item{trace}{Report progress updates. Can be among \code{0:2} with increasing detail.}
  \item{...}{Further arguments passed to \code{\link{ctmm.fit}}.}
}

%\details{
%\code{ctmm.boot} can leverage multiple estimators via the \code{method} argument (see \code{\link{ctmm.fit}}) and as described in Fleming et al (2019), though generally this is only useful if the specified estimators deviate substantially from each other realtive to the target \code{error}.
%}

\value{A model fit object with relatively unbiased estimates of location covariance, and autocorrelation timescales (and more accurate CIs than \code{ctmm.fit}). If \code{AICc=TRUE}, then, in addition to an updated \code{AICc} slot, the model fit object will also contain a \code{VAR.AICc} slot quantifying the numerical variance in the \code{AICc} estimate. This variance can be decreased by decreasing argument \code{error}.}

\references{
C. H. Fleming, M. J. Noonan, E. P. Medici, J. M. Calabrese,
``Overcoming the challenge of small effective sample sizes in home-range estimation'',
Methods in Ecology and Evolution 10:10, 1679-1689 (2019) \doi{10.1111/2041-210X.13270}.
}

\author{ C. H. Fleming. }

\note{The bootstrapped COV estimates tend to be far more noisy than the bootstrapped point estimates. \code{clamp} can fix the bootstrapped COV/CoV estimate to the initial COV/CoV estimate until the point estimates obtain higher numerical precision.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{ctmm.fit}}. }

\examples{\donttest{
# Load package and data
library(ctmm)
data(gazelle)
DATA <- gazelle[[3]]

GUESS <- ctmm.guess(DATA,interactive=FALSE)
FIT <- ctmm.select(DATA,GUESS)

# some human-readable information
summary(FIT)

# in general, you will want to set iterate=TRUE,trace=TRUE
BOOT <- ctmm.boot(DATA,FIT,iterate=FALSE,trace=FALSE)

# compare to the previous estimate
summary(BOOT)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
