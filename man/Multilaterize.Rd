\name{Multilaterize}
\alias{Multilaterize}

\title{Estimating animal location from UHF and VHF telemetry data}

\description{
Estimates animal locations from received signal strength (RSSI) measurements
collected by multiple UHF or VHF receivers using a gamma-GLM calibration model
and Newton-Raphson algorithm. The function combines simultaneous RSSI
observations from multiple receivers to estimate the most likely animal location
and its associated uncertainty.
}

\usage{
Multilaterize(data, RSSIpar, Start = NULL, precision = 1/2, max_iter = 1e+05, ...)
}

\arguments{
  \item{data}{
  A data frame in which each row represents a single RSSI detection from one receiver.
  Multiple rows with the same timestamp (or localization time) and
  tag ID should correspond to simultaneous detections from different receivers.
 A data frame must have these columns: projected X and Y UTM coordinates of the
 receiver which detected the tag (\code{UTM.x} and \code{UTM.y}), corresponding \code{timestamp},
 recorded signal strength (\code{tag_rssi}), tag ID (\code{tag.ID}).
 Note, a minimum number of receivers should be at least three (3) for efficient multilateration.
 }
 \item{RSSIpar}{
  A data frame containing a range of RSSI measurements and the corresponding Gamma
  distribution parameters. The data frame must contain the following columns:

  \describe{
    \item{\code{RSSI}}{Integer RSSI value.}
    \item{\code{Shape}}{Estimated Gamma shape parameter corresponding to the RSSI value.}
    \item{\code{Lambda}}{Estimated Gamma scale parameter corresponding to the RSSI value.}
  }
}
  \item{Start}{ An optional initial estimate of the animal's \code{(x, y)}
  location used to initialize the Newton-Raphson optimization. If \code{NULL},
  an initial location is computed automatically.
   If supplied it should be a numeric vector of length 2.}
  \item{precision}{Optimization stops when successive location estimates differ by less than convergence precision threshold \code{(default = 1/2)}.}
  \item{max_iter}{Maximum number of Newton-Raphson iterations allowed before optimization terminates \code{(default = 100000)}.}
  \item{\dots}{Additional arguments passed to methods.}
}

\details{
This function is intended to be used after fitting a gamma GLM to calibration data to
establish a relationship between square distances and the integer RSSI values.
The gamma GLM models squared distance rather than distance because this
provides a convenient likelihood with analytic derivatives for optimization.
The fitted model is used to infer Gamma distribution parameters the \code{Shape} and \code{Lambda}
for each row and over the range of integer RSSI values. These parameters are supplied to \code{RSSIpar} object and
provides distribution over probable distances for each RSSI value. The gamma-GLM is in the following form:

\deqn{p_i(D_i^2)= \frac{\lambda_i^{\alpha_i}}{\Gamma(\alpha_i)} D_i^{2(\alpha_i-1)} e^{-\lambda_i \, D_i^2}}

where \eqn{D_i^2 = (\mathbf{r}-\mathbf{r}_i)^2 = {(x-x_i)^2+(y-y_i)^2}},
\eqn{\mathbf{r}=(x,y)} is the unknown location of the animal,
and where \eqn{\alpha_i} and \eqn{\lambda_i} are the shape and
scale parameters of the distance distribution of the RSSI observation \code{i}.
Because estimating the location of the animal at a particular time requires
multiple observations to be recorded simultaneously from
different receivers, we can write the joint log-density of
multiple observations given the gamma-GLM model above as:

\deqn{\log p(\mathbf{r})= \sum_i \left( \alpha'_i \log (\mathbf{r}-\mathbf{r}_i)^2 - \lambda_i (\mathbf{r}-\mathbf{r}_i)^2 \right) + \mathrm{constant}}

where \eqn{\alpha'_i = \alpha_i-\frac{3}{2}} after a change of coordinates.

The location of an animal can then be estimated by numerically maximizing the log-likelihood or
using the fast Newton-Raphson algorithm with the analytic derivatives (gradient vector and Hessian matrix).

}

\value{Returns a \code{telemetry} object: a data frame with Tag ID (\code{tag.local.identifier}),
date and time (\code{timestamp}), estimated locations (\code{UTM.x} and \code{UTM.y}),and
the elements of the estimated covariance matrix of the location estimate obtained
from the inverse Hessian ( \code{COV.xx} and \code{COV.yy},
and \code{COV.xy}).
}

\references{
%% ~put references to the literature/web site here ~
}
\author{M. M. Masolele and C. H. Fleming}

\note{
The smoothing parameters shown in the vignette example produce a monotonic relationship between RSSI and expected squared distance.
Other values may be appropriate for larger calibration datasets provided the fitted relationship remains monotonic.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{\donttest{

# Load package and data
library(ctmm)
data(Owl)

# See vignette for calibration model

# Detection from multiple receivers at the same timestamp
Obsv<-Owl[Owl$Uid==1,]

# Estimate location using multilateration
Results<-Multilaterize(data=Obsv, RSSIpar=Paramdf)




}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
