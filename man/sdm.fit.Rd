\name{sdm.fit}
\alias{sdm.fit}
\alias{sdm.select}
\alias{sdm.integrate}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit species distribution models (SDMs) [IN DEVELOPMENT]}
\description{This function fits species distribution models, sampling density models, and integrated SDMs.
}
\usage{
sdm.fit(data,R=list(),formula=NULL,area=NULL,reference="auto",standardize=TRUE,
        integrator="MonteCarlo",error=0.01,max.mem="1 Gb",interpolate=TRUE,trace=TRUE,...)

sdm.select(data,R=list(),formula=NULL,area=NULL,verbose=FALSE,IC="AICc",trace=TRUE,...)

sdm.integrate(biased=NULL,bias=NULL,unbiased=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A \code{telemetry} object.}
  \item{R}{A named list of rasters or time-varying raster stacks [NOT TESTED] to fit Poisson regression coefficients to (under a log link).}
  \item{formula}{Formula object for \eqn{\log(\lambda)} referencing the elements of \code{R} and columns of \code{data} (see Details below). If not specified, a linear term will be included for every element of \code{R}.}
    \item{area}{A spatial polygon object defining the extent of the SDM. If left \code{NULL}, an integrated Gaussian model will be used to define the extent of the SDM, which can be a very bad model for geographic ranges.}
  \item{reference}{When expanding categorical predictors into indicator variables, \code{reference="auto"} will choose the most common predictor to be the reference category. Otherwise, the reference category can be specified by this argument.}
  \item{standardize}{For numerical stability, predictors are \emph{internally} standardized, if \code{rescale=TRUE} and no \code{formula} is specified. (The final outputs are not standardized.) \bold{Otherwise, users are responsible for standardizing their predictors.}}
  \item{integrator}{Numerical integrator used for likelihood evaluation. Can be \code{"MonteCarlo"} or \code{"Riemann"} (IN TESTING).}
  \item{error}{Relative numerical error threshold for the parameter estimates and log-likelihood.}
  \item{max.mem}{Maximum amount of memory to allocate for availability sampling.}
  \item{interpolate}{Whether or not to interpolate raster values during extraction.}
  \item{trace}{Report progress on convergence (see Details).}
  \item{verbose}{Returns all candidate models if \code{TRUE}. Otherwise, only the \code{IC}-best model is returned.}
  \item{IC}{Model selection criterion. Can be AIC, AICc, or BIC.}
  \item{...}{Arguments passed to \code{rsf.fit} or \code{\link{optimizer}}.}
  \item{biased}{A biased SDM calculated from occurrence records with non-uniform sampling density.}
  \item{bias}{An ``SDM'' calculated from data representative of the above sampling density.}
  \item{unbiased}{An unbiased SDM or list of RSFs.}
}

\details{
If \code{data} does not contain a \code{count} column, then a presence-only inhomogeneous Poisson point process (iPPP) model is fit to the data.
Whereas, if \code{data} contains a \code{count} column, then a gamma-mixture Poisson (negative binomial) model is fit to the data.
In either case, the \code{formula} object corresponds to \eqn{\log(\lambda)}, where \eqn{\lambda(x,y)} is the intensity function and expected density.

Instead of specifying a number of ``available'' points to sample and having an unknown amount of numerical error to contend with, \code{sdm.fit} specifies an estimation target \code{error} and the number of ``available'' points is increased until this target is met. Moreover, the output log-likelihood is that of the continuous Poisson point process, which does not depend on the number of ``available'' points that were sampled, though the numerical variance estimate is recorded in the \code{VAR.loglike} slot of the fit object.

When \code{trace=TRUE}, a number of convergence estimates are reported, including the standard deviation of the numerical error of the log-likelihood, SD[\eqn{\log(\ell)}], the most recent log-likelihood update, d\eqn{\log(\ell)}, and the most recent (relative) parameter estimate updates d\eqn{\hat{\beta}/}SD[\eqn{\hat{\beta}}].

The \code{formula} object determines \eqn{\log(\lambda)} and can reference static rasters in \code{R}, time-dependent raster stacks in \code{R} [NOT TESTED], and time-dependent effect modifiers in the columns of \code{data}, such as provided by \code{\link{annotate}}.
Any \code{offset} terms are applied under a log transformation (or multiplicatively to \eqn{\lambda}), and can be used to enforce hard boundaries,
where \code{offset(raster)=TRUE} denotes accesible points and \code{offset(raster)=FALSE} denotes inaccessible points [NOT TESTED].
Intercept terms are ignored, as they generally do not make sense for individual Poisson point process models.
This includes terms only involving the columns of \code{data}, as they lack spatial dependence.

Categorical raster variables are expanded into indicator variables, according to the \code{reference} category argument.
Upon import via \code{\link[raster]{raster}}, categorical variables may need to be assigned with \code{\link[raster]{as.factor}}, or else they may be interpreted as numerical variables.
}

%\value{}

%\section{Warnings}{}

\references{
J. M. Alston, C. H. Fleming, R. Kays, J. P. Streicher, C. T. Downs, T. Ramesh, B. Reineking, & J. M. Calabrese, ``Mitigating pseudoreplication and bias in resource selection functions with autocorrelation-informed weighting'', Methods in Ecology and Evolution 14:2  643--654 (2023) \doi{10.1111/2041-210X.14025}.
}

\author{ C. H. Fleming }

\note{It is much faster to calculate all predictors ahead of time and specifying them in the \code{R} list than to reference then in the \code{formula} argument, which will calculate them as needed, saving memory.

AIC and BIC values for \code{integrated=FALSE} models do not include any penalty for the estimated location and shape of the available area, and so their AIC and BIC values  are expected to be \emph{worse} than reported.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{rsf.fit}}, \code{\link{optimizer}}, \code{\link{summary.ctmm}}. }

%\examples{}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
