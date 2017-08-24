\name{emulate}
\alias{emulate}
\alias{emulate.ctmm}
\alias{emulate.telemetry}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Draw a random model-fit from the sampling distribution}
\description{
This function generates random model-fit statistics from the sampling distribution of a given \code{ctmm} movement model and sampling schedule.
If \code{fast=FALSE}, the results are exact, though slow to evaluate.
Else if \code{fast=TRUE}, the central-limit theorem is invoked.
}
\usage{
emulate(object,...)

\S3method{emulate}{ctmm}(object,data=NULL,fast=FALSE,...)

\S3method{emulate}{telemetry}(object,CTMM,fast=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{telemetry} data or \code{ctmm} model object.}
  \item{CTMM}{A \code{ctmm} movement-model object.}
  \item{data}{Optional \code{telemetry} object for exact results.}
  \item{fast}{Whether or not to invoke the central-limit theorem.}
  \item{...}{Arguments passed to \code{\link{ctmm.fit}}.}
}
\details{
Given \code{fast=FALSE}, which requires the \code{data} argument specified, new data are simulated from the \code{CTMM} movement model with the same sampling schedule and error structure as \code{data}. A new model, of the same structure as \code{CTMM}, is then fit to the simulated data and returned.

Given \code{fast=TRUE}, a model-fit object is sampled from the central-limit distribution, using the covariance estimates within \code{CTMM}.
Strictly positive parametes, such as area, are log-transformed prior to the normal approximation.
Note that this faster method does not adjust for bias.
}

\value{ A \code{ctmm} movement model with the same structure as \code{CTMM}.}

%\references{}


\author{ C. H. Fleming. }

%\note{ }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{ctmm.fit}}, \code{\link{simulate.ctmm}} }

%\examples{}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
