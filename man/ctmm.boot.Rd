\name{ctmm.boot}
\alias{ctmm.boot}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parametric bootstrap continuous-time movement models}
\description{This function allows the point estimates and confidence intervals of an initial estimated movement model to be improved by parametric boostrap,  as described in Fleming et al (2018).
}
\usage{
ctmm.boot(data,CTMM,method=c('HREML','pREML','pHREML'),error=0.01,mc.cores=1,trace=TRUE,
          ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ Timeseries data represented as a \code{telemetry} object. }
  \item{CTMM}{ A \code{ctmm} movement-model object from the output of \code{ctmm.fit} containing the initial parameter estimates.}
  \item{method}{Fitting method(s) to use: \code{"ML"}, \code{"HREML"}, \code{"pREML"}, \code{"pHREML"}, or \code{"REML"}. See "Description" below.}
  \item{error}{Relative standard error target for bootstrap ensemble estimates and nonlinear iterations.}
  \item{mc.cores}{Number of processes to run in parallel (UNIX only). \code{mc.cores=NULL} will use all present cores.}
  \item{trace}{Report progress updates. Can be among \code{0:2} with increasing detail.}
  \item{...}{Further arguments passed to \code{ctmm.fit}.}
}

\details{
\code{ctmm.boot} can leverage multiple estimators via the \code{method} argument (see \code{\link{ctmm.fit}}) and as described in Fleming et al (2018), though generally this is only useful if the specified estimators deviate substantially from each other realtive to the prescribed \code{error}.
}

\value{A model fit object with relatively unbiased estimates of home-range area, diffusion rate, and RMS speed (and more accurate CIs than \code{ctmm.fit}).}

%\references{}

\author{ C. H. Fleming. }

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{ctmm.fit}}. }

%\examples{\donttest{
%# Load package and data
%library(ctmm)
%data(buffalo)
%Cilla <- buffalo$Cilla
%
%GUESS <- ctmm.guess(Cilla,interactive=FALSE)
%FIT <- ctmm.fit(Cilla,GUESS)
%
%# some human-readable information
%summary(FIT)
%}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
