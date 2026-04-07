\name{index}
\alias{index}
\encoding{UTF-8}
%\alias{overlap.ctmm}
%\alias{overlap.telemetry}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Calculate the expectation distance of an indicator variable}

\description{This function calculates a useful and numerically stable measure of distance between a spatial indicator variable (e.g., HFI or a habitat category) and the utilization distribution of an animal.}

\usage{ index(x,R,...)
%- maybe also 'usage' for other objects documented here.
}
\arguments{
  \item{x}{A \code{UD} or \code{ctmm} fit object.}
  \item{R}{A raster object with values between zero and one.}
  \item{...}{Not currently used.}
}

\details{This function calculates the index \eqn{\mathcal{I} = \sqrt{-\log \mathrm{E} [R(x,y)]}} with respect to the provided utilization distribution,
where \eqn{0 \leq R(x,y) \leq 1} is a spatial indicator covariate.
The expected value of the indicator variable is given by \eqn{\mathrm{E} [R(x,y)] = e^{-\mathcal{I}^2}},
but this calculation suffers from numerical underflow if the covariate, \eqn{R(x,y)}, is outside of the utilization distribution.
When \eqn{R(x,y)} is a patch outside of the support of the utilization distribution, then \eqn{\mathcal{I}} is the relative distance to the covariate.
}

\value{A numeric point estimate.}

%% ~put references to the literature/web site here ~
%\references{}

\author{C. H. Fleming and N. Hirama}

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{akde}} }

%\examples{\donttest{}}
