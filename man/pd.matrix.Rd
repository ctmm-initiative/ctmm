\name{pd.solve}
\alias{pd.solve}
\alias{pd.logdet}
\alias{pd.sqrtm}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Postive-definite matrix operations}
\description{These functions provide matrix operations for positive-definite and positive-semidefinite matrices (e.g., covariance matrices) that can better handle ill-conditioned matrices, which is a common problem in likelihood estimation with covariance models.}
\usage{
pd.solve(M,sym=TRUE,semi=TRUE,...)

pd.logdet(M,sym=TRUE,semi=TRUE,...)

pd.sqrtm(M,semi=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{A square matrix.}
  \item{sym}{Assume the matrix to be symmetric.}
  \item{semi}{Assume the matrix to only be positive semidefinite (variances can be zero), rather than strictly positive definite (variances must be positive).}
  \item{...}{Optional arguments to other functions, such as \code{\link{qr.solve}}.}
}

\details{If \code{semi=FALSE}, all true variances are assumed to be positive and any numerically estimated variances that fall below machine precision are extrapolated from the numerically estimated variances that fall above machine precision.

Infinite variances can be exactly handled, as long as they are not correlated with finite variances.
}

\value{\code{pd.solve} returns the matrix inverse, \code{pd.logdet} returns the logarithm of the determinant, and \code{pd.sqrtm} returns the square-root matrix.}

%\references{}

\author{ C. H. Fleming. }

%\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{qr.solve}}, \code{\link{det}} }

%\examples{}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
