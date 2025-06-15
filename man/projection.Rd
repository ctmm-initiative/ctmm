\name{projection}
\alias{projection}
\alias{projection,telemetry-method}
\alias{projection<-,telemetry-method}
\alias{projection,ctmm-method}
\alias{projection<-,ctmm-method}
\alias{projection,UD-method}
\alias{projection,list-method}
\alias{projection<-,list-method}
\alias{projection,NULL-method}
\alias{median}
\alias{median.telemetry}
\alias{compass}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Projection}
\description{Functions to manipulate the coordinate reference system (CRS) of ctmm objects}

\usage{
\S4method{projection}{telemetry}(x,asText=TRUE)

\S4method{projection}{ctmm}(x,asText=TRUE)

\S4method{projection}{UD}(x,asText=TRUE)

\S4method{projection}{list}(x,asText=TRUE)

\S4method{projection}{NULL}(x,asText=TRUE)

\S4method{projection}{telemetry}(x) <- value

\S4method{projection}{ctmm}(x) <- value

\S4method{projection}{list}(x) <- value

\S3method{median}{telemetry}(x,na.rm=FALSE,...)

compass(loc=NULL,cex=3,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A \code{telemetry}, \code{ctmm}, or \code{UD} object.}
  \item{asText}{ If \code{TRUE}, the projection is returned as text. Otherwise a \code{\link[sp]{CRS}} object is returned.}
  \item{value}{Projection to apply. Can also be a data.frame of longitude-latitude foci.}
  \item{na.rm}{Not used.}
  \item{...}{Arguments passed to \code{\link[Gmedian]{Gmedian}} or \code{\link[graphics]{text}}.}
  \item{loc}{Optional two-dimensional coordinates (in meters) at which to draw a north-facing compass needle.}
  \item{cex}{Relative size of compass.}
}

\details{\code{projection(x)} returns the projection information from ctmm object \code{x}, while \code{projection(x) <- value} applies the projection \code{value} to object \code{x}.
\code{median(x)} returns the ellipsoidal geometric median of a telemetry object.
\code{compass(c(x,y))} plots a north-pointing compass needle at the coordinates \eqn{(x,y)}.}

%\references{}

\author{C. H. Fleming}

%\note{Plotting UTF-8 chracters in a PDF, like the compass needle "➢", requires specifying a compatible font family. For example:
\note{Plotting UTF-8 chracters in a PDF, like the compass needle, requires specifying a compatible font family. For example:

\code{library(ctmm)} \cr
\code{data(buffalo)} \cr
\code{cairo_pdf(file="buffalo.pdf",family="DejaVu Sans")} \cr
\code{plot(buffalo[[1]])} \cr
\code{compass()} \cr
\code{dev.off()}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{as.telemetry}}. }

\examples{
# Load package and data
library(ctmm)
data(buffalo)

# Apply a 1-point projection that preserves North==up
projection(buffalo) <- median(buffalo)
plot(buffalo)
compass()

# Apply a 2-point projection safer for elongated disributions
projection(buffalo) <- median(buffalo,k=2)
# This is the default projection for ctmm
plot(buffalo)
compass()
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
