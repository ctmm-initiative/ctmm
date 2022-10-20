\name{Unit conversion}
\alias{\%#\%}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert dimensionful quantities to and from SI units}
\description{
This function takes a number in some specified units and converts that number to SI units,
or from SI units to the specified units.
Internally, all \code{ctmm} objects are specified in SI units, and so this is a utility function to facilitate working with \code{ctmm} objects.
}
\usage{
x \%#\% y
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A numeric quantity specified in \code{y} character labeled units, or a character unit label to convert a numeric quantity \code{y} that is specified in SI units.}
  \item{y}{A unit character label for the quantity \code{x} to be converted to SI units, or a numeric quantity in SI units to be converted into unit label \code{x}.}
}

\details{If \code{x} is a number and \code{y} is a character unit label, then \code{x} is converted from units \code{y} to SI units. If \code{x} is a character unit label and \code{y} is a number, then \code{y} is converted from SI units to units \code{x}.

The default non-SI units include the mean solar \code{'day'}, mean synodic \code{'month'} and mean tropical \code{'year'}. These defaults can be changed to conventional calendar units via \code{options(time.units='calendar')}.
}

\value{Returns a numeric in SI units or units specified by character label \code{x}.}

%\references{
%C. H. Fleming and W. F. Fagan and T. Mueller and K. A. Olson and P. Leimgruber and J. M. Calabrese (2015).
%Rigorous home-range estimation with movement data: A new autocorrelated kernel-density estimator.
%\href{http://www.esajournals.org/doi/abs/10.1890/14-2010.1}{Ecology, 96(5), 1182-1188}.

%D. Marcotte. (1996). Fast variogram computation with FFT. \emph{Computers and Geosciences} 22(10), 1175-1186.
%}

\author{ C. H. Fleming. }

%\note{  }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{unit}} }

\examples{
# one yard -> meters
1 \%#\% "yard"

# one meter -> yards
"yard" \%#\% 1

# 1 month -> days
"day" \%#\% 1 \%#\% "month"

# 6 miles per hour -> meters per second
"hour" \%#\% 6 \%#\% "mile"

# the same conversion in one step
6 \%#\% "mph"
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
