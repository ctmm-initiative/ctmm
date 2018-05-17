\name{periodogram}
\alias{periodogram}
\alias{plot.periodogram}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the Lomb-Scargle periodogram of animal-tracking data}
\description{
This function calculates isotropic Lomb-Scargle periodogram (LSP, Scargle, 1982) from a telemetry object. One of two algorithms is used. The slow \eqn{O(n^2)} algorithm vectorizes the exact relations of Scargle (1982), while the fast \eqn{O(n \log n)} algorithm uses the FFT method described in Péron & Fleming et al (2016). The latter method is exact if the data are evenly scheduled, permitting gaps, and otherwise it can be made arbitrarily precise via the \code{res.time} option.
}
\usage{
periodogram(data,CTMM=NULL,dt=NULL,res.freq=1,res.time=1,fast=NULL,axes=c("x","y"))

\S3method{plot}{periodogram}(x,max=FALSE,diagnostic=FALSE,col="black",transparency=0.25,grid=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ \code{telemetry} data object or list of such objects. }
  \item{CTMM}{ An optional \code{ctmm} model object for specifying the mean. }
  \item{dt}{Sampling interval for frequency cutoff.}
  \item{res.freq}{Multiplier to inflate the frequency resolution.}
  \item{res.time}{Integer multiplier to inflate the temporal resolution. Useful when \code{fast}>0 and the sampling rate is variable.}
  \item{fast}{Use the exact algorithm if FALSE, the FFT algorithm if TRUE, and further inflate the frequency resolution to a power of two sample size if \code{fast=2}.}
  \item{axes}{Array of axes to calculate an average (isotropic) variogram for.}
  \item{x}{Output object of \code{periodogram}.}
  \item{max}{Plot only the local maxima of the periodogram. Use only with \code{res>1}.}
  \item{diagnostic}{Plot the sampling schedule's periodogram to check for spurious periodicities.}
  \item{col}{Color of periodogram.}
  \item{transparency}{Adds transparency to clustered data if greater than zero. Should be less than one.}
  \item{grid}{Whether or not to plot gridlines at common periodicities.}
  \item{...}{Optional arguments fed to \code{\link{plot}}.}
}
\details{
If no \code{dt} is specified, the median sampling interval is used. This is typically a good assumption for most data, even when there are gaps and this choice corresponds to the discrete Fourier transform (DFT) periodogram for evenly-sampled data.
%If no \code{T} is specified, the full sampling period is used.

At default resolution the frequency grid interval is given by \code{1/(2*(range(data$t)+dt))} and the frequency cutoff is given by \code{1/(2*dt)}, both in accordance with the DFT periodogram. Increasing \code{res.freq} beyond \code{res.freq=1} will make for a smooth periodogram, but sequential frequencies will be highly correlated. The \code{max=TRUE} option to \code{plot.periodogram} may be useful for \code{res.freq>1}. Increasing \code{res.time} beyond \code{res.time=1} is helpful if there is variability in the sampling rate and \code{fast>0}.

If a \code{CTMM} argument is provided, the ML mean will be detrended from the data prior to calculating the periodogram. Otherwise, the sample mean will be detrended.

If a list of \code{telemetry} objects are fed into \code{periodogram}, then a mean \code{periodogram} object will be returned with the default \code{dt} and base frequency resolution selected on a worst case basis according to the method described by Péron & Fleming et al (2016).
}

\value{
Returns a periodogram object (class \code{periodogram}) which is a dataframe containing the frequency, \code{f} and the Lomb-Scargle periodogram at that frequency, \code{LSP}.
}
\references{
J. D. Scargle.
Studies in astronomical time-series analysis. II. Statistical aspects of spectral analysis of unevenly-sampled data.
The Astrophysical Journal, 263, 835-853 (1952).

G. Péron, C. H. Fleming,  R. C. de Paula, J. M. Calabrese.
Uncovering periodic patterns of space use in animal tracking data with periodograms, including a new algorithm for the Lomb-Scargle periodogram and improved randomization tests.
\href{https://movementecologyjournal.biomedcentral.com/articles/10.1186/s40462-016-0084-7}{Movement Ecology, 4:19, DOI:10.1186/s40462-016-0084-7 (2016)}.
}

\author{ C. H. Fleming and G. Péron }

\note{ The LSP is totally inappropriate if you in any way alter the sampling rate within the dataset. Stick with variograms in that case. There is a \code{diagnostic} option in \code{plot.periodogram} that can check for spurious periodicities that result from an autocorrelated sampling schedule. This plot will not contain any periodicities if the LSP is appropriate.

\code{res.time>1} relies on Lagrange interpolation of the sinusoids (not the data), which can suffer from Runge's phenomena. \code{periodogram} tests for an invalid result and can fail with an error message. For whatever reason, this more frequently seems to happen when \code{res.time=3}. }

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{ \code{\link{plot.periodogram}}. }

\examples{
#Load package and data
library(ctmm)
data(wolf)

#Extract movement data for a single animal
Tay <- wolf$Tay

#Calculate periodogram (fast==2 for a speedy example)
#There is some variability in the sampling frequency, so we increase res.time
LSP <- periodogram(Tay,fast=2,res.time=2)

#Plot the periodogram
plot(LSP,max=TRUE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
