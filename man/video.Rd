\name{video}
\alias{video}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Video record animated telemetry objects.}
\description{
Produces an MP4 video file by animating telemetry objects.
}
\usage{
video(x,ext=extent(x),fps=60,dt=NULL,ghost=0,timestamp=FALSE,file="ctmm.mp4",res=720,
      col="red",pch=1,cex=NULL,lwd=1,par.list=list(),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{telemetry} object or list of \code{telemetry} objects.}
  \item{ext}{Plot extent for all frames.}
  \item{fps}{Frames per viewed second.}
  \item{dt}{Tracked time per frame (not per viewed second). By default, the median timestep will be used.}
  \item{ghost}{Timescale over which image retention (ghosting) decays.}
  \item{timestamp}{Display timestamps on title.}
  \item{file}{File name for MP4 file to save. The full path can also be specified. Otherwise the working directory will be used.}
  \item{res}{Pixel resolution for square videos or pixel \code{c(width,height)} for rectangular videos.}
  \item{col}{Color option for telemetry data. Can be an array or list of arrays.}
  \item{pch}{Plotting symbol. Can be an array or list of arrays.}
  \item{cex}{Relative size of plotting symbols. Only used when errors are missing.}
  \item{lwd}{Line widths of \code{telemetry} points.}
  \item{par.list}{List of additional arguments passed to \code{\link{par}} within \code{animate} that do not work outside of \code{animate}, like \code{mar}.}
  \item{\dots}{ Additional options passed to \code{\link{plot.telemetry}}.}
}

\details{This function does not interpolate locations to make smooth animations. For that, please use \code{\link{predict}} or \code{\link{simulate}} outputs instead of a raw tracking data.}

\value{Saves an MP4 file named \code{file} to the working directory.}

%\references{}

\author{C. H. Fleming.}

\note{Further \code{animation} and ffmpeg options can be set via \code{\link[animation]{ani.options}}.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{plot}}, \code{\link{plot.telemetry}}, \code{\link[animation]{ani.options}} }

\examples{\donttest{
# Load package and data
library(ctmm)
data(coati)

# temporary file to store videos for CRAN compliance
FILE <- tempfile("ctmm",fileext=".mp4")
# you will likely want to save your video elsewhere
# the working directory is the default location

# create guess object
GUESS <- ctmm.guess(coati[[2]],interactive=FALSE)
# in general, use ctmm.select instead of ctmm.fit
FIT <- ctmm.fit(coati[[2]],GUESS)

# consider a few hours of consecutive sampling, at 1 minute per frame
t <- seq(coati[[2]]$t[19],coati[[2]]$t[27],by=60)

# tau[velocity] is a natural scale to demonstrate persistance of motion
ghost <- FIT$tau[2]

# predicted locations each minute
PRED <- predict(coati[[2]],FIT,t=t)

# most likely path
video(PRED,error=FALSE,pch=16,ghost=ghost,file=FILE)

# prediction (distribution)
video(PRED,error=3,file=FILE)

# conditional simulations
SIMS <- lapply(1:6,function(i){simulate(coati[[2]],FIT,t=t)})

# random paths
video(SIMS,pch=16,ghost=ghost,file=FILE)
}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
