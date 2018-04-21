\name{ctmm-FAQ}
\alias{ctmm-FAQ}
\alias{ctmm-faq}
\encoding{UTF-8}

\description{Frequently asked questions for the \code{ctmm} package.}

\title{ctmm FAQ}

\details{
\bold{General recommendations}
\enumerate{
\item{Never edit your CSV in Microsoft Excel. The dates will be reformatted incorrectly.}
\item{Do not save workspaces between sessions. They can become corrupted.}
\item{Upgrade \code{R} to the latest version and update all of your packages.}
\item{Work through the vignettes \code{vignette("variogram")} and \code{vignette("akde")}.}
\item{Check the \href{http://www2.physics.umd.edu/~hfleming/}{\code{ctmm} website} to see if there is a beta release that fixes any bugs.}
\item{Check the \href{https://groups.google.com/forum/#!forum/ctmm-user}{\code{ctmm} user's group}.}
}

\bold{Help installing packages on Linux}
\itemize{
\item{These are the packages I needed in Ubuntu: \code{sudo apt install fftw3 fftw3-dev libgdal1-dev libproj-dev libgeos-dev r-base-core}}
\item{See Xianghui Dong's installation tips \href{https://github.com/ctmm-initiative/ctmm-webapp/blob/master/README-install.md}{here}.}
}

\bold{\code{ctmm.select} only returns one model with \code{verbose=TRUE}}

{The newer version of \code{ctmm.select} is much more intelligent, necessitated by the exploding number of new model permutations. Any simpler models were judged to be so unlikely that they were not attempted. You can set \code{level=1} to force alternative model fitting, but the additional models should lose badly in model selection (please contact me if I am wrong).}


\bold{\code{as.telemetry} reports abnormal sampling intervals and speeds}

{Make sure that you have the correct \code{timezone} and \code{timeformat} arguments specified. Also, see \code{\link{outlie}}.}

\bold{rdb database corruption, "could not find function", "cannot coerce class", and other weird errors}

{\code{R} might not have installed or loaded the package correctly---e.g., some files may have failed to overwrite previous versions---or the workspace/session might be corrupted. Uninstall \code{ctmm}, restart \code{R} without saving the workspace/session, and install \code{ctmm} again.}


\bold{Infinite recursion and stack overflow errors}

{\code{ctmm} has no recursive functions, so I am not exactly sure what causes this error, but it only occurs with certain versions of \code{R} on certain computer architectures. There are several solutions that have worked for people, including restarting \code{R} in a fresh session and updating their software. Alternatively:
\enumerate{
\item{Reboot your computer.}
\item{Increase the allowed number of nested expressions within \code{R} via \code{options(expressions=10000)} or some other large number.}
\item{Try a different computer.}
}}


\bold{\code{manipulate::isAvailable} is not found}

{You probably have an outdated copy of the \code{manipulate} package installed. Update \code{R} to the latest version and then update all of your packages. This seems to happen frequently with the MacOS release of \code{R}.}


\bold{\code{manipulate} panel does not popup in \code{\link{variogram.fit}} or \code{\link{zoom}}}

{Click the gear icon in the upper-left corner of the plot window.}


\bold{\code{plot} complains about the datatype or has weird errors}

{Namespace collision sometimes occurs from \code{raster} or \code{sp} against \code{move} and \code{ctmm}. Restart \code{R} and only load the \code{ctmm} package.}


\bold{\code{\link{variogram.fit}} has no save button}

{Maximize the plot window and/or increase your screen resolution.}

}


\author{C. H. Fleming}
