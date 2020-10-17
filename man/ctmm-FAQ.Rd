\name{ctmm-FAQ}
\alias{ctmm-FAQ}
\alias{ctmm-faq}
\encoding{UTF-8}

\description{Frequently asked questions for the \code{ctmm} package.}

\title{ctmm FAQ}

\details{
\bold{General recommendations}
\enumerate{
\item{Work through the vignettes \code{vignette("variogram")} and \code{vignette("akde")}. Also, see the \code{help} file for the method of interest, and its \code{example}.}
\item{Do not save workspaces between sessions. They will become corrupted over time. In RStudio, go to \code{Tools: Global Options: Workspace}, uncheck ``Restore'' and set ``Save'' to ``Never''.}
\item{Never edit or save your CSV in Microsoft Excel. The dates will be reformatted incorrectly and inconsistently.}
\item{If using Windows, make sure to have the suggested version of ``Rtools'' installed. If using MacOS, make sure to have ``Xcode'' installed. If using Ubuntu, make sure to have ``build-essential'' installed. Otherwise, you can sometimes run into problems when trying to update packages.}
\item{Upgrade \code{R} to the latest version and update all of your packages.}
\item{The development build can be installed via \code{devtools::install_github("ctmm-initiative/ctmm")}.}
\item{Stable beta releases between the CRAN release are published \href{http://www2.physics.umd.edu/~hfleming/}{here} on request.}
\item{The \href{https://groups.google.com/forum/#!forum/ctmm-user}{\code{ctmm} user's group} is a good place to find and ask for help.}
\item{Bug reports and feature requests can be raised at the \href{https://github.com/ctmm-initiative/ctmm}{Github project page}.}
}

\bold{Help installing packages on Linux}

{These are the packages I needed in Ubuntu:

\code{sudo apt install fftw3 fftw3-dev libgdal1-dev libproj-dev libgeos-dev libgsl-dev libv8-dev r-base-core}}


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

\bold{\code{plot} complains about the datatype or has weird errors}

{Namespace collision sometimes occurs from \code{raster} or \code{sp} against \code{move} and \code{ctmm}. Restart \code{R} and only load the \code{ctmm} package.}

\bold{North is no longer up after importing data}

{The default projection in \code{ctmm} does not preserve the direction of North, but better preserves distances for elongated distributions. See the \code{projection} argument in \code{\link{as.telemetry}} and the example in \code{\link[ctmm]{projection}}. The \code{\link{compass}} function is also useful for pointing north.}

\bold{\code{\link{ctmm.guess}} has no save button}

{Maximize the plot window and/or increase your screen resolution.}

\bold{\code{manipulate} panel does not popup in \code{\link{ctmm.guess}} or \code{\link{zoom}}}

{Click the gear icon in the upper-left corner of the plot window.}

\bold{Gear icon missing in \code{\link{ctmm.guess}} or \code{\link{zoom}}}

{Recent versions of \code{manipulate} and/or RStudio seem to have some issues. Sometimes the gear icon does not render unless you re-run the function 2-5 times.}

\bold{\code{manipulate::isAvailable} is not found}

{You probably have an outdated copy of the \code{manipulate} package installed. Update \code{R} to the latest version and then update all of your packages. This seems to happen frequently with the MacOS release of \code{R}.}
}


\author{C. H. Fleming}
