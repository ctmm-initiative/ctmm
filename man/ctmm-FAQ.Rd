\name{ctmm-FAQ}
\alias{ctmm-FAQ}
\encoding{UTF-8}

\description{Frequently asked questions for the \code{ctmm} package.}

\title{ctmm FAQ}

\details{
\bold{General recommendations:}
\enumerate{
\item{Do not save workspaces between sessions. They can become corrupted.}
\item{Upgrade \code{R} to the latest version and update all of your packages.}
\item{Check the \href{http://www2.physics.umd.edu/~hfleming/}{\code{ctmm} website} to see if there is a beta release that fixes any bugs.}
}
\bold{Specific issues:}
\describe{
\item{Database corruption, missing functions, and other weird errors}{\code{R} likely did not install the package correctly. I.e., some files may have failed to overwrite previous versions. Uninstall \code{ctmm}, restart \code{R}, and install \code{ctmm} again.}
\item{\code{manipulate::isAvailable} is not found}{You probably have an outdated copy of the \code{manipulate} package installed. Update \code{R} to the latest version and then update all of your packages. This seems to happen frequently with the MacOS release of R.}
\item{The \code{manipulate} panel does not popup in \code{\link{variogram.fit}} or \code{\link{zoom}}}{Click the gear icon in the upper-left corner of the plot window.}
\item{There is no save button in \code{\link{variogram.fit}}}{Maximize the plot window and/or increase your screen resolution.}
\item{\code{plot} complains about the datatype or has weird errors}{Namespace collision sometimes occurs from \code{raster} or \code{sp} against \code{move} and \code{ctmm}. Restart \code{R} and only load the \code{ctmm} package.}
\item{Infinite recursion error}{\code{ctmm} has no recursive functions, so I am not exactly sure what causes this error, but it only occurs with certain versions of \code{R} on certain computer architectures. There are several solutions that have worked for people, including restarting \code{R} in a fresh session and updating their software:
\enumerate{
\item{Reboot your computer.}
\item{Increase the allowed number of nested expressions within \code{R} via \code{options(expressions=10000)} or some other large number.}
\item{Try a different computer.}
}}

}
}

\author{Chris H. Fleming}