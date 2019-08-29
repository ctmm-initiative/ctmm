\name{ctmm-package}
\alias{ctmm-package}
\encoding{UTF-8}
\docType{package}
\title{Continuous-time movement modeling}
\description{
Functions for identifying, fitting, and applying continuous-space, continuous-time stochastic movement models to animal tracking data.
The package is described in Calabrese & Fleming (2016) <doi:10.1111/2041-210X.12559> and its methods are based on those introduced in
Fleming & Calabrese et al (2014) <doi:10.1086/675504>,
Fleming et al (2014) <doi:10.1111/2041-210X.12176>,
Fleming et al (2015) <doi:10.1890/14-2010.1>,
Fleming et al (2016) <doi:10.1890/15-1607>,
Péron & Fleming et al (2016) <doi:10.1186/s40462-016-0084-7>,
Fleming & Calabrese (2016) <doi:10.1111/2041-210X.12673>,
Péron et al (2017) <doi:10.1002/ecm.1260>,
Fleming et al (2017) <doi:10.1016/j.ecoinf.2017.04.008>,
Fleming et al (2018) <doi:10.1002/eap.1704>,
and
Winner & Noonan et al (2018) <doi:10.1111/2041-210X.13027>.
}
\details{
\tabular{ll}{
Package: \tab ctmm\cr
Type: \tab Package\cr
Version: \tab 0.5.7 \cr
Date: \tab 2019-08-29 \cr
License: \tab GPL-3 \cr
}

\itemize{
\item{ \href{http://biology.umd.edu/movement.html}{CTMM Initiative} }
\item{ \href{https://CRAN.R-project.org/package=ctmm}{CRAN package} }
\item{ \href{https://github.com/ctmm-initiative/ctmm}{Github project} }
\item{ \href{http://www2.physics.umd.edu/~hfleming/}{Source packages} }
\item{ \href{https://ctmm-initiative.github.io/ctmm/}{Github reference} }
\item{ \href{https://groups.google.com/forum/#!forum/ctmm-user}{Google group} }
\item{ \link{ctmm-FAQ} }
}

}

\author{
Christen H. Fleming and Justin M. Calabrese

Maintainer: Christen H. Fleming <flemingc@si.edu>
}

\references{
C. H. Fleming, J. M. Calabrese, T. Mueller, K.A. Olson, P. Leimgruber, W. F. Fagan.
From fine-scale foraging to home ranges: A semi-variance approach to identifying movement modes across spatiotemporal scales.
\href{http://www.jstor.org/discover/10.1086/675504}{The American Naturalist, 183:5, E154-E167 (2014)}.

C. H. Fleming, J. M. Calabrese, T. Mueller, K. A. Olson, P. Leimgruber, W. F. Fagan.
Non-Markovian maximum likelihood estimation of autocorrelated movement processes.
\href{http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12176/abstract}{Methods in Ecology and Evolution, 5:5 462-472 (2014)}.

C. H. Fleming, Y. Subasi, J. M. Calabrese.
A maximum-entropy description of animal movement.
\href{http://journals.aps.org/pre/abstract/10.1103/PhysRevE.91.032107}{Physical Review E, 91, 032107 (2015)}.

C. H. Fleming, W. F. Fagan, T. Mueller, K. A. Olson, P. Leimgruber, J. M. Calabrese.
Rigorous home-range estimation with movement data: A new autocorrelated kernel-density estimator.
\href{https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/14-2010.1}{Ecology, 96:5, 1182-1188 (2015)}.

J. M. Calabrese, C. H. Fleming.
ctmm: an R package for analyzing animal relocation data as a continuous-time stochastic process.
\href{http://onlinelibrary.wiley.com/wol1/doi/10.1111/2041-210X.12559/abstract}{Methods in Ecology and Evolution, 7:9, 1124-1132 (2016)}.

C. H. Fleming, W. F. Fagan, T. Mueller, K. A. Olson, P. Leimgruber, J. M. Calabrese.
Estimating where and how animals travel: An optimal framework for path reconstruction from autocorrelated tracking data.
\href{https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1890/15-1607.1}{Ecology, 97:3, 576-582 (2016)}.

G. Péron, C. H. Fleming,  R. C. de Paula, J. M. Calabrese.
Uncovering periodic patterns of space use in animal tracking data with periodograms, including a new algorithm for the Lomb-Scargle periodogram and improved randomization tests.
\href{https://movementecologyjournal.biomedcentral.com/articles/10.1186/s40462-016-0084-7}{Movement Ecology, 4:19, DOI:10.1186/s40462-016-0084-7 (2016)}.

C. H. Fleming, J. M. Calabrese.
A new kernel-density estimator for accurate home-range and species-range area estimation.
\href{http://onlinelibrary.wiley.com/wol1/doi/10.1111/2041-210X.12673/abstract}{Methods in Ecology and Evolution, 8:5, 571-579 (2016)}.

G. Péron, C. H. Fleming, R. C. de Paula, N. Mitchell, M. Strohbach, P. Leimgruber, J. M. Calabrese.
Periodic continuous-time movement models uncover behavioral changes of wild canids along anthropization gradients.
\href{http://onlinelibrary.wiley.com/wol1/doi/10.1002/ecm.1260/abstract}{Ecological Monographs 87:3, 442-456 (2017)}

C. H. Fleming, D. Sheldon, E. Gurarie, W. F. Fagan, S. LaPoint, J. M. Calabrese.
Kálmán filters for continuous-time movement models.
\href{http://www.sciencedirect.com/science/article/pii/S1574954117301115}{Ecological Informatics, 40, 8-21 (2017)}.

C. H. Fleming, D. Sheldon, W. F. Fagan, P. Leimgruber, T. Mueller, D. Nandintsetseg, M. J. Noonan, K. A. Olson, E. Setyawan, A. Sianipar, J. M. Calabrese.
Correcting for missing and irregular data in home-range estimation.
\href{https://esajournals.onlinelibrary.wiley.com/doi/10.1002/eap.1704}{Ecological Applications, 28:4, 1003-1010 (2018)}.

K. Winner, M. J. Noonan, C. H. Fleming, K. Olson, T. Mueller, D. Sheldon, J. M. Calabrese.
Statistical inference for home range overlap.
\href{https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13027}{Methods in Ecology and Evolution, DOI:10.1111/2041-210X.13027 (2018)}.
}

%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
%\keyword{ package }

%\note{}

%\seealso{}
%\examples{
%~~ simple examples of the most important functions ~~
%}
