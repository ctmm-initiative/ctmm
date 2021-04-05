\name{ctmm-package}
\alias{ctmm-package}
\encoding{UTF-8}
\docType{package}
\title{Continuous-time movement modeling}
\description{
Functions for identifying, fitting, and applying continuous-space, continuous-time stochastic movement models to animal tracking data.
The package is described in Calabrese & Fleming (2016) <doi:10.1111/2041-210X.12559> and its models and methods are based on those introduced in
Fleming & Calabrese et al (2014) <doi:10.1086/675504>,
Fleming et al (2014) <doi:10.1111/2041-210X.12176>,
Fleming et al (2015) <doi:10.1103/PhysRevE.91.032107>,
Fleming et al (2015) <doi:10.1890/14-2010.1>,
Fleming et al (2016) <doi:10.1890/15-1607>,
Péron & Fleming et al (2016) <doi:10.1186/s40462-016-0084-7>,
Fleming & Calabrese (2017) <doi:10.1111/2041-210X.12673>,
Péron et al (2017) <doi:10.1002/ecm.1260>,
Fleming et al (2017) <doi:10.1016/j.ecoinf.2017.04.008>,
Fleming et al (2018) <doi:10.1002/eap.1704>,
Winner & Noonan et al (2018) <doi:10.1111/2041-210X.13027>,
Fleming et al (2019) <doi:10.1111/2041-210X.13270>,
Noonan & Fleming et al (2019) <doi:10.1186/s40462-019-0177-1>,
Fleming et al (2020) <doi:10.1101/2020.06.12.130195>,
and
Noonan et al (2020) <doi:10.1101/2020.08.24.261628>.
}
\details{
\tabular{ll}{
Package: \tab ctmm\cr
Type: \tab Package\cr
Version: \tab 0.6.1 \cr
Date: \tab 2021-04-04 \cr
License: \tab GPL-3 \cr
}

\itemize{
\item{ \href{http://biology.umd.edu/movement.html}{CTMM Initiative} }
\item{ \href{https://movementoflife.si.edu/analytical-tools/}{Movement of Life} }
\item{ \href{https://CRAN.R-project.org/package=ctmm}{CRAN package} }
\item{ \href{https://github.com/ctmm-initiative/ctmm}{Github project} }
\item{ \href{http://www2.physics.umd.edu/~hfleming/}{Source packages} }
\item{ \href{https://ctmm-initiative.github.io/ctmm/}{Github reference} }
\item{ \href{https://groups.google.com/g/ctmm-user}{Google group} }
\item{ \link{ctmm-FAQ} }
}

}

\author{
Christen H. Fleming and Justin M. Calabrese

Maintainer: Christen H. Fleming <flemingc@si.edu>
}

\references{
C. H. Fleming, J. M. Calabrese, T. Mueller, K. A. Olson, P. Leimgruber, W. F. Fagan,
``From fine-scale foraging to home ranges: A semi-variance approach to identifying movement modes across spatiotemporal scales'',
The American Naturalist, 183:5, E154-E167 (2014) \doi{10.1086/675504}.

C. H. Fleming, J. M. Calabrese, T. Mueller, K. A. Olson, P. Leimgruber, W. F. Fagan,
``Non-Markovian maximum likelihood estimation of autocorrelated movement processes'',
Methods in Ecology and Evolution, 5:5 462-472 (2014) \doi{10.1111/2041-210X.12176}.

C. H. Fleming, Y. Subaşı, J. M. Calabrese,
``A maximum-entropy description of animal movement'',
Physical Review E, 91, 032107 (2015) \doi{10.1103/PhysRevE.91.032107}.

C. H. Fleming, W. F. Fagan, T. Mueller, K. A. Olson, P. Leimgruber, J. M. Calabrese,
``Rigorous home-range estimation with movement data: A new autocorrelated kernel-density estimator'',
Ecology, 96:5, 1182-1188 (2015) \doi{10.1890/14-2010.1}.

J. M. Calabrese, C. H. Fleming,
``ctmm: an R package for analyzing animal relocation data as a continuous-time stochastic process'',
Methods in Ecology and Evolution, 7:9, 1124-1132 (2016) \doi{10.1111/2041-210X.12559}.

C. H. Fleming, W. F. Fagan, T. Mueller, K. A. Olson, P. Leimgruber, J. M. Calabrese,
``Estimating where and how animals travel: An optimal framework for path reconstruction from autocorrelated tracking data'',
Ecology, 97:3, 576-582 (2016) \doi{10.1890/15-1607.1}.

G. Péron, C. H. Fleming,  R. C. de Paula, J. M. Calabrese,
``Uncovering periodic patterns of space use in animal tracking data with periodograms, including a new algorithm for the Lomb-Scargle periodogram and improved randomization tests'',
Movement Ecology, 4:19 (2016) \doi{10.1186/s40462-016-0084-7}.

C. H. Fleming, J. M. Calabrese,
``A new kernel-density estimator for accurate home-range and species-range area estimation'',
Methods in Ecology and Evolution, 8:5, 571-579 (2017) \doi{10.1111/2041-210X.12673}.

G. Péron, C. H. Fleming, R. C. de Paula, N. Mitchell, M. Strohbach, P. Leimgruber, J. M. Calabrese,
``Periodic continuous-time movement models uncover behavioral changes of wild canids along anthropization gradients'',
Ecological Monographs 87:3, 442-456 (2017) \doi{10.1002/ecm.1260}

C. H. Fleming, D. Sheldon, E. Gurarie, W. F. Fagan, S. LaPoint, J. M. Calabrese,
``Kálmán filters for continuous-time movement models'',
Ecological Informatics, 40, 8-21 (2017) \doi{10.1016/j.ecoinf.2017.04.008}.

C. H. Fleming, D. Sheldon, W. F. Fagan, P. Leimgruber, T. Mueller, D. Nandintsetseg, M. J. Noonan, K. A. Olson, E. Setyawan, A. Sianipar, J. M. Calabrese,
``Correcting for missing and irregular data in home-range estimation'',
Ecological Applications, 28:4, 1003-1010 (2018) \doi{10.1002/eap.1704}.

K. Winner, M. J. Noonan, C. H. Fleming, K. Olson, T. Mueller, D. Sheldon, J. M. Calabrese.
``Statistical inference for home range overlap'',
Methods in Ecology and Evolution, 9:7, 1679-1691 (2018) \doi{10.1111/2041-210X.13027}.

C. H. Fleming, M. J. Noonan, E. P. Medici, J. M. Calabrese,
``Overcoming the challenge of small effective sample sizes in home-range estimation'',
Methods in Ecology and Evolution 10:10, 1679-1689 (2019) \doi{10.1111/2041-210X.13270}.

M. J. Noonan, C. H. Fleming, T. S. Akre, J. Drescher-Lehman, E. Gurarie, A.-L. Harrison, R. Kays, Justin Calabrese,
``Scale-insensitive estimation of speed and distance traveled from animal tracking data'',
Movement Ecology, 7:35 (2019) \doi{10.1186/s40462-019-0177-1}.
}

%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
%\keyword{ package }

%\note{}

%\seealso{}
%\examples{
%~~ simple examples of the most important functions ~~
%}
