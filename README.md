ctmm: Continuous-Time Movement Modeling
=======================================

ctmm is an R package for analyzing animal tracking data as a continuous-time stochastic processes. Package features include FFT variogram and Lomb-Scargle periodogram analysis, perturbative REML estimation, Kriging, simulation, and autocorrelated kernel-density home-range estimation.

New user resources:

* [Methods in Ecology and Evolution paper](https://doi.org/10.1111/2041-210X.12559)
* [GitHub reference](https://ctmm-initiative.github.io/ctmm/)
  - [Variogram vignette](https://ctmm-initiative.github.io/ctmm/articles/variogram.html)
  - [AKDE vignette](https://ctmm-initiative.github.io/ctmm/articles/akde.html)
* [AniMove lectures (2022)](https://streaming.uni-konstanz.de/talks-und-events/2022/animove-2022/animove-2022-09-16/)
* [ctmmlearn workshop materials](https://github.com/ctmm-initiative/ctmmlearn)
* [Google group](https://groups.google.com/g/ctmm-user)

More links:

* [The ctmm initiative](https://biology.umd.edu/movement)
* [Movement of Life](https://movementoflife.si.edu/analytical-tools/)
* [GitHub project](https://github.com/ctmm-initiative/ctmm)

Installation
============

Installing the latest stable release from CRAN:
```r
install.packages("ctmm")
```

Installing the latest development release from Github:
```r
remotes::install_github("ctmm-initiative/ctmm")
```
or
```r
devtools::install_github("ctmm-initiative/ctmm")
```
