ctmm: Continuous-Time Movement Modeling
=======================================

ctmm is an R package for analyzing animal tracking data as a continuous-time stochastic processes. Package features include FFT variogram and Lomb-Scargle periodogram analysis, perturbative REML estimation, Kriging, simulation, and autocorrelated kernel-density home-range estimation. First-time users can start with [the Methods in Ecology and Evolution paper](https://doi.org/10.1111/2041-210X.12559) and then move on to the up-to-date vignettes `vignette('variogram')` and `vignette('akde')`.

* [The ctmm initiative](https://biology.umd.edu/movement)
* [Movement of Life](https://movementoflife.si.edu/analytical-tools/)
* [Beta packages](http://www2.physics.umd.edu/~hfleming/)
* [Google group](https://groups.google.com/g/ctmm-user)
* [GitHub project](https://github.com/ctmm-initiative/ctmm)
* [GitHub reference](https://ctmm-initiative.github.io/ctmm/)

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
