0.3.6 2017-02-22
================

  * Better UD plot boxes with new extent methods

  * Reinstated column name preferences in as.telemetry
  
  * as.telemetry faster with fastPOSIXct
  
  * new trace option for ctmm.fit
  
  * more robust CIs for pREML, REML
  
  * chi-square CIs (area, semi-variance, etc.) more robust for DOF<1

0.3.5 2017-02-01
================

  * added a FAQ page to the documentation help("ctmm-FAQ")

  * bugfix in occurrence method for BM & IOU models
  
  * unit conversion can now be disabled in summary with units=FALSE argument
  
  * added trace option to ctmm.select & bandwidth/akde
  
  * improved telemetry error support in summary.ctmm and plot.variogram
  
  * as.telemetry more robust to alternative column label capitalizations
  
  * ctmm.loglike & ctmm.fit more robust when tau_velocity ~ tau_position

  * Kalman filter & smoother upgraded to Joseph form covariance updates

0.3.4 2016-11-28
================

  * weighted AKDE implemented, fast option, covered in vignette

  * overlap arguments & ouput changed/generalized

  * method akde.bandwidth renamed to bandwidth inline with S3 standards

  * predict now returns covariance estimates
  
  * occurrence distributions now exportable

  * AKDE overlap bugfixes
  
  * summary.ctmm now returns correct RMS speed
  
  * bugfix for eccentricity errors

  * variogram CIs fixed for odd dimensions
  
  * variogram.fit can now accept OU models
  
  * periodogram rare index bugfix
  
  * fixed missing lag in dt-argumented variogram
  
  * as.telemetry column identification more robust
  
  * as.telemetry defined for MoveStack objects

0.3.3 2016-09-05
================

  * improved import of 'move' objects

  * preliminary 3D AKDE support, debiased

  * new method predict for ctmm objects
  
  * akde now supports smoothing errors

  * variogram.fit and plot.variogram now support telemetry error

  * UERE fitting now possible simultaneous with tracking data

  * tag.local.identifier now used as backup to individual.local.identifier in as.telemetry
  
  * multiple bug fixes in uere
  
  * res.space fixed in occurrence

0.3.2 2016-05-12
================

  * new function overlap for stationary Gaussian distributions and KDEs

  * new function uere calculates UERE from calibration data
  
  * akde debias argument removes most bias from area estimtes, now default
  
  * akde CIs further improved

  * variogram, periodogram generalized to arbitrary dimensions
  
  * periodic mean function option for ctmm, ctmm.fit, ctmm.select, plot.variogram, summary (not yet documented)
  
  * new method residuals for ctmm objects

  * ctmm.select now only considers likely model modifications
  
  * DOFs now returned in summary

  * new methods [.telemetry, [.variogram, [.periodogram, subset.periodogram
  
  * methods for zoom, raster, writeShapefile now properly assigned to generics

  * new plot.periodogram option max
  
  * new periodogram option res.time (with Lagrange interpolation). Old option res renamed to res.freq.

  * akde res argument is now relative to the bandwidth
  
  * occurrence res.space argument is now relative to the average diffusion

  * plot.telemetry with data error now uses level.UD for error radius instead of one standard deviation
  
  * gridding function for fast=TRUE variogram and periodogram now always fast
  
  * bad location removed from buffalo "Pepper"

0.3.1 2016-02-23
================

  * variogram.fit now stores global variables of any name

  * variogram.fit sliders now use pretty units

  * variogram.fit range argument depreciated in favor of a more general CTMM prototype argument

  * akde UD CIs significantly improved for high quality datasets

  * akde bugfix: subscript out of bounds

  * circulatory model introduced via circle ctmm argument

  * oscillatory CPF model introduced via CPF ctmm argument

  * as.telemetry now imports GPS.HDOP columns with a UERE argument

  * summary now works on arbitrary lists of ctmm objects

  * ctmm.fit now tries to make sense of ML parameters that lie on boundaries
  
  * occurrence() now works when some timesteps are tiny
  
0.3.0 2015-11-26
================

  * new function "occurrence" to estimate occurrence distributions
  
  * "akde" & "occurrence" class objects generalized to "UD" class

  * alpha & alpha.HR arguments simplified and generalized to level & level.UD

  * AKDE= and *.HR= arguments generalized to UD= and *.UD=

  * new basic telemetry error functionality in ctmm, ctmm.fit

  * new function ctmm.select

  * new methods subset.telemetry and subset.variogram

  * fixed a bug in the uncertainty report of uncorrelated processes
  
  * ctmm.fit is now much faster by specifying a reasonable parscale for optim 

  * ctmm.fit now has a backup for when Brent fails

0.2.9 2015-10-13
================

  * fixed a rare condition in ctmm.fit where solve would fail on correlated errors

  * multiscale variogram and mean variogram example in vignette

  * new data example Mongolian gazelle

  * new fast option for periodogram
  
  * improvements in plot.periodogram

  * bugfix in as.telemetry for numeric indentifiers

  * bugfix in dt array option of variogram
  
  * new resolution option and better estimation algorithms in akde
  
  * alpha, alpha.HR, res arguments made consistent across all functions

0.2.8 2015-08-25
================

  * efficiency gains in as.telemetry with multiple animals

  * bugfix in plot.telemetry for multiple Gaussian PDFs

  * bugfix in variogram for rare condition when fast=TRUE

0.2.7 2015-07-27
================

  * CRAN check compliance achieved.
  
  * all methods (plot, mean, summary, simulate) can and must be run without class extensions
  
  * argument names no longer clash with function names and are more explicit about their object class

0.2.6 2015-07-17
================

  * export bugfixes

0.2.5 2015-07-14
================

  * IOU bug fixes in ctmm.fit and plot.variogram
  
0.2.4 2015-06-28
================

  * cleaned up and labeled tau parameter arrays

  * implemented Workaround for when subset demotes S4 objects to S3 objects

  * plot.telemetry now enforces asp=1 even with xlim/ylim arguments

0.2.3 2015-06-19
================

  * new function summary.telemetry

  * bugfix in as.telemetry for data$t
  
  * bugfix in ctmm.loglike for some cases with numeric underflow

  * periodogram and plot.periodogram can now check for spurious periodicities
  
  * minimal support for BM and IOU motion

0.2.2 2015-05-21
================

  * new functions periodogram, plot.periodogram

0.2.1 2015-05-08
================

  * new function SpatialPoints.telemetry returns SpatialPoints object from telemetry data

  * new function SpatialPolygonsDataFrame.akde returns akde home-range contour SpatialPolygons objects from akde object

  * new function writeShapefile.akde writes akde home-range contours to ESRI shapefile

  * new function raster.akde returns akde pdf raster object

  * new function summary.akde returns HR area of AKDE

  * fixed bad CI in plot.telemetry model option

  * as.telemetry now takes a timezone argument for as.POSIXct and defaults to UTC

  * telemetry, ctmm, and akde objects now have idenification and projection information slotted, with consistent naming throughout

0.2.0 2015-04-27
================

  * vignettes "variogram" and "akde"

  * new function as.telemetry imports MoveBank formatted csv data and returns telemetry objects

  * new function variogram.zoom plots a variogram with zoom slider

  * variogram.fit and variogram.zoom default to a logarithmic-scale zoom slider, which requires much less fiddling

  * plot.variogram now takes multiple variogram, model, and color options

  * plot.telemetry now takes multiple data, model, akde, and color options

  * plot.telemetry can now make probability density plots for both Gaussian model and akde data

  * ctmm.fit no longer screws up results with initial sigma guesstimates. ML parameter estimates now match closely with published Mathematica results. CIs are improved.

  * ks-package was producing incorrect home-range contours and has been replaced with custom code. ML home ranges now match published Mathematica results. CIs should be improved.