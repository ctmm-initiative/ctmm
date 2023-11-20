ctmm 1.2.1 (2023-11-19)
================
* periodic mean functions now estimate frequency
* periodic mean summary now uses beta CIs
* new mean.ctmm() argument formula for functional response estimation
* new grid argument dr.fn in akde(), occurrence(), pkde()
* as.telemetry() now imports ATLAS error ellipses
* new ctmm.boot() argument clamp
* mean() of occurrence() now time weighted by default
* new plot method: plot.ctmm()
* plot() argument 'col.DF' renamed to 'col.UD'
* improvements DOF[area] calculations in mean() and pkde()
* improvement to optimizer() when initial guess is on a boundary and hessian is bad
* bugfix in rsf.select with functions in formulas
* bugfix in raster factor expansion with more than one raster factor
* bugfix in mean.ctmm() potential endless loop
* bugfix in normal meta-analysis REML correction

ctmm 1.2.0 (2023-09-22)
================
* new function names: cde() and encounter() replacing encounter() and rates()
* new functions rsf.select(), intensity()
* new functions sdm.fit(), sdm.select()
* new function writeVector(), depreciating function writeShapefile()
* new function funnel() for funnel plots
* new function midpoint()
* new population covariance models and improved model selection in mean.ctmm()
* new argument 'sqrt' in distance()
* new argument 'dt.hot' in as.telemetry()
* new argument 'variable' in Log()
* new argument 'compute' in ctmm.loglike()
* new argument 't' in proximity()
* as.telemetry() now supports GBIF format data
* as.telemetry() datum argument now works on UTM import, and is no longer to a be a complete PROJ string
* as.telemetry() timeformat='auto' now default
* as.telemetry(), plot.telemetry(), rsf.fit() updated from sp to sf transforms
* distance() can now take location arguments
* plot.telemetry() col.DF & col.level arguments can now be color() lists
* suitability() now produces a raster stack corresponding to the CIs
* suitability() on population RSFs now outputs the population suitability
* suitability() extrapolation disabled
* bugfix in tbind for conflicting location classes
* bugfix in suitability()
* bugfix in distance() method="Euclidean", debias=TRUE
* bugfix in rates() debias=TRUE
* bugfix in summary() of population mean location DOF
* bugfix in distances() for 0/0
* bugfix in UD polygon export for tiny areas
* as.telemetry() UTM import updated to new PROJ specification
* mean.ctmm() improved convergence, numerical stability, and covariance selection
* meta() stability improvements for tiny DOF estimates, and OUf support
* overlap() and meta() can now extract object names
* pkde(...) -> akde(...) -> bandwidth(...) -> mean(...) arguments now passed
* rsf.fit() AICc formula improved

ctmm 1.1.0 (2022-11-03)
================
* new function pkde() for population kernel density estimates
* new functions difference(), distances(), proximity() for estimating distances between individuals
* new functions Log(), Exp() to log transform parameter estimates and their uncertainties for meta-analytic regression
* new functions dimfig(), sigfig() to represent quantities with concise units and significant digits
* new argument 'sample' in mean()
* new argument 'interpolate' in rsf.fit()
* new arguments 'xlim', 'ylim' to plot.outlie()
* numerical stability improvements in rsf.fit optimization and hessian calculations
* numerical convergence improvements in location error fitting
* numerical convergence improvements in AKDE weight optimization
* plot.telemetry() can now subset and reproject rasters
* bugfix in sp::polygon derived areas (used since v1.0.0 for summary, plot, meta)
* bugfix in agde(), suitability(), akde() when reprojecting onto the same raster
* bugfix in mean() when averaging isotropic and anisotropic models together
* bugfix in speeds() without telemetry object
* bugfix in cluster() with 0/0 bias correction error
* bugfix in occurrence() with multiple error classes
* bugfix in chi dof computation
* bugfix in outlie() for error ellipses
* summary() now works on mean.ctmm() outputs from different input model structures (OUF & OUO)
* fixed log-chi^2 bias correction in mean.ctmm()

ctmm 1.0.0 (2022-07-07)
================
* new function rsf.fit() to fit integrated resource selection functions (iRSFs) with autocorrelation-adjusted weighted likelihood
* new function mean.ctmm() to calculate population average movement models
* new function revisitation() to calculate the distribution of revisitations
* new function npr() to calculate non-parametric spatial regressions
* new function agde() to calculate autocorrelated Gaussian distribution estimates, with RSF support
* new function suitability() to calculate suitability rasters from RSF fit objects
* new function rates() to calculate relative encounter rates
* new function dt.plot() to inspect sampling intervals
* akde() and occurrence() now support RSF-informed kernels and boundary-respecting kernels
* summary.ctmm() now outputs diffusion rate estimates
* new argument variable for meta() to estimate population diffusion rates, mean speeds, and autocorrelation timescales
* new arguments R and SP in plot.telemetry() and plot.UD() for plotting raster and shapefile base layers
* new option method="Encounter" in overlap()
* mean.UD() now propagates uncertainties
* mean.UD() now functions on occurrence distributions
* new convex argument to UD summary(), plot(), and export functions
* plot() and raster() now work on 3D UDs
* plot.outlie() now works on lists of outlie objects
* speed() output now includes DOF estimate for use with meta()
* tbind() now works correctly with different projections and calibrations
* %#% unit conversion operator can now interpret products and ratios
* summary() timescale confidence intervals are now gamma/inverse-gamma more inline with meta()
* progress bar added to optimizer() when trace=1
* bugfix in IID area CIs
* bugfix in ctmm.loglike() when fitting multiple error classes, where some are zero
* bugfix in ctmm.boot() when bias estimate exceeds variance parameter
* bugfixes in 3D akde()
* bugfix in time gridding code when dt is coarse
* bugfix in SpatialPoints.telemetry for single individuals

ctmm 0.6.1 (2021-07-26)
================
* ctmm.fit() can now fit multiple UERE parameters and update uncertain calibration parameter estimates
* new function cluster()
* new function video()
* new function as.sf()
* new function tbind()
* new argument VMM in simulate(), predict()
* new argument timeformat="auto" in as.telemetry() 
* new argument verbose in meta()
* uere()<- can now assign posterior/updated error estimates from ctmm model objects
* bugfix in ctmm.loglike() for circle!=0 and REML
* bugfixes in optimzer()
* bugfix in ctmm.fit() for 1D processes
* bugfix in variogram.fit() for 1D processes
* bugfixes in simulate(), predict for 1D processes
* bugfix in ctmm.fit() with zero variance models
* bugfix in meta() colors when sort=TRUE
* bugfixes in ctmm.guess(), ctmm.fit(), speed() for tiny amounts of data
* bugfixes in occurrence(), Kalman smoother for IOU process
* ctmm.select() now stores IC and MSPE information for summary()
* extent objects now include missing columns
* extent object longitudes can now cross the international date line
  
ctmm 0.6.0 (2021-01-08)
================
* new function meta() for meta-analysis of home-range areas
* new function encounter() for the conditional distribution of encounters (CDE)
* new function distance() to calculate square Bhattacharyya, Mahalanobis, and Euclidean distances
* new function compass() to plot a north-pointing compass
* new argument 't' in speed()
* new argument 'axes' in outlie()
* as.telemetry() now accepts most tibble objects
* akde() on multiple individuals is now more memory efficient
* bugfix in ctmm.fit() for IOU model
* bugfix in occurrence() with repeated timestamps
* bugfix in summary.ctmm() rowname droped for single parameter CIs
* bugfix in outlie() with list input
* bugfixes in plot.outlie with zero error
* bugfix in variogram() with res>1 and CI="Gauss"
* bugfix in ctmm.select() if stepping OU->OUf->OUF
* bugfix in as.telemetry() for Move objects with empty idData slot
* bugfix in as.telemetry(), median() when importing single location estimate
* bugfix in plot.telemery() with add=TRUE and non-SI units
* bugfix in speed() for ctmm objects (no data), where CIs were incorrect
* bugfix in median() with >=50% repeating observations
* bugfix in summary() for periodic models with tau[velocity]==0
* bugfix in occurrence() for PDclamp() error
* bugfix in ctmm.select() giving incorrect model names when run indirectly
* bugfix in occurrence() with IID autocorrelation model
* workaround in export functions where sp objects change timezones
* workaround in as.telemetry() when Move idData() names are dropped
* workaround in plot.UD() when image() has alpha overflow
* improvements to akde(), occurrence() grid argument when incomplete
* improvements to overlap() Wishart approximation in bias correction
* improvements to cleave()

ctmm 0.5.10 (2020-05-04)
================
* as.telemetry() location class code improved
* as.telemetry() message for marked outliers
* jaguar data in sync with ctmmweb

ctmm 0.5.9 (2020-03-23)
================
* new argument CI="Gauss" in variogram()
* new argument weights in mean.UD()
* new argument datum in as.telemetry() -- input and ouput datums can now differ
* new data 'jaguar'
* bugfix in ctmm.select() for infinte loop
* improvements in ctmm.select, ctmm.loglike for collapsing variance/error estimates
* rewrite of optimizer's line search to be more exact & reliable
* improvements in optimizer for degenerate likelihood surfaces
* improvements in optimization for bad covariance estimates---fit object structure changed
* bugfix in uere.fit with multiple location classes in different orders
* bugfix in speed when fast=FALSE and sampled models lose features
* bugfix in IID pREML CIs
* bugfix in ctmm.guess with large errors causing eigen() to fail
* bugfix in optimizer expansion search step size not increasing
* bugfix in as.telemetry() -- MoveStack objects are given a common projection if not projected

ctmm 0.5.8 (2019-12-09)
================
* improvements to ctmm.select() stepwise selection, especially with error and/or circulation
* improvements to ctmm.fit() for nearly linear home ranges
* improvements to %#% operator -- units of speed supported
* bugfix in ctmm.loglike() for BM/IOU models with error
* new argument units in plot.outlie()
* new options(time.units='mean') and options(time.units='calendar') for %#% operator and display units
* ctmm.select() no longer warns when model features are not supported (ctmm.fit does)
* compatibility fix for R version 4
  
ctmm 0.5.7 (2019-10-06)
================
* new function optimizer()
* new function SpatialPolygonsDataFrame.telemetry() for location estimate error circles/ellipses
* 'pNewton' now the default optimization method
* 'pHREML' now the default estimator & all CI names updated
* grid argument now supported in akde and occurrence methods
* outlie() output now includes CIs with plot method
* error-adjusted variogram implemented when fast=FALSE
* LOOCV now supported in ctmm.select(), summary()
* new buffer argument in occurrence()
* head(), tail() methods for telemetry objects
* str() method for ctmm objects
* new data object 'pelican'
* SpatialPointsDataFrame now includes timestamp
* uere(data) <- numeric now overrides all location classes
* improved support for ARGOS-GPS hybrid data
* missing DOP values now correctly treated as separate location class
* bugfix in conditional simulations with dt argument
* bugfix in plot.UD gridlines
* bugfix in as.telemetry timeout argument when datasets lack timed-out values
* stability fixes in ctmm.fit() for BM/IOU models
* further stability enhancements in ctmm.loglike() and optimizer
* bugfix in simultaneously fit RMS UERE CIs
* AICc formulas fixed for tiny n
* reduced Z^2 now exactly normalized in UERE objects
* minor enhancements to cleave() function
* as.telemetry() no longer automatically calibrates e-obs errors (inconsistent with newer devices)
* as.telemetry() no longer complains on reverse-time-ordered files

ctmm 0.5.6 (2019-05-14)
================
* new functions lasso, marquee, and cleave
* new functions annotate and color
* summary can now compare joint UERE objects to lists of individualized UERE objects
* support for UTM locations in as.telemetry
* support for GPS-ARGOS hybrid data in as.telemetry & uere.fit
* new plot option ext for extent objects
* increased numerical precision in ctmm.loglike for 0 < dt << tau, including the limit OU/OUF -> BM/OU
* BM/IOU model likelihoods are now exact limits of OU/OUF likelihoods modulo a constant
* covariance matrices can now take arbitrary eccentricty and scale
* ctmm.boot new argument iterate=FALSE and bugfixes for iterate=TRUE
* ctmm.boot now debiases the covariance matrix directly (linearly)
* occurrence default dt.max & cor.min arguments now tighter
* periodogram functionality restored for one-dimensional data
* bugfix in IID ctmm.fit with elliptical errors
* bugfix in occurrence when projection origin is far from the mean location
* bugfix in akde.list where location errors were not smoothed
* bugfix in ctmm.guess/variogram.fit for BM/IOU models
* bugfix in speed for IOU models
* e-obs calibration cross checked and fixed
* ctmm.loglike now returns -Inf when movement and error variance are zero
* stability improvements to base R optimizer usage
* bugfix in mark.rm argument of as.telemetry
* cores option added to ctmm.select
* only physical cores now counted by cores arguments
* cores option now used in Windows when appropriate
* improvements to speed, speeds, ctmm.select for short tracks of data

ctmm 0.5.5 (2019-02-11)
================
* bugfix in summary where timescale CIs were always (0,Inf)
* ctmm.select default now level=1
* summary on UERE lists now works with more than one axis
* R dependency increased to >=3.5 for parallel functions
  
ctmm 0.5.4 (2019-02-07)
================
* bugfix in ctmm.select where OU was not considered over the new OUO/OUf models introduced in v0.5.3
* bugfix in ctmm.boot for heteroskedastic errors
* multiplicative option depreciated from ctmm.boot

ctmm 0.5.3 (2019-01-29)
================
* oscillatory (and critically damped) OUO/OUf models now supported, starting with omega option of ctmm()
* summary() now works on lists of UERE objects for error model selection
* MSPE slots & arguments restructured and fully utilized in both summary and ctmm.select
* new method speeds() for estimating instantaneous speeds
* speed() more efficient on very coarse data, slightly improved CIs
* new complete argument in simulate() and predict() to calculate timestamps and geographic coordinates
* now avoiding fastPOSIXct timezone and epoch issues in as.telemetry
* outlie() now works on lists of telemetry objects
* bugfixes in overlap() CIs
* overlap() now robust to bad model fits
* new as.telemetry() argument mark.rm to delete marked outliers
* bugfix in predict() & occurrence() where eccentricity was dropped from covariances
* projection information in Move & MoveStack objects now preserved if possible
* identities preserved with newer MoveStack objects
* ctmm.boot() better handles parameter estimation near boundaries
* e-obs data with missing error/speed/altitude now importing correctly in as.telemetry
* correlogram plots now cap estimates to appropriate range
* beta optimizer now more aggressive in searching along boundaries
* bugfix in ctmm.fit with duplicate timestamps and IID processes without error
* bugfix in ctmm.select with pREML & error
* summary() on telemetry lists no longer fails on length-1 timeseries
* years updated to tropical years and calendar days updated to stellar days
  
ctmm 0.5.2 (2018-09-10)
================
* location classes (multiple UEREs) now supported by uere.fit() and uere()<-
* uere() forked into separate uere() and uere.fit() methods
* AICc slot included in UERE objects for error model selection
* overlap() telemetry and CTMM arguments depreciated
* fixed bug in as.telemetry() when importing ARGOS error ellipses
* e-obs error calibration updated
* numerical stability increased in ctmm.fit when distance scales are extreme

ctmm 0.5.1 (2018-08-06)
================
* Units of measurement down to microns and microseconds now supported
* ctmm.select() now builds up autocovariance features stepwise to help with fit convergence
* residuals() can now be calculated from (calibrated) calibration data---diagnostic argument removed from uere()
* summary.ctmm() now returns DOF[speed] information on individuals
* MSPE of ctmm objects was previously w.r.t. in-sample times and is now time averaged
* summary.list.ctmm() now returns MSPE when useful
* new speed() argument robust for coarse data
* options multiplicative & robust added to ctmm.boot to help with parameters near boundaries
* E-OBS errors adjusted by empirical results of Scott LaPoint's calibration data
* Telonics Gen4 errors estimates imported with results of Patricia Medici's calibration data --- Quick Fixes not yet fully supported
* fixed critical bug in speed()
* fixed bug in as.telemetry with projection argument
* fixed bugs in ctmm.loglike when !isotropic && error && circle
* fixed bug in emulate when fast=FALSE and error=TRUE
* fixed bug in new variogram error calculations (v0.5.0) used for plotting
* simultaneously fitted UERE's from ctmm slot "error" can now be assigned to data for plotting
  
ctmm 0.5.0 (2018-05-15)
================
* Extensive re-write of the Kalman filter & smoother, now supporting an arbitrary number of spatial dimensions, necessary for ARGOS error ellipse support. (Previously, all multi-dimensional problems were transformed into multiple one-dimensional problems.) Many new models will be supported going forward, based on the v0.5.0 code.
* telemetry error vignette "error"
* ARGOS error ellipse support in ctmm.fit() and simulate()
* plotted variogram errors now estimated from HDOP and no longer assumed to be homoskedastic
* as.telemetry() default projections now use robust ellipsoidal statistics
* new median.telemetry() method for help with projecting data
* (anisotropic & circulation & error) models now exact with 2D Kalman filter & smoother
* simulate() & predict() velocities now correct with mean="periodic"
* units argument in speed()
* REML and related methods fixed from 0.4.X 1/2 bug
* ctmm.loglike COV[mu] bugfix for circular error & elliptical movement
* summary() rotation % bugfix with circle=TRUE
* parameter boundary bugfix in ctmm.fit() and ctmm.loglike()
* fixed bandwidth() bug when weights=TRUE on IID process
* variogram.fit() manipulate more appropriate with calibrated errors
* fixed bug in plot.variogram for isotropic model fits
* fixed bug in ctmm.fit with fitted errors and any(diff(t)==0)
* fixed bug in plot.variogram() from stats::qchisq() with k<<1

ctmm 0.4.2 (2018-02-12)
================
* new speed() method
* new ctmm.boot() method
* new outlie() method
* new export functionality for telemetry class
* overlap debias=TRUE option (approximate)
* pHREML, pREML, HREML ctmm.fit methods implemented and documented
* IID pREML & REML AICc values implemented
* MSPE values implemented
* new uere()<- assignment method
* velocity esimtates now included in predict() [fitting one model to multiple behaviors can result in wildly optimistic confidence intervals]
* velocities now included in simulate()
* simulate precompute option
* as.telemetry drop=TRUE option
* as.telemetry will no longer drop individuals with missing data columns
* as.telemetry will try to approximate DOP values
* as.telemetry imports velocity vectors
* as.telemetry default projection orientation now robust with GmedianCOV
* plot.UD resolution grid less obnoxious, NA/FALSE contour label option
* plot.telemetry error=0:3 options for data with recorded error circles/ellipses
* plot.telemetry velocity=TRUE option for data with recorded velocities
* plot.variogram bugfixes with telemetry errors
* fixed AIC bug in new parameterization code (0.4.0-0.4.1) where isotropic=TRUE model would never be selected
* fixed rare endless loop in akde/bandwidth with weights=TRUE
* outlier removed from buffalo$Cilla

ctmm 0.4.1 (2017-08-30)
================
* projection method for ctmm objects

ctmm 0.4.0 (2017-08-29)
================
* periodigram vignette
* new utility function %#% for unit conversions
* new model-fit sampling function "emulate"
* summary now works on lists of telemetry objects
* new extent method for variogram objects
* bugfixes in plot.variogram with fit UERE, tau==0
* bugfixes with ctmm.fit/select/summary near boundaries
* resetting Polak–Ribiere formula in weighted AKDE conjugate gradient routine
* read.table fallback in as.telmetry
* R 3.4 compatibility fixes
* various improvements to plot.variogram
* plot.UD & export can now accept multiple level.UD values
* increased numerical precision in ctmm.loglike
* SI speeds & diffusion fixed with units=FALSE

ctmm 0.3.6 (2017-04-23)
================
* AICc formulas updated from univariate to multivariate
* ctmm.select more aggressive on small sample sizes where AICc >> AIC
* new residuals and correlogram functions
* ctmm.fit now has unified options controling optimization & differentiation
* ctmm.fit Hessian and pREML calculations 2x faster
* new writeRaster method for UD objects
* better UD plot boxes with new extent methods
* variogram fast=TRUE less biased for irregular data with new res>1 option
* variogram fast=FALSE more robust to irregularity
* akde() can now handle duplicate times (with an error model)
* plot.variogram bugfix for fixed error models [still not quite correct]
* Column name preferences in as.telemetry
* as.telemetry faster with fread & fastPOSIXct
* new trace option for ctmm.fit
* new labels option for plot.UD
* more robust CIs for pREML, REML
* chi-square CIs (area, semi-variance, etc.) more robust when DOF<1

ctmm 0.3.5 (2017-02-01)
================
* added a FAQ page to the documentation help("ctmm-FAQ")
* bugfix in occurrence method for BM & IOU models
* unit conversion can now be disabled in summary with units=FALSE argument
* added trace option to ctmm.select & bandwidth/akde
* improved telemetry error support in summary.ctmm and plot.variogram
* as.telemetry more robust to alternative column label capitalizations
* ctmm.loglike & ctmm.fit more robust when tau_velocity ~ tau_position
* Kalman filter & smoother upgraded to Joseph form covariance updates

ctmm 0.3.4 (2016-11-28)
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

ctmm 0.3.3 (2016-09-05)
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

ctmm 0.3.2 (2016-05-12)
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

ctmm 0.3.1 (2016-02-23)
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
  
ctmm 0.3.0 (2015-11-26)
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

ctmm 0.2.9 (2015-10-13)
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

ctmm 0.2.8 (2015-08-25)
================
* efficiency gains in as.telemetry with multiple animals
* bugfix in plot.telemetry for multiple Gaussian PDFs
* bugfix in variogram for rare condition when fast=TRUE

ctmm 0.2.7 (2015-07-27)
================
* CRAN check compliance achieved.
* all methods (plot, mean, summary, simulate) can and must be run without class extensions
* argument names no longer clash with function names and are more explicit about their object class

ctmm 0.2.6 (2015-07-17)
================
* export bugfixes

ctmm 0.2.5 (2015-07-14)
================
* IOU bug fixes in ctmm.fit and plot.variogram
  
ctmm 0.2.4 (2015-06-28)
================
* cleaned up and labeled tau parameter arrays
* implemented Workaround for when subset demotes S4 objects to S3 objects
* plot.telemetry now enforces asp=1 even with xlim/ylim arguments

ctmm 0.2.3 (2015-06-19)
================
* new function summary.telemetry
* bugfix in as.telemetry for data$t
* bugfix in ctmm.loglike for some cases with numeric underflow
* periodogram and plot.periodogram can now check for spurious periodicities
* minimal support for BM and IOU motion

ctmm 0.2.2 (2015-05-21)
================
* new functions periodogram, plot.periodogram

ctmm 0.2.1 (2015-05-08)
================
* new function SpatialPoints.telemetry returns SpatialPoints object from telemetry data
* new function SpatialPolygonsDataFrame.akde returns akde home-range contour SpatialPolygons objects from akde object
* new function writeShapefile.akde writes akde home-range contours to ESRI shapefile
* new function raster.akde returns akde pdf raster object
* new function summary.akde returns HR area of AKDE
* fixed bad CI in plot.telemetry model option
* as.telemetry now takes a timezone argument for as.POSIXct and defaults to UTC
* telemetry, ctmm, and akde objects now have idenification and projection information slotted, with consistent naming throughout

ctmm 0.2.0 (2015-04-27)
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
