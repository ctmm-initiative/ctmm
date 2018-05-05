## ------------------------------------------------------------------------
library(ctmm)
fisher <- as.telemetry(system.file("extdata","leroy.csv.gz",package="move"))
names(fisher) # imported column names
plot(fisher,error=2) # fisher plot with 95% error discs

## ------------------------------------------------------------------------
library(ctmm)
data(turtle)
names(turtle)
plot(turtle[1:2],col=rainbow(2)) # calibration data only

## ------------------------------------------------------------------------
UERE <- uere(turtle[1:2]) # only using calibration data
UERE

## ------------------------------------------------------------------------
uere(turtle) <- UERE
plot(turtle[[3]],error=2) # turtle plot with 95% error discs

## ------------------------------------------------------------------------
outlie(turtle[[3]]) -> OUT

## ------------------------------------------------------------------------
plot(OUT)

## ------------------------------------------------------------------------
BAD <- which.max(OUT$speed)
turtle[[3]] <- turtle[[3]][-BAD,]
outlie(turtle[[3]]) -> OUT

## ----  fig.show='hold', echo=FALSE---------------------------------------
# ARGOS type errors
curve(1+x,0,5,xlab="Short time lag",ylab="Semi-variance",ylim=c(0,6))
points(c(0,0),c(0,1))
title("ARGOS")
# detector array type errors (qualitatively only)
curve((1-exp(-2*x))/(1-exp(-2/4)),0,1/4,xlab="Short time lag",ylab="Semi-variance",ylim=c(0,6),xlim=c(0,5),add=FALSE)
curve(3/4+x,1/4,5,xlab="Short time lag",ylab="Semi-variance",ylim=c(0,6),add=TRUE,xlim=c(0,5))
points(1/4,1)
title("Detector Array")

## ------------------------------------------------------------------------
# ctmm beta optimizer is more reliable here
# control <- list(method='pNewton',cores=-1) # use all but 1 core
control <- list(method='pNewton',cores=2) # CRAN policy limits us to 2 processes

# default model guess without error
GUESS <- ctmm.guess(fisher,interactive=FALSE)
# first fit without telemetry error
FIT <- ctmm.fit(fisher,GUESS,control=control)
summary(FIT)

# second fit based on first, but with telemetry error activated
GUESS <- FIT
GUESS$error <- TRUE
FIT <- ctmm.fit(fisher,GUESS,control=control)
summary(FIT)

## ------------------------------------------------------------------------
# automated guestimates with circular covariance and calibrated errors
GUESS <- ctmm.guess(turtle[[3]],CTMM=ctmm(error=TRUE,isotropic=TRUE),interactive=FALSE)
# first fit circular model
FIT <- ctmm.fit(turtle[[3]],GUESS,control=control)
summary(FIT)

# create elliptical model guess
GUESS <- FIT
GUESS$isotropic <- FALSE
# populate eccentricity and orientation with guesstimates
GUESS <- ctmm.guess(turtle[[3]],CTMM=GUESS,interactive=FALSE)
# second fit elliptical model
FIT <- ctmm.fit(turtle[[3]],GUESS,control=control)
# velocity does not appear supported
summary(FIT)

# proceed with model selection
FIT <- ctmm.select(turtle[[3]],FIT,control=control)
# velocity was not supported by this data
summary(FIT)

## ------------------------------------------------------------------------
# delete UERE information
uere(turtle[[3]]) <- NULL

## ------------------------------------------------------------------------
# cheat and use previous fit as initial guess
GUESS <- FIT
GUESS$error <- 10 # 10 meter error guess
# fit parameter estimates
FIT <- ctmm.fit(turtle[[3]],GUESS,control=control)
summary(FIT)

