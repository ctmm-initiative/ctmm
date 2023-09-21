## -----------------------------------------------------------------------------
library(ctmm)
data(turtle)
names(turtle[[1]]) # data are not yet calibrated, but HDOP and location class is present
names(turtle) # two calibration datasets and two turtle datasets
plot(turtle[1:2],col=rainbow(2)) # calibration data only

## -----------------------------------------------------------------------------
UERE <- uere.fit(turtle[1:2]) # only using calibration data
summary(UERE)

## -----------------------------------------------------------------------------
uere(turtle) <- UERE
summary(uere(turtle[[3]])) # this should be the same as summary(UERE)
plot(turtle[[3]],error=2) # turtle plot with 95% error discs

## -----------------------------------------------------------------------------
data(pelican)
names(pelican)
names(pelican$argos) # error ellipse information (COV and VAR) already present
plot(pelican$argos) # pelican Argos plot with 95% error ellipses

## -----------------------------------------------------------------------------
t.noHDOP  <- lapply(turtle,function(t){ t$HDOP  <- NULL; t })
t.noclass <- lapply(turtle,function(t){ t$class <- NULL; t })
t.nothing <- lapply(turtle,function(t){ t$HDOP  <- NULL; t$class <- NULL; t })

## -----------------------------------------------------------------------------
UERE.noHDOP  <- uere.fit(t.noHDOP[1:2])
UERE.noclass <- uere.fit(t.noclass[1:2])
UERE.nothing <- uere.fit(t.nothing[1:2])

## -----------------------------------------------------------------------------
summary(list(HDOP.class=UERE,class=UERE.noHDOP,HDOP=UERE.noclass,homoskedastic=UERE.nothing))

## -----------------------------------------------------------------------------
UERES <- lapply(turtle[1:2],uere.fit)
summary(list(joint=UERE,individual=UERES))

## -----------------------------------------------------------------------------
outlie(turtle[[3]]) -> OUT

## -----------------------------------------------------------------------------
plot(OUT,units=FALSE)

## -----------------------------------------------------------------------------
BAD <- OUT$speed>0.08 # not appropriate for other species!
turtle[[3]] <- turtle[[3]][!BAD,]
outlie(turtle[[3]]) -> OUT

## ----fig.show='hold', echo=FALSE----------------------------------------------
# Argos type errors
curve(1+x,0,5,xlab="Short time lag",ylab="Semi-variance",ylim=c(0,6))
points(c(0,0),c(0,1))
title("Argos")
# detector array type errors (qualitatively only)
curve((1-exp(-2*x))/(1-exp(-2/4)),0,1/4,xlab="Short time lag",ylab="Semi-variance",ylim=c(0,6),xlim=c(0,5),add=FALSE)
curve(3/4+x,1/4,5,xlab="Short time lag",ylab="Semi-variance",ylim=c(0,6),add=TRUE,xlim=c(0,5))
points(1/4,1)
title("Detector Array")

## -----------------------------------------------------------------------------
# automated guesstimate for calibrated data
GUESS <- ctmm.guess(turtle[[3]],CTMM=ctmm(error=TRUE),interactive=FALSE)
# stepwise fitting # CRAN policy limits us to 2 cores
FIT <- ctmm.select(turtle[[3]],GUESS,trace=TRUE,cores=2)
# if you get errors on your platform, then try cores=1
summary(FIT)

## -----------------------------------------------------------------------------
# delete UERE information
uere(turtle) <- NULL

## -----------------------------------------------------------------------------
# automated guesstimate for uncalibrated data (with 10 meter RMS UERE guess)
GUESS <- ctmm.guess(turtle[[3]],CTMM=ctmm(error=10),interactive=FALSE)
# fit and select models # CRAN policy limits us to 2 cores
FIT <- ctmm.select(turtle[[3]],GUESS,trace=TRUE,cores=2)
# if you get errors on your platform, then try cores=1
summary(FIT)

## -----------------------------------------------------------------------------
# these data have two location classes: 2D & 3D
summary(uere(turtle))
# assign 20-meter 2D RMS UERE and 10-meter 3D RMS UERE
uere(turtle) <- c(20,10)
# for one location class, the above would likely be unnecessary, but would look like
# uere(turtle) <- 10

# the default uncertainty is none for numerical assignments
UERE <- uere(turtle)
summary(UERE)
# this is because the degrees-of-freedom are set to Inf
UERE$DOF
# here I set the DOF to a smaller value
UERE$DOF[] <- 2
# which now gives plausible credible intervals
summary(UERE)
# assign the prior to the data
uere(turtle) <- UERE
# automated guesstimate for calibrated data
GUESS <- ctmm.guess(turtle[[3]],CTMM=ctmm(error=TRUE),interactive=FALSE)
# stepwise fitting # CRAN policy limits us to 2 cores
FIT <- ctmm.select(turtle[[3]],GUESS,trace=TRUE,cores=2)
# if you get weird errors on your platform, then try cores=1
summary(FIT)

