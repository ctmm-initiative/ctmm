## ------------------------------------------------------------------------
library(ctmm)
data('coati')
names(coati[[1]]) # imported column names
plot(coati[[1]],col=rainbow(2),error=2,trans=0.4) # coati plot with 95% error discs

## ------------------------------------------------------------------------
data(turtle)
names(turtle[[1]]) # data are not yet calibrated
names(turtle) # two calibration datasets and two turtle datasets
plot(turtle[1:2],col=rainbow(2)) # calibration data only

## ------------------------------------------------------------------------
UERE <- uere.fit(turtle[1:2]) # only using calibration data
summary(UERE)

## ------------------------------------------------------------------------
uere(turtle) <- UERE
names(turtle[[3]]) # now the data are calibrated
plot(turtle[[3]],error=2) # turtle plot with 95% error discs

## ------------------------------------------------------------------------
squirtle <- lapply(turtle,function(t){ t$HDOP <- NULL ; t })

## ------------------------------------------------------------------------
UERE2 <- uere.fit(squirtle[1:2])

## ------------------------------------------------------------------------
summary(list(HDOP=UERE,homo=UERE2))

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
# automated guestimates with circular covariance and calibrated errors
GUESS <- ctmm.guess(turtle[[3]],CTMM=ctmm(error=TRUE),interactive=FALSE)
# the beta optimizer is more reliable than the default optimizer
# control <- list(method='pNewton',cores=-1) # will use all but one cores
control <- list(method='pNewton',cores=2) # CRAN policy limits us to 2 cores
# stepwise fitting
FIT <- ctmm.select(turtle[[3]],GUESS,control=control,trace=TRUE)
summary(FIT)

## ------------------------------------------------------------------------
# delete UERE information
uere(turtle[[3]]) <- NULL

## ------------------------------------------------------------------------
# cheat and use previous fit as initial guess
GUESS$error <- 10 # 10 meter error guess
# fit parameter estimates
FIT <- ctmm.select(turtle[[3]],GUESS,control=control)
summary(FIT)

