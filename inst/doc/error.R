## ------------------------------------------------------------------------
library(ctmm)
data(turtle)
names(turtle[[1]]) # data are not yet calibrated, but HDOP and location class is present
names(turtle) # two calibration datasets and two turtle datasets
plot(turtle[1:2],col=rainbow(2)) # calibration data only

## ------------------------------------------------------------------------
UERE <- uere.fit(turtle[1:2]) # only using calibration data
summary(UERE)

## ------------------------------------------------------------------------
uere(turtle) <- UERE
names(turtle[[3]]) # now the data are calibrated, as VAR is present
plot(turtle[[3]],error=2) # turtle plot with 95% error discs

## ------------------------------------------------------------------------
data(coati)
names(coati[[1]]) # VAR already present
plot(coati[[1]],error=2) # coati plot with 95% error discs

## ------------------------------------------------------------------------
t.noHDOP  <- lapply(turtle,function(t){ t$HDOP  <- NULL; t })
t.noclass <- lapply(turtle,function(t){ t$class <- NULL; t })
t.nothing <- lapply(turtle,function(t){ t$HDOP  <- NULL; t$class <- NULL; t })

## ------------------------------------------------------------------------
UERE.noHDOP  <- uere.fit(t.noHDOP[1:2])
UERE.noclass <- uere.fit(t.noclass[1:2])
UERE.nothing <- uere.fit(t.nothing[1:2])

## ------------------------------------------------------------------------
summary(list(HDOP.class=UERE,class=UERE.noHDOP,HDOP=UERE.noclass,homoskedastic=UERE.nothing))

## ------------------------------------------------------------------------
UERES <- lapply(turtle[1:2],uere.fit)
summary(list(joint=UERE,individual=UERES))

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
control <- list(method='pNewton')
# stepwise fitting # CRAN policy limits us to 2 cores
FIT <- ctmm.select(turtle[[3]],GUESS,control=control,trace=TRUE,cores=2)
summary(FIT)

## ------------------------------------------------------------------------
# delete UERE information
uere(turtle) <- NULL
# toss out 2D locations for now
turtle <- lapply(turtle,function(t){ t[t$class=="3D",] })

## ------------------------------------------------------------------------
# cheat and use previous fit as initial guess
GUESS$error <- 10 # 10 meter error guess
# fit parameter estimates
FIT <- ctmm.select(turtle[[3]],GUESS,control=control,trace=TRUE,cores=2)
summary(FIT)

