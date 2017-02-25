## ----  fig.show='hold'---------------------------------------------------
library(ctmm)
data("buffalo")
cilla <- buffalo[[1]]
plot(cilla)
title("1 Buffalo")
plot(buffalo,col=rainbow(length(buffalo)))
title("5 Buffalo")

## ----  fig.show='hold'---------------------------------------------------
SVF <- variogram(cilla)
plot(SVF,fraction=0.005)
title("zoomed in")
plot(SVF,fraction=0.65,level=c(0.5,0.95))
title("zoomed out")

## ----  fig.show='hold'---------------------------------------------------
m0 <- ctmm(sigma=23*1000^2) # 23 km^2 in m^2
m1 <- ctmm(sigma=23*1000^2,tau=6*24*60^2) # and 6 days in seconds
plot(SVF,CTMM=m0,fraction=0.65,level=c(0.5,0.95),col.CTMM="red")
title("m0: IID")
plot(SVF,CTMM=m1,fraction=0.65,level=c(0.5,0.95),col.CTMM="purple")
title("m1: OU")

## ----  fig.show='hold'---------------------------------------------------
m2 <- ctmm(sigma=23*1000^2,tau=c(6*24*60^2,1*60^2)) # and 1 hour in seconds
plot(SVF,CTMM=m1,fraction=0.002,col.CTMM="purple")
title("m1: OU")
plot(SVF,CTMM=m2,fraction=0.002,col.CTMM="blue")
title("m2: OUF")

## ----  fig.show='hold'---------------------------------------------------
plot(SVF,CTMM=m1,fraction=0.65,level=c(0.5,0.95),col.CTMM="purple")
title("m1: OU")
plot(SVF,CTMM=m2,fraction=0.65,level=c(0.5,0.95),col.CTMM="blue")
title("m2: OUF")

## ----  fig.show='hold'---------------------------------------------------
# simulate fake buffalo with the same sampling schedule
willa <- simulate(m2,t=cilla$t)
plot(willa)
title("simulation")
# now calculate and plot its variogram
SVF2 <- variogram(willa)
plot(SVF2,CTMM=m2,fraction=0.65,level=c(0.5,0.95),col.CTMM="blue")
title("simulation")

## ---- fig.show='hold'----------------------------------------------------
data("gazelle")
SVF3 <- variogram(gazelle[[18]])
plot(SVF3,fraction=0.85,level=c(0.5,0.95))
title("Default method")
# 1, 5, 25 hour sampling intervals
dt <- 60*60*c(1,5,25)
SVF3 <- variogram(gazelle[[18]],dt=dt)
plot(SVF3,fraction=0.85,level=c(0.5,0.95))
title("Multi method")

## ---- fig.show='hold'----------------------------------------------------
# 1 hour sampling intervals
dt = 60*60
# buffalo 4 is bad
SVF4 <- lapply(buffalo[-4],function(b){ variogram(b,dt=dt) })
SVF4 <- mean(SVF4)
plot(SVF4,fraction=0.35,level=c(0.5,0.95))
title("Population variogram")

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
DATA <- as.telemetry(system.file("extdata","leroy.csv.gz",package="move"))
# 1 hour and 1 day autocorrelation timescales
GUESS <- ctmm(tau=c(1/4,24)*60^2)
# first fit without telemetry error
FITS <- list()
FITS$NOERR <- ctmm.fit(DATA,GUESS)
# second fit based on first with telemetry error
GUESS <- FITS$NOERR
GUESS$error <- TRUE
FITS$ERROR <- ctmm.fit(DATA,GUESS)
# model improvement
summary(FITS)

## ----  fig.show='hold'---------------------------------------------------
M0 <- ctmm.fit(cilla,m0)
summary(M0)

## ----  fig.show='hold'---------------------------------------------------
M1 <- ctmm.fit(cilla,m1)
summary(M1)

## ----  fig.show='hold'---------------------------------------------------
M2 <- ctmm.fit(cilla,m2)
summary(M2)

## ----  fig.show='hold'---------------------------------------------------
FITS <- list(IID=M0,OU=M1,OUF=M2)
summary(FITS)

## ------------------------------------------------------------------------
FITZ <- ctmm.select(cilla,m2,verbose=TRUE,level=1)
summary(FITZ)

## ----  fig.show='hold'---------------------------------------------------
plot(SVF,CTMM=FITS,col.CTMM=c("red","purple","blue"),fraction=0.65,level=0.5)
title("zoomed out")
plot(SVF,CTMM=FITS,col.CTMM=c("red","purple","blue"),fraction=0.002,level=0.5)
title("zoomed in")

