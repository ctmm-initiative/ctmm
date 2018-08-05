## ----  fig.show="hold"---------------------------------------------------
library(ctmm)
data(wolf)
Gamba <- wolf$Gamba
plot(Gamba)

## ------------------------------------------------------------------------
LSP <- periodogram(Gamba,fast=2,res.time=2)

## ----  fig.show="hold"---------------------------------------------------
plot(LSP,max=TRUE,diagnostic=TRUE,cex=0.5)

## ------------------------------------------------------------------------
PROTO <- ctmm(mean="periodic",period=c(1 %#% "day",1 %#% "month"),circle=TRUE)

## ------------------------------------------------------------------------
SVF <- variogram(Gamba,res=3)
GUESS <- ctmm.guess(Gamba,PROTO,variogram=SVF,interactive=FALSE)

## ------------------------------------------------------------------------
# ctmm beta optimizer is more reliable here
# control <- list(method="pNewton",cores=-1) # use all but 1 core
control <- list(method="pNewton",cores=2) # CRAN policy limits to 2 processes
FITS <- ctmm.select(Gamba,GUESS,verbose=TRUE,control=control)

## ------------------------------------------------------------------------
summary(FITS,units=FALSE)

## ------------------------------------------------------------------------
# sorted by MSPE
summary(FITS[1:7])
# sorted by IC
summary(FITS[c(5,8:10)])

## ------------------------------------------------------------------------
summary(FITS,MSPE=FALSE)

## ------------------------------------------------------------------------
"hour" %#% stats::median(diff(Gamba$t))

## ------------------------------------------------------------------------
summary(FITS[[7]]) # harmonic 3 0 # selected by IC
summary(FITS[[1]]) # harmonic 2 0 # selected by MSPE

## ---- fig.show='hold'----------------------------------------------------
xlim <- c(0,1/2) %#% "month"
plot(SVF,CTMM=FITS[[7]],xlim=xlim)
title("3 Harmonics")
plot(SVF,CTMM=FITS[[1]],xlim=xlim)
title("2 Harmonics")

