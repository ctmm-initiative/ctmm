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
# circle=TRUE takes longer to fit and gives optim more problems
# PROTO <- ctmm(mean="periodic",period=c(1 %#% "day",1 %#% "month"),circle=TRUE)
PROTO <- ctmm(mean="periodic",period=c(1 %#% "day",1 %#% "month"))

## ------------------------------------------------------------------------
SVF <- variogram(Gamba,res=3)
GUESS <- ctmm.guess(Gamba,PROTO,variogram=SVF,interactive=FALSE)

## ------------------------------------------------------------------------
control <- list(method="pNewton") # beta optimizer is much faster here
control <- list(method="pNewton",cores=2) # CRAN policy limits to 2 processes
FITS <- ctmm.select(Gamba,GUESS,verbose=TRUE,control=control)

## ------------------------------------------------------------------------
summary(FITS)

## ------------------------------------------------------------------------
"hour" %#% stats::median(diff(Gamba$t))

## ------------------------------------------------------------------------
summary(FITS[[1]]) # harmonic 3 0
summary(FITS[[4]]) # harmonic 2 0

## ---- fig.show='hold'----------------------------------------------------
xlim <- c(0,1/2) %#% "month"
plot(SVF,CTMM=FITS[[1]],xlim=xlim)
title("3 Harmonics")
plot(SVF,CTMM=FITS[[4]],xlim=xlim)
title("2 Harmonics")

## ------------------------------------------------------------------------
summary(FITS,IC="BIC")

## ------------------------------------------------------------------------
summary(FITS,IC="MSPE")

