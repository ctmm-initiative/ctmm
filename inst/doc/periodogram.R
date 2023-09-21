## ----fig.show="hold"----------------------------------------------------------
library(ctmm)
data(wolf)
Gamba <- wolf$Gamba
plot(Gamba)

## -----------------------------------------------------------------------------
LSP <- periodogram(Gamba,fast=2,res.time=2)

## ----fig.show="hold"----------------------------------------------------------
plot(LSP,max=TRUE,diagnostic=TRUE,cex=0.5)

## -----------------------------------------------------------------------------
PROTO <- ctmm(mean="periodic",period=c(24 %#% "hours",1 %#% "month"),circle=TRUE)

## -----------------------------------------------------------------------------
SVF <- variogram(Gamba,res=3)
GUESS <- ctmm.guess(Gamba,PROTO,variogram=SVF,interactive=FALSE)

## -----------------------------------------------------------------------------
# CRAN policy limits to 2 processes (cores)
FITS <- ctmm.select(Gamba,GUESS,verbose=TRUE,cores=2)
# if you get errors on your platform, then try cores=1

## -----------------------------------------------------------------------------
summary(FITS)

## -----------------------------------------------------------------------------
# these are sorted by MSPE
SUB <- grepl("OUF anisotropic harmonic",names(FITS))
summary(FITS[SUB])
# these are sorted by IC
SUB <- grepl('0 0',names(FITS))
summary(FITS[SUB])

## -----------------------------------------------------------------------------
# sorting by IC only
summary(FITS,MSPE=NA)

## -----------------------------------------------------------------------------
"hour" %#% stats::median(diff(Gamba$t))

## -----------------------------------------------------------------------------
SUB <- rownames(summary(FITS,MSPE=NA))[1]
summary(FITS[[SUB]]) # harmonic 3 0 # selected by IC
summary(FITS[[1]])   # harmonic 2 0 # selected by IC/MSPE

## ----fig.show='hold'----------------------------------------------------------
xlim <- c(0,1/2) %#% "month"
plot(SVF,CTMM=FITS[[SUB]],xlim=xlim)
title("3 Harmonics")
plot(SVF,CTMM=FITS[[1]],xlim=xlim)
title("2 Harmonics")

