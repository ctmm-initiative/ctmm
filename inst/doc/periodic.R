## ----  fig.show="hold"---------------------------------------------------
library(ctmm)
data(wolf)
Tay <- wolf$Tay
plot(Tay)

## ------------------------------------------------------------------------
LSP <- periodogram(Tay,fast=2,res.time=2)

## ----  fig.show="hold"---------------------------------------------------
plot(LSP,max=TRUE,diagnostic=TRUE,transparency=0.2,cex=0.5)

## ------------------------------------------------------------------------
PROTO <- ctmm(mean="periodic",period=c(1 %#% "month",1 %#% "week",1 %#% "day"),circle=TRUE)

