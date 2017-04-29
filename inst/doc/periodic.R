## ----  fig.show="hold"---------------------------------------------------
library(ctmm)
data(wolf)
Tay <- wolf$Tay
plot(Tay)

## ------------------------------------------------------------------------
LSP <- periodogram(Tay,fast=2,res.time=2)

## ----  fig.show="hold"---------------------------------------------------
plot(LSP, max=TRUE, diagnostic=TRUE)

