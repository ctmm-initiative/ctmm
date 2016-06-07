## ----  fig.show='hold'---------------------------------------------------
library(ctmm)
data(buffalo)
cilla <- buffalo[[1]]
M0 <- ctmm.fit(cilla) # no autocorrelation timescales
m2 <- ctmm(tau=c(6*24,1)*60^2) # ~ 6 day and 1 hour autocorrelation timescales
M2 <- ctmm.fit(cilla,m2) 

## ----  fig.show='hold', results = "hide"---------------------------------
UD0 <- akde(cilla,M0)
UD2 <- akde(cilla,M2)

## ----  fig.show='hold'---------------------------------------------------
plot(cilla,UD=UD0,ylim=c(-14,12)*1000)
title("M0")
plot(cilla,UD=UD2,ylim=c(-14,12)*1000)
title("M2")

## ----  fig.show='hold'---------------------------------------------------
summary(UD0)
summary(UD2)

