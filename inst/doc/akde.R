## ----  fig.show='hold'---------------------------------------------------
library(ctmm)
data(buffalo)
pepper <- buffalo$Pepper
M0 <- ctmm.fit(pepper) # no autocorrelation timescales
m2 <- ctmm(tau=c(6*24,1)*60^2) # ~ 6 day and 1 hour autocorrelation timescales
M2 <- ctmm.fit(pepper,m2)

## ----  fig.show='hold', results = "hide"---------------------------------
UD0 <- akde(pepper,M0)
UD2 <- akde(pepper,M2)
UD2w <- akde(pepper,M2,weights=TRUE)
# calculate one extent for all UDs
EXT <- extent(list(UD0,UD2,UD2w),level=0.95)

## ----  fig.show='hold'---------------------------------------------------
plot(pepper,UD=UD0,xlim=EXT$x,ylim=EXT$y)
title(expression("IID KDE"["C"]))
plot(pepper,UD=UD2,xlim=EXT$x,ylim=EXT$y)
title(expression("OUF AKDE"["C"]))
plot(pepper,UD=UD2w,xlim=EXT$x,ylim=EXT$y)
title(expression("weighted OUF AKDE"["C"]))

## ----  fig.show='hold'---------------------------------------------------
summary(UD0)
summary(UD2w)

