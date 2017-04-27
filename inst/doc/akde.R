## ----  fig.show='hold'---------------------------------------------------
library(ctmm)
data(buffalo)
Pepper <- buffalo$Pepper
M.IID <- ctmm.fit(Pepper) # no autocorrelation timescales
m.ouf <- ctmm.guess(Pepper,interactive=FALSE) # automated model guess
M.OUF <- ctmm.fit(Pepper,m.ouf)

## ----  fig.show='hold', results = "hide"---------------------------------
UD0 <- akde(Pepper,M.IID)
UD2 <- akde(Pepper,M.OUF)
UD2w <- akde(Pepper,M.OUF,weights=TRUE)
# calculate one extent for all UDs
EXT <- extent(list(UD0,UD2,UD2w),level=0.95)

## ----  fig.show='hold'---------------------------------------------------
plot(Pepper,UD=UD0,xlim=EXT$x,ylim=EXT$y)
title(expression("IID KDE"["C"]))
plot(Pepper,UD=UD2,xlim=EXT$x,ylim=EXT$y)
title(expression("OUF AKDE"["C"]))
plot(Pepper,UD=UD2w,xlim=EXT$x,ylim=EXT$y)
title(expression("weighted OUF AKDE"["C"]))

## ----  fig.show='hold'---------------------------------------------------
summary(UD0)
summary(UD2w)

