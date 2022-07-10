## ----  fig.show='hold'--------------------------------------------------------
library(ctmm)
data(buffalo)
Pepper <- buffalo$Pepper
M.IID <- ctmm.fit(Pepper) # no autocorrelation timescales
GUESS <- ctmm.guess(Pepper,interactive=FALSE) # automated model guess
M.OUF <- ctmm.fit(Pepper,GUESS) # in general, use ctmm.select instead

## ----  fig.show='hold', results = "hide"--------------------------------------
KDE <- akde(Pepper,M.IID) # KDE
AKDE <- akde(Pepper,M.OUF) # AKDE

## ----  fig.show='hold'--------------------------------------------------------
# plot all sampling intervals
dt.plot(Pepper)

## ----  fig.show='hold', results = "hide"--------------------------------------
# Option 1) exact calculation - slow on larger datasets
# UD2w <- akde(Pepper,M.OUF,weights=TRUE,fast=FALSE,PC='direct')

# Option 2) use tiny dt
dt <- min(diff(Pepper$t))
wAKDE <- akde(Pepper,M.OUF,weights=TRUE,dt=dt)

# Option 3) remove tiny sampling intervals and use smaller dt
# dt <- 1 %#% 'hr'
# ...
# wAKDE <- akde(Pepper,M.OUF,weights=TRUE,dt=dt)

## ----  fig.show='hold', results = "hide", eval=FALSE--------------------------
#  # calculate one extent for all UDs
#  EXT <- extent(list(KDE,AKDE,wAKDE),level=0.95)
#  
#  plot(Pepper,UD=KDE,xlim=EXT$x,ylim=EXT$y)
#  title(expression("IID KDE"["C"]))
#  plot(Pepper,UD=AKDE,xlim=EXT$x,ylim=EXT$y)
#  title(expression("OUF AKDE"["C"]))
#  plot(Pepper,UD=wAKDE,xlim=EXT$x,ylim=EXT$y)
#  title(expression("weighted OUF AKDE"["C"]))

## ----  fig.show='hold', results = "hide", echo=FALSE--------------------------
# calculate one extent for all UDs
EXT <- extent(list(KDE,AKDE,wAKDE),level=0.95)
# sampling intervals in hours
col <- "hr" %#% diff(Pepper$t)
# minimum adjacent sampling interval
col <- pmin(c(Inf,col),c(col,Inf))
# sampling intervals under 1.5 hours
col <- (col < 1.5)
# red (low-frequency) or yellow (high-frequency)
col <- grDevices::rgb(1,col,0)

plot(Pepper,UD=KDE,xlim=EXT$x,ylim=EXT$y,col=col)
title(expression("IID KDE"["C"]))
plot(Pepper,UD=AKDE,xlim=EXT$x,ylim=EXT$y,col=col)
title(expression("OUF AKDE"["C"]))
plot(Pepper,UD=wAKDE,xlim=EXT$x,ylim=EXT$y,col=col)
title(expression("weighted OUF AKDE"["C"]))

## ----  fig.show='hold'--------------------------------------------------------
summary(KDE)
summary(wAKDE)

