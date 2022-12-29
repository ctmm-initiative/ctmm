funnel <- function(x,y,variable="area",precision="t",level=0.95,level.UD=0.95,...)
{
  variable <- canonical.name(variable)
  variable <- match.arg(variable,c("area","diffusion","speed","tauposition","tauvelocity","distance"))
  effect <- variable

  if(precision=="N") { precision <- "DOF" }
  precision <- canonical.name(precision)
  precision <- match.arg(precision,c("n","dof","t"))

  # sort x and y
  if(class(x)[1]=="telemetry")
  {
    z <- x
    x <- y
    y <- z
    rm(z)
  }

  effect.dat <- x
  precision.dat <- y

  R <- import.variable(x,variable=variable)
  ID <- R$ID
  VAR <- R$AREA
  DOF <- R$DOF

  alpha <- 1-level

  for(i in 1:length(effect.dat))
  {
    # extract effect information
    if(class(effect.dat[[i]])[1]=="list")
    {
      # could be in non-SI units
      CI <- effect.dat[[i]]$CI
      CI <- t( (VAR[i]/CI[,2]) * t(CI) )
    }
    else
    { CI <- summary(effect.dat[[i]],level=level,level.UD=level.UD,units=FALSE)$CI }

    NAMES <- rownames(CI)
    ROW <- grepl(effect,NAMES)
    CI <- CI[ROW,]
    effect.dat[[i]] <- CI

    if(precision=="n")
    { precision.dat[[i]] <- nrow(precision.dat[[i]]) }
    else if(precision=="dof")
    { precision.dat[[i]] <- DOF[i] }
    else if(precision=="t")
    {
      if(effect %in% c("area","tauposition","diffusion"))
      { precision.dat[[i]] <- last(precision.dat[[i]]$t) - first(precision.dat[[i]]$t) }
      else if(effect %in% c("speed","tauvelocity"))
      { precision.dat[[i]] <- stats::quantile(diff(precision.dat[[i]]$t),probs=c(alpha/2,0.5,1-alpha/2)) }
    }
  }

  # some individuals may be missing model features
  GOOD <- sapply(effect.dat,length) > 0
  effect.dat <- effect.dat[GOOD]
  precision.dat <- precision.dat[GOOD]

  precision.dat <- sapply(precision.dat,identity)
  effect.dat <- sapply(effect.dat,identity) # [est,ind]

  if(effect=="area")
  {
    xlab <- paste0("Area")
    UNITS <- unit(effect.dat,"area")
  }
  else if(effect=="speed")
  {
    xlab <- "Speed"
    UNITS <- unit(effect.dat,"speed")
  }
  else if(effect=="diffusion")
  {
    xlab <- "Diffusion Rate"
    UNITS <- unit(effect.dat,"diffusion")
  }
  else if(effect=="tauposition")
  {
    xlab <- "\u03C4[position]"
    UNITS <- unit(effect.dat,"time")
  }
  else if(effect=="tauvelocity")
  {
    xlab <- "\u03C4[velocity]"
    UNITS <- unit(effect.dat,"time")
  }
  xlab <- paste0(xlab," (",UNITS$name,")")
  effect.dat <- effect.dat/UNITS$scale

  if(precision=="n")
  { ylab <- "Sample Size" }
  else if(precision=="dof")
  { ylab <- "Effective Sample Size" }
  else if(precision=="t")
  {
    if(effect %in% c("area","tauposition","diffusion"))
    {
      ylab <- "Sampling Period"
      UNITS <- unit(precision.dat,"time")
    }
    else if(effect %in% c("speed","tauvelocity"))
    {
      precision <- "f"
      ylab <- "Sampling Frequency"
      precision.dat <- 1/precision.dat[3:1,]
      UNITS <- unit(precision.dat,"frequency")
    }
    precision.dat <- precision.dat/UNITS$scale
    ylab <- paste0(ylab," (",UNITS$name,")")
  }

  x <- effect.dat
  y <- precision.dat

  xlim <- range(x[x>0 | x<Inf])
  ylim <- range(y)

  # base plot
  if(length(dim(y))) { y <- y[2,] }
  plot(x[2,],y,xlim=xlim,ylim=ylim,pch=19,xlab=xlab,ylab=ylab,log='x',...)

  # ERROR BAR PLOT
  # hack: we draw arrows but with very special "arrowheads"
  # (c) Laryx Decidua
  # very annoying warning when zero length --- cannot be suppressed with suppressWarnings()
  if(length(dim(x)))
  {
    SUB <- x[3,]-x[1,] > .Machine$double.eps # still does not avoid annoying warning
    suppressWarnings( graphics::arrows(x[1,SUB],y[SUB],x[3,SUB],y[SUB],length=0.05,angle=90,code=3,...) )
  }

  x <- effect.dat[2,]
  y <- precision.dat
  if(length(dim(y)))
  {
    SUB <- y[3,]-y[1,] > .Machine$double.eps # still does not avoid annoying warning
    suppressWarnings( graphics::arrows(x[SUB],y[1,SUB],x[SUB],y[3,SUB],length=0.05,angle=90,code=3,...) )
  }
  # will need to switch from arrows to segment to just avoid annoying warning...
}
