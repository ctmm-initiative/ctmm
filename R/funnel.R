funnel <- function(data,object,estimate="area",precision="t",level=0.95,level.UD=0.95,...)
{
  precision <- match.arg(precision,c("n","N","t"))
  if(class(object)[1]=="speed") { estimate <- "speed" }

  if(estimate=="tau position") { estimate <- "position" }
  else if(estimate=="tau velocity") { estimate <- "velocity" }

  alpha <- 1-level

  for(i in 1:length(object))
  {
    # extract estimate information
    CI <- summary(object[[i]],level=level,level.UD=level.UD,units=FALSE)
    DOF <- CI$DOF
    CI <- CI$CI
    NAMES <- rownames(CI)
    ROW <- grepl(estimate,NAMES)
    CI <- CI[ROW,]
    object[[i]] <- CI

    if(precision=="n")
    { data[[i]] <- nrow(data[[i]]) }
    else if(precision=="N")
    {
      if(estimate=="area")
      { data[[i]] <- DOF['area'] }
      else if(estimate %in% c("speed","velocity"))
      { data[[i]] <- DOF['speed'] }
      else if(estimate %in% c("diffusion","position"))
      { data[[i]] <- DOF['diffusion'] }
    }
    else if(precision=="t")
    {
      if(estimate %in% c("area","position","diffusion"))
      { data[[i]] <- last(data[[i]]$t) - first(data[[i]]$t) }
      else if(estimate %in% c("speed","velocity"))
      { data[[i]] <- stats::quantile(diff(data[[i]]$t),probs=c(alpha/2,0.5,1-alpha/2)) }
    }
  }

  # some individuals may be missing model features
  GOOD <- sapply(object,length) > 0
  object <- object[GOOD]
  data <- data[GOOD]

  data <- sapply(data,identity)
  object <- sapply(object,identity) # [est,ind]

  if(estimate=="area")
  {
    xlab <- paste0("Area")
    UNITS <- unit(object,"area")
  }
  else if(estimate=="speed")
  {
    xlab <- "Speed"
    UNITS <- unit(object,"speed")
  }
  else if(estimate=="diffusion")
  {
    xlab <- "Diffusion Rate"
    UNITS <- unit(object,"diffusion")
  }
  else if(estimate=="position")
  {
    xlab <- "\u03C4[position]"
    UNITS <- unit(object,"time")
  }
  else if(estimate=="velocity")
  {
    xlab <- "\u03C4[velocity]"
    UNITS <- unit(object,"time")
  }
  xlab <- paste0(xlab," (",UNITS$name,")")
  object <- object/UNITS$scale

  if(precision=="n")
  { ylab <- "Sample Size" }
  else if(precision=="N")
  { ylab <- "Effective Sample Size" }
  else if(precision=="t")
  {
    if(estimate %in% c("area","position","diffusion"))
    {
      ylab <- "Sampling Period"
      UNITS <- unit(data,"time")
    }
    else if(estimate %in% c("speed","velocity"))
    {
      precision <- "f"
      ylab <- "Sampling Frequency"
      data <- 1/data[3:1,]
      UNITS <- unit(data,"frequency")
    }
    data <- data/UNITS$scale
    ylab <- paste0(ylab," (",UNITS$name,")")
  }

  x <- object
  y <- data

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

  x <- object[2,]
  y <- data
  if(length(dim(y)))
  {
    SUB <- y[3,]-y[1,] > .Machine$double.eps # still does not avoid annoying warning
    suppressWarnings( graphics::arrows(x[SUB],y[1,SUB],x[SUB],y[3,SUB],length=0.05,angle=90,code=3,...) )
  }
  # will need to switch from arrows to segment to just avoid annoying warning...
}
