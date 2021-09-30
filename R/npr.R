revisitation <- function(data,UD,error=0.001,...)
{ npr(data,UD,variable="revisitation",error=error,...) }


# TODO Occurrence NPRs ##
# TODO NPR debias
npr <- function(data,UD,variable="revisitation",error=0.001,...)
{
  # arbitrary variables allowed
  # variable <- match.arg(variable,c('revisitation','speed'))

  info <- attr(data,"info")
  type <- UD@type

  CTMM <- UD@CTMM
  axes <- CTMM$axes
  res <- 10

  if(length(CTMM$tau)<2)
  {
    if(variable=='revisitation')
    {
      warning("Revisitation estimation requires a continuous-velocity model.")
      return(UD)
    }
    else if(variable=='speed')
    { stop("Speed estimation requires a continuous-velocity model.") }
  }

  if(type=="occurrence")
  { return(occurrence(data,CTMM,variable=variable,error=error,grid=UD,...)) }

  data <- predict(CTMM,data=data,t=data$t)
  if(variable %in% c("revisitation","speed")) # append a debiased speed column
  { data <- speeds.fast(data,append=TRUE) }
  CTMM$error <- FALSE # smoothed error model (approximate)

  GRID <- kde.grid(data,H=UD$H,axes=axes,alpha=error,grid=UD)

  weights <- UD$weights
  UD <- as.list(UD)
  UD$PDF <- UD$CDF <- UD$MISE <- NULL

  if(variable=="revisitation")
  {
    weights <- weights * data$speed
    weight <- sum(weights)

    UD$axes <- UD$r <- UD$dr <- NULL

    UD <- c(UD,kde(data=data,H=UD$H,W=weights,alpha=error,grid=GRID,bias=UD$bias))
    UD$weight <- weight # needed for averaging over individuals
    UD <- new.UD(UD,info=info,type='range',variable=variable,CTMM=CTMM)
  }
  else
  {
    stop()

    UD$RS <- kde(data=data,H=UD$H,W=weights,alpha=error,grid=GRID,variable=variable)

    NAMES <- names(UD) # why is this being erased?
    #UD <- new.RS(UD,info=attr(data,"info"),type='range',variable=variable,CTMM=CTMM)
    names(UD) <- NAMES # why is this necessary?
  }

  return(UD)
}
