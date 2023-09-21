revisitation <- function(data,UD,debias=TRUE,error=0.001,...)
{ npr(data,UD,variable="revisitation",normalize=TRUE,debias=debias,error=error,...) }


# TODO Occurrence NPRs ##
# TODO NPR debias
npr <- function(data,UD,variable="speed",normalize=FALSE,debias=TRUE,error=0.001,...)
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

  VAR <- data[[variable]]
  data <- predict(CTMM,data=data,t=data$t)

  if(variable %in% c("revisitation","speed")) # append a debiased speed column
  { data <- speeds_fast(data,append=TRUE) }
  else
  {
    data[[variable]] <- VAR
    rm(VAR)

    if(normalize && any(data[[variable]]<0))
    { stop("Negative variables and cannot provide a weighted distribution.") }
  }

  CTMM$error <- FALSE # smoothed error model (approximate)

  GRID <- kde.grid(data,H=UD$H,axes=axes,alpha=error,grid=UD)

  weights <- UD$weights
  UD <- as.list(UD)

  if(debias)
  { bias <- UD$bias } # consider recalculating for normalize=TRUE
  else
  { bias <- FALSE }

  if(normalize) # generate a new disribution
  {
    if(variable=="revisitation")
    { weights <- weights * data$speed }
    else
    { weights <- weights * data[[variable]] }
    # total weight for means
    weight <- sum(weights)

    UD$PDF <- UD$CDF <- UD$MISE <- NULL
    UD$axes <- UD$r <- UD$dr <- NULL

    UD <- c(UD,kde(data=data,H=UD$H,W=weights,alpha=error,grid=GRID,bias=bias))

    UD$weight <- weight # needed for averaging over individuals
    UD <- new.UD(UD,info=info,type='range',variable=variable,CTMM=CTMM)
  }
  else # append to current distribution
  {
    UD[[variable]] <- kde(data=data,H=UD$H,W=weights,alpha=error,grid=GRID,variable=variable,normalize=FALSE)

    #NAMES <- names(UD) # why is this being erased?
    #UD <- new.RS(UD,info=attr(data,"info"),type='range',variable=variable,CTMM=CTMM)
    #names(UD) <- NAMES # why is this necessary?
  }

  return(UD)
}
