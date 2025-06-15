revisitation <- function(data,UD,error=0.001,...)
{ npr(data,UD,variable="revisitation",normalize=TRUE,error=error,...) }


# TODO Occurrence NPRs ##
npr <- function(data,UD,variable="speed",normalize=FALSE,error=0.001,...)
{
  # arbitrary variables allowed
  # variable <- match.arg(variable,c('revisitation','speed'))

  info <- attr(data,"info")
  type <- UD@type

  CTMM <- UD@CTMM
  axes <- CTMM$axes
  # res <- 10

  if(length(CTMM$tau)<2)
  {
    if(variable=='revisitation')
    {
      warning("Revisitation estimation requires a continuous-velocity model.")
      if(normalize)
      { return(UD) }
      else
      { stop("Unnormalized revisitation rate cannot be estimated.") }
    }
    else if(variable=='speed')
    { stop("Speed estimation requires a continuous-velocity model.") }
  }

  # !!! FIX THIS
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

  if(variable=="revisitation") # weighted distribution
  {
    weights <- pi * weights * data$speed
    weight <- sum(weights)

    NEW <- kde(data=data,H=UD$H,W=weights,alpha=error,grid=GRID) # normalized UD object

    UD$weight <- weight # total weight for means
    UD$rate <- weight * sum(UD$PDF*NEW$PDF) * prod(NEW$dr) # <dF/dR>

    UD$PDF <- NEW$PDF
    UD$CDF <- NEW$CDF
  }
  else  # NR regression surface
  {
    NEW <- kde(data=data,H=UD$H,W=weights,alpha=error,grid=GRID,variable=variable,normalize=normalize)

    if(normalize) # generated a new distribution
    {
      UD$PDF <- NEW$NPR * (NEW$P/sum(NEW$NPR)) # normalize surface
      UD$CDF <- pmf2cdf(UD$PDF)
      UD$PDF <- UD$PDF / prod(NEW$dr)
    }
    else # append expectation layer to current object
    { UD[[variable]] <- NEW$NPR }
  }

  return(UD)
}
