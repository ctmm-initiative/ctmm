# location difference vector
# data is a list of two telemetry objects
# CTMM is a list of two ctmm fit objects corresponding to data
difference <- function(data,CTMM,t=NULL,...) { combine(data,CTMM,t=t,method="diff",...) }
midpoint <- function(data,CTMM,t=NULL,complete=FALSE,...) { combine(data,CTMM,t=t,complete=complete,method="mean",...) }

combine <- function(data,CTMM,t=NULL,complete=FALSE,method="diff",...)
{
  check.projections(data)
  INFO <- mean_info(data)
  INFO$identity <- paste0(method,"(",data[[1]]@info$identity,",",data[[2]]@info$identity,")")

  if(is.null(t))
  {
    # shared time range
    t1 <- max(data[[1]]$t[1],data[[2]]$t[1])
    t2 <- min(last(data[[1]]$t),last(data[[2]]$t))

    t <- c(t1,t2)
  }

  if(length(t)==2)
  {
    t1 <- t[1]
    t2 <- t[2]

    # shared times
    t <- c( data[[1]]$t , data[[2]]$t )
    t <- t[t>=t1]
    t <- t[t<=t2]
    t <- sort(t)
    t <- unique(t)

    if(!length(t))
    {
      warning("No overlapping times.")
      data <- data.frame(t=numeric(),x=numeric(),y=numeric())
      data <- new.telemetry(data,info=INFO)
      return(data)
    }
  }

  # predict over fine grid
  data[[1]] <- predict(data[[1]],CTMM[[1]],t=t)
  data[[2]] <- predict(data[[2]],CTMM[[2]],t=t)

  axes <- CTMM[[1]]$axes
  for(z in axes)
  {
    if(method=="diff")
    { data[[1]][[z]] <- data[[1]][[z]] - data[[2]][[z]] }
    else if(method=="mean")
    { data[[1]][[z]] <- data[[1]][[z]] + data[[2]][[z]] }
  }

  VAR <- DOP.LIST$horizontal$VAR
  COV <- DOP.LIST$horizontal$COV

  VARS <- COVS <- rep(FALSE,2)
  for(i in 1:2)
  {
    VARS[i] <- VAR %in% names(data[[1]])
    COVS[i] <- all(COV %in% names(data[[2]]))
  }

  if(!all(COVS)) # use only variance for speed later
  {
    v <- VAR
    # var of difference
    data[[1]][[v]] <- data[[1]][[v]] + data[[2]][[v]]

    COLS <- c('t',axes,VAR)
  }
  else # use covariance
  {
    for(i in 1:2)
    {
      if(!COVS[i]) # missing cov
      {
        data[[i]][[COV[1]]] <- data[[i]][[COV[3]]] <- data[[i]][[VAR]]
        data[[i]][[COV[2]]] <- 0
      }

      if(!VARS[i]) # missing var
      { data[[i]][[VAR]] <- (data[[i]][[COV[1]]] + data[[i]][[COV[3]]])/2 }
    }

    for(v in c(VAR,COV)) # var/cov of difference
    { data[[1]][[v]] <- data[[1]][[v]] + data[[2]][[v]] }

    COLS <- c('t',axes,VAR,COV)
  }
  data <- data.frame(data[[1]][,COLS])

  if(method=="mean")
  {
    data[,axes] <- data[,axes]/2
    data[,VAR] <- data[,VAR]/4
    if(all(COV %in% COLS)) { data[,COV] <- data[,COV]/4 }
  }

  # make this a telemetry object
  data <- new.telemetry(data,info=INFO)

  TYPE <- DOP.match(axes)

  # make sure UERE is correct
  data@UERE$UERE[] <- 1
  data@UERE$DOF[] <- Inf
  data@UERE$AICc[] <- -Inf
  data@UERE$Zsq[] <- 0
  data@UERE$VAR.Zsq[] <- 0
  data@UERE$N[] <- Inf

  colnames(data@UERE$UERE) <- TYPE
  colnames(data@UERE$DOF) <- TYPE
  names(data@UERE$AICc) <- TYPE
  names(data@UERE$Zsq) <- TYPE
  names(data@UERE$VAR.Zsq) <- TYPE
  names(data@UERE$N) <- TYPE

  rownames(data@UERE$UERE) <- "all"
  rownames(data@UERE$DOF) <- "all"

  if(complete) { data <- pseudonymize(data,tz=INFO$timezone,proj=INFO$projection,origin=EPOCH) }

  return(data)
}


# simple correlation test
proximity <- function(data,CTMM,t=NULL,GUESS=ctmm(error=TRUE),debias=TRUE,level=0.95,...)
{
  # difference vector with uncertainties
  data.diff <- difference(data,CTMM,t=t,...)

  if(!nrow(data.diff))
  {
    INF <- c(0,1,Inf)
    names(INF) <- NAMES.CI
    return(INF)
  }

  GUESS <- ctmm.guess(data.diff,CTMM=GUESS,interactive=FALSE)
  # # fit an autocorrelation model to the difference
  # GUESS <- CTMM
  # GUESS[[1]]$error <- GUESS[[2]]$error <- TRUE
  # # reflection from minus
  # GUESS[[2]]$sigma[1,2] <- GUESS[[2]]$sigma[2,1] <- -GUESS[[2]]$sigma[1,2]
  # GUESS[[2]]$sigma <- covm(GUESS[[2]]$sigma@.Data,axes=GUESS[[2]]$axes,isotropic=GUESS[[2]]$isotropic)
  M.diff <- ctmm.select(data.diff,GUESS,...)

  SIG <- var.covm(CTMM[[1]]$sigma,ave=TRUE) + var.covm(CTMM[[2]]$sigma,ave=TRUE)
  VAR.SIG <- axes2var(CTMM[[1]])['variance','variance'] + axes2var(CTMM[[2]])['variance','variance']

  SIG.diff <- var.covm(M.diff$sigma,ave=TRUE)
  VAR.SIG.diff <- axes2var(M.diff)['variance','variance']

  N <- 2*SIG^2/VAR.SIG
  if(debias)
  {
    DEN <- N/max(N-2,1) * 1/SIG
    VAR.DEN <- 2*DEN^2/max(N-4,1)
  }
  else
  {
    DEN <- 1/SIG
    VAR.DEN <- 2*DEN^2/N
  }
  # numerator-denominator covariance missing

  F.CI(SIG.diff,VAR.SIG.diff,DEN,VAR.DEN,level=level)
}


distances <- function(data,CTMM,t=NULL,level=0.95,...)
{
  timezone <- attr(data[[1]],'info')$timezone

  data <- difference(data,CTMM,t=t,...)
  n <- nrow(data)
  t <- data$t

  # estimate distances
  DISTS <- abs_data(data)
  M1 <- DISTS$M1
  M2 <- DISTS$M2
  DOF <- DISTS$DOF
  VAR <- DISTS$VAR

  # if no level, return point estimate and DOF
  if(is.null(level))
  { DISTS <- cbind(speed=M1,DOF=DOF,VAR=VAR) } # output v and chi DOF
  else
  {
    DISTS <- vapply(1:n,function(i){ chisq.ci(M2[i],DOF=DOF[i],level=level) },numeric(3)) # (3,n)
    DISTS <- sqrt(t(DISTS)) # (n,3)
    DISTS[,2] <- M1
  }

  DISTS <- as.data.frame(DISTS)
  DISTS$t <- t
  # include timestamps if possible
  if(!is.null(timezone))
  { DISTS$timestamp <- as.POSIXct(DISTS$t,tz=timezone,origin=EPOCH) }

  return(DISTS)
}
