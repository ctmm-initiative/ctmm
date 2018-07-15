# Try to fit a couple of different ways to ensure convergence to MLE
ctmm.refit <- function(data,CTMM,method="ML",COV=TRUE,control=list(),trace=FALSE,...)
{
  method <- match.arg(method,c("ML","pREML","pHREML","HREML","REML"))
  REML <- (method=="REML")
  int.method <- if(REML) { "REML" } else { "ML" }

  trace2 <- if(trace) { trace-1 } else { 0 }

  cores <- control$cores
  if(is.null(cores)){ cores <- 1 }
  cores <- resolveCores(cores,fast=TRUE)

  axes <- CTMM$axes
  MINS <- telemetry.mins(data,axes)
  dt <- MINS$dt
  df <- MINS$df
  dz <- MINS$dz

  # error variance per time
  ERROR <- get.error(data,CTMM,circle=TRUE)
  UERE <- attr(ERROR,"flag") # UERE flag

  # give partial fit a bit of parameter to make sure of subsequent fit attempt
  # unaccounted for error or excess coarsening can turn these parameters off
  perturb.ctmm <- function(M)
  {
    if(length(CTMM$tau) && (!length(M$tau) || M$tau[1]<=0)) { M$tau[1] <- dt/2 }
    if(length(CTMM$tau)>1 && (length(M$tau)<2 || M$tau[2]<=0)) { M$tau[2] <- min(M$tau[1],dt/2)/2 }

    if(CTMM$error && !M$error) { M$error <- dz/2 }

    if(CTMM$circle && !M$circle) { M$circle <- sign(CTMM$circle) * df/2 }

    return(M)
  }

  # just try to fit straight
  fit.all <- function(M,COV=FALSE,method=int.method)
  {
    if(trace) { message("Attempting direct fit.") }

    M <- ctmm.fit(data,M,method=method,COV=COV,control=control,trace=trace2)

    return(M)
  }

  # try to fit without error, but otherwise the same as M
  fit.no.error <- function(M)
  {
    if(trace) { message("Attempting to coarsen data and fit without errors.") }

    # semi-variance per time difference
    GUESS <- M
    GUESS$error <- FALSE

    svf.func <- svf.func(GUESS,moment=TRUE)$svf
    SVF <- svf.func(diff(data$t))
    # maximum semi-variance to/from a time
    SVF <- pmax(c(SVF[1],SVF),c(SVF,last(SVF)))

    # LOOP THIS???

    # quality index
    Q <- ERROR/SVF
    BAD <- which(Q>1) # subset of bad times to (potentially) remove
    QBAD <- Q[BAD] # bad quality values
    IND <- sort(QBAD,index.return=TRUE)$ix
    IND <- BAD[IND] # sorted indices of bad Q
    for(i in IND)
    {
      if(Q[i]<=1) { break }
      # quality index after deletion
      ###############!!!!!!!!!!!!!!!!!!!!!!
    }

    MAX <- which.max(Q)
    # thin data until errors are relatively small
    while(Q[MAX]>=1)
    {
      data <- data[-MAX,]
      Q <- Q[-MAX]
      ERROR <- ERROR[-MAX]

      SVF <- if(MAX==1) { diff(data$t[1:2]) }
      else if (MAX==length(data$t)) { max(diff(data$t[length(data$t) + -1:1])) }
      else { diff(data$t[MAX + 0:1]) }
      SVF <- svf.func(SVF)
      #
    }

    # error catch for above removing all data

    # coarsen data
    data <- (ERROR < SVF)
    data <- data[data,]

    # if circulation && !isotropic, then refit those
    if(CTMM$circle && !CTMM$isotropic)
    { FIT <- ctmm.refit(data,GUESS,method=int.method,COV=FALSE,control=control,trace=trace) }
    else
    { FIT <- ctmm.fit(data,GUESS,method=int.method,COV=FALSE,control=control,trace=trace2) }

    # find best combination of parameters - re-including error
    if(trace) { message("Combining fits.") }
    M <- permute.ctmm(data,list(M,FIT),UERE=UERE,cores=cores,REML=REML)

    # enforce minimum parameter so they still get fit attempt - in case original error estimate was too large or too small (tau->0)
    M <- perturb.ctmm(M)

    # re-fit with errors
    M <- ctmm.fit(data,M,method=int.method,COV=FALSE,control=control,trace=trace2)

    return(M)
  }

  # try to fit without circulation
  fit.no.circle <- function(M)
  {
    if(trace) { message("Attempting to fit without circulation.") }

    GUESS <- M
    GUESS$circle <- 0

    FIT <- ctmm.fit(data,GUESS,method=int.method,COV=FALSE,control=control,trace=trace2)
    M <- permute.ctmm(data,list(M,FIT),UERE=UERE,cores=cores,REML=REML)

    # re-fit with circulation
    M <- ctmm.fit(data,M,method=int.method,COV=FALSE,control=control,trace=trace2)

    return(M)
  }

  # try to fit without anisotropy
  fit.no.anisotropy <- function(M)
  {
    if(trace) { message("Attempting to fit without anisotropy.") }

    # setup isotropic version
    GUESS <- M
    GUESS$isotropic <- TRUE
    GUESS$sigma <- covm(GUESS$sigma,isotropic=TRUE,axes=M$axes)

    # fit without anisotropy
    FIT <- ctmm.fit(data,GUESS,method=int.method,COV=FALSE,control=control,trace=trace2)
    # re-include anisotropy
    FIT$sigma <- covm( c( attr(FIT$sigma,"par")[1] , attr(M$sigma,"par")[2:3] ) , isotropic=FALSE, axes=M$axes )
    # find best combination of parameters
    M <- permute.ctmm(data,list(M,FIT),UERE=UERE,cores=cores,REML=REML)

    # re-fit with anisotropy
    M <- ctmm.fit(data,M,method=int.method,COV=FALSE,control=control,trace=trace2)

    return(M)
  }

  ############
  FITS <- list()
  FITS[[1]] <- fit.all(CTMM)

  if(CTMM$error)
  { FITS[[2]] <- fit.no.error(CTMM) }
  else if(CTMM$circle && !CTMM$isotropic)
  { FITS[[2]] <- fit.no.circle(CTMM) }
  else
  {
    FIT <- ctmm.fit(data,FITS[[1]],method=method,COV=COV,control=control,trace=trace2)
    return(FIT)
  }

  if(!CTMM$isotropic) { FITS[[3]] <- fit.no.anisotropy(CTMM) }

  # permute all results
  if(trace) { message("Combining all fits.") }
  FIT <- permute.ctmm(data,FITS,UERE=UERE,cores=cores,REML=REML)
  # final attempt
  if(trace) { message("Final fit attempt.") }
  FIT <- ctmm.fit(data,FIT,method=method,COV=COV,control=control,trace=trace2)
  return(FIT)
}



# consider every combination of model parameters and return the best combination
permute.ctmm <- function(data,CTMM,UERE,cores=cores,REML=FALSE)
{
  P <- id.parameters(CTMM[[1]],UERE=UERE)
  NAMES <- P$NAMES
  n <- length(NAMES)

  m <- length(CTMM)

  # all possible parameter values
  P <- sapply(CTMM,function(M){ get.parameters(M,NAMES) }) # (n,m)
  # unique parameter values
  P <- lapply(1:n,function(i){ unique(P[i,]) })
  M <- sapply(P,length)
  m <- max(M)

  FITS <- list()
  # iterate through all m^n combinations
  # not sure of a faster way to do this than base-m counting
  for(i in 1:prod(M))
  {
    p <- P[,1]
    # resolve each digit in n-digit base-m number
    for(j in 1:n)
    {
      # last digit
      k <- i %% m

      FILL <- (k<=M[j])
      if(FILL) { p[j] <- P[[j]][k] }

      # remove last digit and shift to next digit
      i <- (i-k)/m
    }
    if(FILL) { FITS[[length(FITS)+1]] <- set.parameters(CTMM[[1]],p) }
  }

  # use ML or REML
  FITS <- plapply(FITS,function(M){ ctmm.loglike(data,M,verbose=TRUE,REML=REML) },cores=cores,FAST=TRUE)
  MAX <- sapply(FITS,function(M){ M$loglike })
  MAX <- which.max(MAX)
  # best model
  CTMM <- FITS[[MAX]]

  return(CTMM)
}

