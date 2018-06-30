# Try to fit a couple of different ways to ensure convergence to MLE
ctmm.refit <- function(data,CTMM,method="ML",trace=FALSE,control=list(),...)
{
  # if no error, we don't need to do this
  if(!CTMM$error && !(CTMM$circle && !CTMM$isotropic)) { return(ctmm.fit(data,CTMM,method=method,trace=trace,control=control,...)) }

  method <- match.arg(method,c("ML","pREML","pHREML","HREML","REML"))
  if(method!="REML") { REML <- FALSE } else { REML <- TRUE }
  int.method <- if(REML) { "REML" } else { "ML" }

  cores <- control$cores
  if(is.null(cores)){ cores <- 1 }
  cores <- resolveCores(cores,fast=TRUE)

  axes <- CTMM$axes
  MINS <- telemetry.mins(data,axes)
  dt <- MINS$dt
  df <- MINS$df
  dz <- MINS$dz

  # error variance per time
  ERROR <- get.error(data,M,circle=TRUE)
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

  # try to fit without error, but otherwise the same as M
  fit.no.error <- function(M)
  {
    if(trace) { message("Attempting to coarsen data and fit without errors.") }

    # semi-variance per time difference
    GUESS <- M
    GUESS$error <- FALSE

    SVF <- svf.func(GUESS,moment=TRUE)
    SVF <- SVF$svf(diff$data$t)
    # maximum semi-variance to/from a time
    SVF <- pmax(c(SVF[1],SVF),c(SVF,last(SVF)))

    # coarsen data
    DATA <- (ERROR < SVF)
    DATA <- data[DATA,]

    # fit to coarsened data without error
    FIT <- ctmm.fit(DATA,GUESS,method=int.method,COV=FALSE,control=control,trace=trace)
    # re-include error
    FIT$error <- M$error
    # find best combination of parameters
    M <- permute.ctmm(data,list(M,FIT),UERE=UERE,cores=cores,REML=REML)

    # enforce minimum parameter so they still get fit attempt
    M <- perturb.ctmm(M)

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
    FIT <- ctmm.fit(data,GUESS,method=int.method,COV=FALSE,control=control,trace=trace)
    # re-include anisotropy
    FIT$sigma <- covm( c( attr(FIT$sigma,"par")[1] , attr(M$sigma,"par")[2:3] ) , isotropic=FALSE, axes=M$axes )
    # find best combination of parameters
    M <- permute.ctmm(data,list(M,FIT),UERE=UERE,cores=cores,REML=REML)

    return(M)
  }

  # try to fit without circulation
  fit.no.circle <- function(M)
  {
    if(trace) { message("Attempting to fit without circulation.") }

    GUESS <- M
    GUESS$circle <- 0

    FIT <- ctmm.fit(data,GUESS,method=int.method,COV=FALSE,control=control,trace=trace)
    M <- permute.ctmm(data,list(M,FIT),UERE=UERE,cores=cores,REML=REML)

    return(M)
  }

  ############
  # ATTEMPT 1 -- coarsen data and turn error off
  if(CTMM$error)
  {
    CTMM <- fit.no.error(CTMM)
  }

  if(CTMM$circle) { }

  ############
  # ATTEMPT 2 -- consider isotropic model
  if(!CTMM$isotropic) { CTMM <- fit.no.anisotropy(CTMM) }

  # what about circulation
  # circle && !isotropic
  # circle && error
  # circle && error && !isotropic

  # final attempt
  #
  #
}



# consider every combination of model parameters and return the best combination
permute.ctmm <- function(data,CTMM,UERE,cores=cores,REML=FALSE)
{
  P <- id.parameters(CTMM[[1]],UERE=UERE)
  NAMES <- P$NAMES
  n <- length(NAMES)

  m <- length(CTMM)

  # all possible parameter combinations
  P <- sapply(CTMM,function(M){ get.parameters(M,NAMES) }) # (n,m)
  # delete redundant parameters
  RED <- logical(n)
  for(i in 1:n) { RED[i] <- all(P[i,]==P[i,1]) }
  if(all(RED)) { return(CTMM[[1]]) } # all models are the same
  P <- P[!RED,]
  n <- nrow(P)

  FITS <- list()
  # iterate through all m^n combinations
  for(i in 1:(m^n))
  {
    p <- P[,1]
    # resolve each digit in n-digit base-m number
    for(j in 1:n)
    {
      # last digit
      k <- i %% m
      p[j] <- P[j,k]
      # remove last digit and shift to next digit
      i <- (i-k)/m
    }
    FITS[[length(FITS)+1]] <- set.parameters(CTMM[[1]],p)
  }

  # use ML or REML
  FITS <- plapply(FITS,function(M){ ctmm.loglike(data,M,verbose=TRUE,REML=REML) },cores=cores,FAST=TRUE)
  MAX <- sapply(FITS,function(M){ M$loglike })
  MAX <- which.max(MAX)
  # best model
  CTMM <- FITS[[MAX]]

  return(CTMM)
}

