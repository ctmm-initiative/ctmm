# Try to fit a couple of different ways to ensure convergence to MLE
ctmm.refit <- function(data,CTMM,method="ML",trace=FALSE,control=list(),...)
{
  # if no error, we don't need to do this
  if(CTMM$error==FALSE && CTMM$circle==FALSE) { return(ctmm.fit(data,CTMM,method=method,trace=trace,control=control,...)) }

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

  ############
  # ATTEMPT 1 -- coarsen data and turn error off
  if(trace) { message("Attempting to coarsen data and fit without errors.") }

  # error variance per time
  ERROR <- get.error(data,CTMM,circle=TRUE)
  UERE <- attr(ERROR,"flag") # UERE flag

  # semi-variance per time difference
  GUESS <- CTMM
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

  # find best combination of parameters
  CTMM <- permute.ctmm(data,list(CTMM,FIT),UERE=UERE,cores=cores,REML=REML)

  # enforce minimum parameter so they still get fit attempt
  CTMM <- perturb.ctmm(CTMM,df=df,dt=dt,dz=dz)

  ############
  # ATTEMPT 2 -- consider isotropic model

  # find best combination of parameters

  # final attempt
}


# give model a bit of parameter to make sure of subsequent fit attempt
perturb.ctmm <- function(CTMM,df=0,dt=0,dz=10)
{

}


# consider every combination of model parameters and return the best combination
permute.ctmm <- function(data,CTMMS,UERE,cores=cores,REML=FALSE)
{
  P <- id.parameters(CTMM[[1]],UERE=UERE)
  NAMES <- P$NAMES
  n <- length(NAMES)

  m <- length(CTMM)

  # all possible parameter combinations
  P <- sapply(CTMM,function(M){ get.parameters(M,NAMES) }) # (n,m)
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

