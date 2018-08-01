# small sample size adjustment for ctmm.select to be more agressive
alpha.ctmm <- function(CTMM,alpha)
{
  z <- stats::qnorm(alpha)
  z <- sqrt(z^2 + (CTMM$AICc-CTMM$AIC))
  alpha <- 1-stats::pnorm(z)
  return(alpha)
}


###############
# keep removing uncertain parameters until AIC stops improving
ctmm.select <- function(data,CTMM,verbose=FALSE,level=0.99,IC="AICc",trace=FALSE,...)
{
  alpha <- 1-level
  trace2 <- if(trace) { trace-1 } else { 0 }

  UERE <- get.error(data,CTMM,flag=TRUE) # error flag only

  drift <- get(CTMM$mean)
  if(CTMM$mean=="periodic")
  {
    Nyquist <- CTMM$period/stats::median(diff(data$t))/2
    message("Nyquist frequency estimated at harmonic ",paste(Nyquist,collapse=" ")," of the period.")
  }

  # initial guess in case of pREML (better for optimization)
  get.mle <- function(FIT=CTMM)
  {
    MLE <- FIT
    if(!get("EMPTY",pos=MLE.env)) # will have been set from ctmm.fit first run
    {
      MLE <- get("MLE",pos=MLE.env)
      # check that structure is consistent
      if(is.null(MLE) || name.ctmm(MLE)!=name.ctmm(FIT)) { MLE <- FIT }
    }
    return(MLE)
  }

  # function to simplify complexity of models
  simplify <- function(M,par)
  {
    if("eccentricity" %in% par)
    {
      M$isotropic <- TRUE
      M$sigma <- covm(M$sigma,isotropic=TRUE)
    }
    if("circle" %in% par) { M$circle <- FALSE }

    return(M)
  }

  # consider a bunch of new models and update best model without duplication
  iterate <- function(GUESS)
  {
    # name the proposed models
    names(GUESS) <- sapply(GUESS,name.ctmm)
    # remove models alreadt fit
    RM <- which(names(GUESS) %in% names(MODELS))
    if(length(RM)) { GUESS <- GUESS[-RM] }

    # fit every model
    if(trace && length(GUESS)) { message("* Fitting models ",paste(names(GUESS),collapse=", ")) }
    #? should I run select here instead of fit ?
    GUESS <- lapply(GUESS,function(g){ctmm.fit(data,g,trace=trace2,...)})
    MODELS <<- c(MODELS,GUESS)

    # what is the new best model?
    OLD <<- CTMM
    CTMM <<- min.ctmm(MODELS)
  }

  ########################
  # PHASE 1: work our way up to complicated autocorrelation models
  # all of the features we need to fit numerically
  FEATURES <- id.parameters(CTMM,UERE=UERE)$NAMES
  # consider only features unnecessary "compatibility"
  FEATURES <- FEATURES[!(FEATURES=="area")]
  FEATURES <- FEATURES[!(FEATURES=="error")]
  FEATURES <- FEATURES[!grepl("tau",FEATURES)]

  # start with the most basic "compatible" model
  GUESS <- simplify(CTMM,FEATURES)
  if(trace) { message("* Fitting model ",name.ctmm(GUESS)) }
  TARGET <- CTMM
  CTMM <- ctmm.fit(data,GUESS,trace=trace2,...)
  MODELS <- list(CTMM)
  names(MODELS) <- sapply(MODELS,name.ctmm)

  OLD <- ctmm()
  while(!identical(CTMM,OLD))
  {
    GUESS <- list()
    MLE <- get.mle()

    # consider non-zero eccentricity
    if(("eccentricity" %in% FEATURES) && MLE$isotropic)
    {
      GUESS <- c(GUESS,list(MLE))
      n <- length(GUESS)
      GUESS[[n]]$isotropic <- FALSE
      # copy over target angle, but leave eccentricity zero to start (featureless)
      sigma <- attr(GUESS[[n]]$sigma,"par")
      sigma["angle"] <- attr(TARGET$sigma,"par")['angle']
      sigma <- covm(sigma,isotropic=FALSE,axes=TARGET$axes)
      GUESS[[n]]$sigma <- sigma
    }

    # consider circulation
    if(("circle" %in% FEATURES) && !MLE$circle)
    {
      GUESS <- c(GUESS,list(MLE))
      GUESS[[length(GUESS)]]$circle <- 2 * .Machine$double.eps * sign(TARGET$circle)
    }

    # consider no error?
    #

    # consider a bunch of new models and update best model without duplication
    iterate(GUESS)
  }

  #############################
  # PHASE 2: work our way down to simpler autocorrelation models & work our way up to more complex trend models
  OLD <- ctmm()
  # CTMM <- min.ctmm(MODELS)
  while(!identical(CTMM,OLD))
  {
    GUESS <- list()
    MLE <- get.mle()

    beta <- alpha.ctmm(CTMM,alpha)
    # consider if some timescales are actually zero
    CI <- confint.ctmm(CTMM,alpha=beta)

    if(length(CTMM$tau)==2)
    {
      Q <- CI["tau velocity",1]
      if(is.nan(Q) || (Q<=0))
      {
        GUESS <- c(GUESS,list(MLE))
        GUESS[[length(GUESS)]]$tau <- MLE$tau[-length(MLE$tau)]
      }
    }
    else if(length(CTMM$tau)==1)
    {
      Q <- CI["tau position",1]
      if(is.nan(Q) || (Q<=0))
      {
        GUESS <- c(GUESS,list(MLE))
        GUESS[[length(GUESS)]]$tau <- NULL
      }
    }

    # consider if there is no circulation
    if(CTMM$circle)
    {
      Q <- CI["circle",3]
      if(is.nan(Q) || (Q==Inf)) { GUESS <- c(GUESS,list(simplify(MLE,"circle"))) }
    }

    # consider if eccentricity is zero
    if(!CTMM$isotropic)
    {
      Q <- "eccentricity"
      Q <- stats::qnorm(beta/2,mean=CTMM$sigma@par[Q],sd=sqrt(CTMM$COV[Q,Q]))
      if(Q <= 0) { GUESS <- c(GUESS,list(simplify(MLE,"eccentricity"))) }
    }

    # consider if the mean could be more detailed
    GUESS <- c(GUESS,drift@refine(MLE))

    # consider a bunch of new models and update best model without duplication
    iterate(GUESS)
  }

  # return the best or return the full list of models
  if(verbose)
  {
    MODELS <- sort.ctmm(MODELS)
    return(MODELS)
  }
  else
  { return(CTMM) }
}


name.ctmm <- function(CTMM)
{
  # base model
  if(length(CTMM$tau)==2)
  { NAME <- "OUF" }
  else if(length(CTMM$tau)==1)
  { NAME <- "OU" }
  else if(length(CTMM$tau)==0)
  { NAME <- "IID" }

  # isotropy
  if(CTMM$isotropic)
  { NAME <- c(NAME,"isotropic") }
  else
  { NAME <- c(NAME,"anisotropic") }

  # circulation
  if(CTMM$circle)
  { NAME <- c(NAME,"circle") }

  # error
  if(CTMM$error)
  { NAME <- c(NAME,"error") }

  # mean
  drift <- get(CTMM$mean)
  NAME <- c(NAME,drift@name(CTMM))

  NAME <- paste(NAME,sep="",collapse=" ")
  return(NAME)
}

