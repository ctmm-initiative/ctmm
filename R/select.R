# small sample size adjustment for ctmm.select to be more agressive
alpha.ctmm <- function(CTMM,alpha)
{
  z <- stats::qnorm(alpha)
  z <- sqrt(z^2 + (CTMM$AICc-CTMM$AIC))
  alpha <- 1-stats::pnorm(z)
  return(alpha)
}


#########
get.MSPE <- function(CTMM,MSPE="position")
{
  if(!is.na(MSPE)) { MSPE <- CTMM$MSPE[MSPE] }
  else { MSPE <- Inf }
  return(MSPE)
}
########
get.IC <- function(CTMM,IC="AICc")
{
  if(!is.na(IC)) { IC <- CTMM[[IC]] }
  else { IC <- Inf }
  return(IC)
}

###############
# keep removing uncertain parameters until AIC stops improving
ctmm.select <- function(data,CTMM,verbose=FALSE,level=0.99,IC="AICc",MSPE="position",trace=FALSE,...)
{
  IC <- match.arg(IC,c("AICc","AIC","BIC",NA))
  MSPE <- match.arg(MSPE,c("position","velocity",NA))

  alpha <- 1-level
  trace2 <- if(trace) { trace-1 } else { 0 }
  IC <- match.arg(IC,c("AICc","AIC","BIC",NA))
  MSPE <- match.arg(MSPE,c("position","velocity",NA))

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

    if("range" %in% par)
    {
      M$tau[1] <- Inf
      M$range <- FALSE
    }

    return(M)
  }

  # consider a bunch of new models and update best model without duplication
  iterate <- function(DROP,REFINE=list())
  {
    # name the proposed models
    names(DROP) <- sapply(DROP,name.ctmm)
    names(REFINE) <- sapply(REFINE,name.ctmm)

    # remove models already fit
    DROP <- DROP[!(names(DROP) %in% names(MODELS))]
    REFINE <- REFINE[!(names(REFINE) %in% names(MODELS))]

    N <- length(DROP)
    M <- length(REFINE)
    GUESS <- c(DROP,REFINE)

    # fit every model
    if(trace && length(GUESS)) { message("* Fitting models ",paste(names(GUESS),collapse=", ")) }
    #? should I run select here instead of fit ?
    GUESS <- lapply(GUESS,function(g){ctmm.fit(data,g,trace=trace2,...)})

    MODELS <<- c(MODELS,GUESS)

    # check MSPE for improvement in REFINEd models
    # if(M>0 && !is.na(MSPE))
    # {
    #   if(N>0) { DROP <- GUESS[1:N] } else { DROP <- list() }
    #   REFINE <- GUESS[N + 1:M]
    #
    #   GOOD <- sapply(REFINE,function(M){get.MSPE(M,MSPE)}) <= get.MSPE(CTMM,MSPE)
    #   REFINE <- REFINE[GOOD]
    #
    #   GUESS <- c(DROP,REFINE)
    # }

    # what is the new best model?
    OLD <<- CTMM
    CTMM <<- min.ctmm(c(GUESS,list(CTMM)),IC=IC,MSPE=MSPE)
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
    if(length(CTMM$tau)==2 && !is.na(IC))
    {
      Q <- CI["tau velocity",1]
      if(is.nan(Q) || (Q<=0))
      {
        GUESS <- c(GUESS,list(MLE))
        GUESS[[length(GUESS)]]$tau <- MLE$tau[-length(MLE$tau)]
      }
    }
    else if(length(CTMM$tau)==1 && !is.na(IC))
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

    # consider if we can relax range residence (non-likelihood comparison only)
    if(CTMM$range && is.na(IC))
    { GUESS <- c(GUESS,list(simplify(MLE,"range"))) }

    # consider if the mean could be more detailed
    REFINE <- drift@refine(MLE)

    # consider a bunch of new models and update best model without duplication
    iterate(GUESS,REFINE)
  }

  # return the best or return the full list of models
  if(verbose)
  {
    MODELS <- sort.ctmm(MODELS,IC=IC,MSPE=MSPE)
    return(MODELS)
  }
  else
  { return(CTMM) }
}

################
name.ctmm <- function(CTMM,whole=TRUE)
{
  # base model
  if(length(CTMM$tau)==2)
  { if(CTMM$tau[1]<Inf) { NAME <- "OUF" } else { NAME <- "IOU" } }
  else if(length(CTMM$tau)==1)
  { if(CTMM$tau[1]<Inf) { NAME <- "OU" } else { NAME <- "BM" } }
  else if(length(CTMM$tau)==0)
  { NAME <- "IID" }

  # isotropy
  if(CTMM$isotropic)
  { NAME <- c(NAME,"isotropic") }
  else
  { NAME <- c(NAME,"anisotropic") }

  # circulation
  if(CTMM$circle)
  { NAME <- c(NAME,"circulation") }

  # error
  if(CTMM$error)
  { NAME <- c(NAME,"error") }

  # mean
  drift <- get(CTMM$mean)
  DNAME <- drift@name(CTMM)

  NAME <- paste(NAME,sep="",collapse=" ")

  if(whole && !is.null(DNAME))
  { NAME <- paste(NAME,DNAME) }
  else if(!whole)
  {
    if(is.null(DNAME)) { DNAME <- "stationary" }
    NAME <- c(NAME,DNAME)
  }

  return(NAME)
}

########
sort.ctmm <- function(x,decreasing=FALSE,IC="AICc",MSPE="position",flatten=TRUE,...)
{
  if(is.na(MSPE))
  { ICS <- sapply(x,function(m){get.IC(m,IC)}) }
  else if(is.na(IC))
  { ICS <- sapply(x,function(m){get.MSPE(m,MSPE)}) }
  if(is.na(MSPE) || is.na(IC))
  {
    IND <- sort(ICS,index.return=TRUE,decreasing=decreasing)$ix
    x <- x[IND]

    if(flatten) { return(x) }

    # structure the same as below
    if(is.na(IC)) { x <- list(x) }
    if(is.na(MSPE)) { x <- lapply(x,list) }

    return(x)
  }

  # model type names
  NAMES <- sapply(x,function(fit) name.ctmm(fit,whole=FALSE) )
  ACOV <- NAMES[1,]
  MEAN <- NAMES[2,]

  # group by ACOV
  ACOVS <- unique(ACOV)
  MEANS <- unique(MEAN)

  # partition into ACF-identical blocks for MSPE sorting, and then all-identical blocks for likelihood sorting
  y <- list()
  ICS <- numeric(length(ACOVS)) # ICs of best MSPE models
  for(i in 1:length(ACOVS))
  {
    # ACF-identical block
    SUB <- (ACOV==ACOVS[i])
    y[[i]] <- x[SUB]
    MEAN.SUB <- MEAN[SUB]
    MEANS.SUB <- unique(MEAN.SUB)

    z <- list()
    MSPES <- numeric(length(MEANS.SUB))
    for(j in 1:length(MEANS.SUB)) # sort exactly same models by IC
    {
      # all-identical block
      SUB <- (MEAN.SUB==MEANS.SUB[j])
      z[[j]] <- sort.ctmm(y[[i]][SUB],IC=IC,MSPE=NA)
      MSPES[j] <- get.MSPE(z[[j]][[1]],MSPE) # associate block with best's MSPE
    }
    IND <- sort(MSPES,index.return=TRUE,decreasing=decreasing)$ix
    y[[i]] <- do.call(c,z[IND]) # flatten to ACF blocks
    ICS[i] <- get.IC(y[[i]][[1]],IC) # associate block with best's IC
  }
  # sort blocks by IC and flatten
  IND <- sort(ICS,index.return=TRUE,decreasing=decreasing)$ix
  y <- y[IND]

  if(flatten) { y <- do.call(c,y) }

  return(y)
}

############
min.ctmm <- function(x,IC="AICc",MSPE="position",...)
{ sort.ctmm(x,IC=IC,MSPE=MSPE,...)[[1]] }

########
summary.ctmm.list <- function(object, IC="AICc", MSPE="position", units=TRUE, ...)
{
  IC <- match.arg(IC,c("AICc","AIC","BIC",NA))
  MSPE <- match.arg(MSPE,c("position","velocity",NA))

  N <- length(object)
  object <- sort.ctmm(object,IC=IC,MSPE=MSPE,flatten=FALSE)
  M <- length(object)

  # if(N==M) { MSPE <- NA } # don't need to sort MSPE
  # if(M==1) { IC <- NA } # don't need to sort IC
  object <- do.call(c,object)

  if(!is.na(IC))
  {
    ICS <- sapply(object,function(m){get.IC(m,IC)})

    # show relative IC
    ICS <- ICS - ICS[1]
    ICS <- cbind(ICS)
    colnames(ICS) <- paste0("\u0394",IC)
  }
  else { ICS <- NULL }

  if(!is.na(MSPE))
  {
    MSPES <- sapply(object,function(m){get.MSPE(m,MSPE)})

    # convert to meters/kilometers
    CNAME <- paste0("\u0394","RMSPE")
    MSPES <- sqrt(MSPES)
    MSPES <- MSPES - MSPES[1]
    MIN <- min(c(abs(MSPES[MSPES!=0]),Inf))
    UNIT <- unit(MIN,if(MSPE=="position"){"length"}else{"speed"},concise=TRUE,SI=!units)
    MSPES <- MSPES/UNIT$scale
    CNAME <- paste0(CNAME," (",UNIT$name,")")
    MSPES <- cbind(MSPES)
    colnames(MSPES) <- CNAME
  }
  else { MSPES <- NULL }

  ICS <- cbind(ICS,MSPES)
  rownames(ICS) <- names(object)

  DOF <- sapply(object,DOF.mean)
  METH <- sapply(object,function(m){m$method})

  DOF <- cbind(DOF)
  colnames(DOF) <- "DOF[mean]"

  if(FALSE) # only prints correctly in unicode locale (Windows R bug)
  {
    DOF <- data.frame(DOF,METH)
    colnames(DOF) <- c("DOF[mean]","method")
  }

  ICS <- cbind(ICS,DOF)

  return(ICS)
}

