# small sample size adjustment for ctmm.select to be more agressive
alpha.ctmm <- function(CTMM,alpha)
{
  z <- stats::qnorm(alpha)
  z <- sqrt(z^2 + (CTMM$AICc-CTMM$AIC))
  alpha <- 1-stats::pnorm(z)
  if(!length(alpha)) { alpha <- 0 }
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

##################
# function to simplify complexity of models
simplify.ctmm <- function(M,par)
{
  if("minor" %in% par)
  {
    if("COV" %in% names(M))
    {
      M$COV <- axes2var(M)
      NAMES <- colnames(M$COV)
      VAR <- which(NAMES=="variance")
      colnames(M$COV)[VAR] <- rownames(M$COV)[VAR] <- "major"
    }

    M$isotropic <- TRUE
    M$sigma <- covm(M$sigma,isotropic=TRUE,axes=M$axes)
    par <- c(par,'angle')
  }

  if("major" %in% par)
  {
    M$isotropic <- TRUE
    M$sigma <- covm(0,isotropic=TRUE,axes=M$axes)
    M$tau <- NULL
    par <- c(par,c('circle','tau position','tau velocity','tau','omega'))
  }

  if("circle" %in% par)
  { M$circle <- FALSE }

  if("range" %in% par && length(M$tau))
  {
    # convert to diffusion matrix
    if(M$tau[1]>0) { M$sigma <- scale.covm(M$sigma,1/M$tau[1]) }
    M$tau[1] <- Inf
    M$range <- FALSE
  }

  # # reduce time-link parameters
  # if("timelink" %in% par) { M <- timelink.simplify(M) }

  # autocorrelation timescales can't be distinguished
  if("diff.tau" %in% par)
  {
    M$tau <- c(1,1)*mean(M$tau)
    M$omega <- FALSE
  }

  if('tau velocity' %in% par) { M$tau <- M$tau[1] }
  if(any(c('tau position','tau') %in% par)) { M$tau <- NULL }
  if('omega' %in% par) { M$omega <- FALSE }
  if('circle' %in% par) { M$circle <- FALSE }
  TEST <- grepl('error',par)
  if(any(TEST))
  {
    PAR <- par[TEST]
    PAR <- substr(PAR,nchar("error ?"),nchar(PAR))
    M$error[PAR] <- FALSE # why would you ever do this?
  }

  # re-identify features
  UERE.FIT <- paste("error",names(M$error)) %in% M$features
  names(UERE.FIT) <- names(M$error)

  M$features <- id.parameters(M,profile=FALSE,linear=FALSE,UERE.FIT=UERE.FIT)$NAMES

  # don't think we need this
  if("MLE" %in% names(M)) { M$MLE <- simplify.ctmm(M$MLE,par) }

  return(M)
}


#############
# function to make autocorrelation models more complex
complexify.ctmm <- function(M,par,TARGET)
{
  # consider non-zero eccentricity
  if(("minor" %in% par) && M$isotropic)
  {
    M$isotropic <- FALSE
    # copy over target angle, but leave eccentricity zero to start (featureless)
    sigma <- attr(M$sigma,"par")
    sigma["angle"] <- attr(TARGET$sigma,"par")['angle']
    sigma <- covm(sigma,isotropic=FALSE,axes=TARGET$axes)
    M$sigma <- sigma
  }

  # consider circulation
  if(("circle" %in% par) && !M$circle)
  {  M$circle <- 2 * .Machine$double.eps * sign(TARGET$circle) }

  # consider finite range
  if(("range" %in% par) && !M$range)
  {
    # convert from diffusion matrix
    M$tau[1] <- TARGET$tau[1]
    M$sigma <- scale.covm(M$sigma,TARGET$tau[1])
    M$range <- TRUE
  }

  # # increase time-link parameters
  # if("timelink" %in% par) { M <- timelink.complexify(M) }

  if("tau velocity" %in% par) { M$tau[2] <- TARGET$tau[2] }
  if("omega" %in% par) { M$omega <- TARGET$omega }

  # re-identify features
  UERE.FIT <- M$error>0 & (paste('error',names(M$error)) %in% M$features)
  # if(M$error && ('error' %in% M$features)) { UERE <- TRUE } else { UERE <- 3 } # calibrated or not?
  M$features <- id.parameters(M,profile=FALSE,linear=FALSE,UERE.FIT=UERE.FIT)$NAMES

  # don't think we need this
  if("MLE" %in% names(M)) { M$MLE <- complexify.ctmm(M$MLE,par,TARGET) }

  return(M)
}


# initial guess in case of pREML/HREML/pHREML (better for optimization)
get.mle <- function(FIT)
{
  if("MLE" %in% names(FIT))
  {
    axes <- FIT$axes
    MLE <- FIT$MLE
    FIT$MLE <- NULL
    # have parameters been altered after the fact?
    if(MLE$checksum==digest::digest(FIT,algo='md5')) { FIT <- MLE }
    FIT$axes <- axes # goes missing in MLE?
  }
  return(FIT)
}


###############
# keep removing uncertain parameters until AIC stops improving
ctmm.select <- function(data,CTMM,verbose=FALSE,level=1,IC="AICc",MSPE="position",trace=FALSE,cores=1,...)
{
  M <- ctmm.iterate(data,CTMM,verbose=TRUE,level=level,IC=IC,MSPE=MSPE,trace=trace,cores=cores,recurse=FALSE,TRYS=NULL,...)

  for(i in 1:length(M))
  {
    if(grepl('anisotropic',names(M)[i]))
    {
      ISO <- gsub(" anisotropic","",names(M)[i],fixed=TRUE)
      if(ISO %in% names(M)) { M[[i]]$ISO <- M[[ISO]] }
    }
  }

  if(verbose) { return(M) }
  else { return(M[[1]]) }
}

## internal function being called recursively by ctmm.select
# recurse: FALSE is outer-most call, TRUE are inner calls
# TRYS: models that we have tried
ctmm.iterate <- function(data,CTMM,verbose=FALSE,level=1,IC="AICc",MSPE="position",trace=FALSE,cores=1,recurse=FALSE,TRYS=NULL,...)
{
  LIST <- list(...)
  if(!is.null(LIST$control$message)) { message <- LIST$control$message }

  CV <- c("LOOCV","HSCV")
  IC <- match.arg(IC,c("AICc","AIC","BIC",CV,NA))
  MSPE <- match.arg(MSPE,c("position","velocity",NA))

  # CV: cross validation
  # requires expensive post-fit calculation of IC
  # can select between IID/OU/OUF and BM/IOU
  CV <- IC %in% CV

  alpha <- 1-level
  trace2 <- if(trace) { trace-1 } else { 0 }

  # accept list as completed candidates after first
  # for internal recursion, not for end users
  if(class(CTMM)[1]=="ctmm")
  {
    # TRYS <- NULL
    MODELS <- list()
  }
  else if(class(CTMM)[1]=="list")
  {
    TRYS <- c(TRYS,attr(CTMM,"attempted")) # models that we have tried
    MODELS <- CTMM[-1]
    CTMM <- CTMM[[1]]
  }
  CAND <- length(MODELS) # number of extra candidate models included for features
  if(CAND) { names(MODELS) <- sapply(MODELS,name.ctmm) }

  TYPE <- DOP.match(CTMM$axes)
  UERE.DOF <- attr(data,"UERE")$DOF[,TYPE]
  names(UERE.DOF) <- rownames(attr(data,"UERE")$DOF)
  UERE.FIT <- CTMM$error & !is.na(UERE.DOF) & UERE.DOF<Inf # will we be fitting any error parameters?

  # best non-zero variance estimates -- avoid zero collapse in intermediate model propagating to best model
  axes <- CTMM$axes
  AXES <- length(axes)
  ERROR <- CTMM$error # best non-zero error estimate
  VAR <- CTMM$sigma@par[1:AXES]
  fix.vars <- function(M)
  {
    if(length(M))
    {
      for(i in 1:length(M))
      {
        TEST <- M[[i]]$error <= .Machine$double.eps
        if(any(TEST)) { M[[i]]$error[TEST] <- ERROR[TEST] }
        M.VAR <- M[[i]]$sigma@par
        for(j in 1:AXES) { if(M.VAR[j]<=.Machine$double.eps) { M.VAR[j] <- VAR[j] } }
        M[[i]]$sigma <- covm(M.VAR,axes=axes,isotropic=M[[i]]$isotropic)
      }
    }
    return(M)
  }
  update.fix <- function(CTMM)
  {
    TEST <- CTMM$error > .Machine$double.eps
    if(any(TEST)) { ERROR[TEST] <<- CTMM$error[TEST] }
    for(i in 1:AXES) { if(CTMM$sigma@par[i]>.Machine$double.eps) { VAR[i] <<- CTMM$sigma@par[i] } }
  }

  if(CTMM$mean=="periodic")
  {
    Nyquist <- CTMM$period/stats::median(diff(data$t))/2
    message("Nyquist frequency estimated at harmonic ",paste(Nyquist,collapse=" ")," of the period.")
  }

  # consider a bunch of new models and update best model without duplication
  iterate <- function(DROP,REFINE=list(),phase=1)
  {
    # make sure everything tries errors
    DROP <- fix.vars(DROP)
    REFINE <- fix.vars(REFINE)

    # name the proposed models
    names(DROP) <- sapply(DROP,name.ctmm)
    names(REFINE) <- sapply(REFINE,name.ctmm)

    # remove models already attempted
    DROP <- DROP[!(names(DROP) %in% TRYS)]
    REFINE <- REFINE[!(names(REFINE) %in% TRYS)]

    # copy over best initial parameter guess for refined drops
    if(!drift.is.stationary(CTMM) && length(DROP))
    {
      for(i in 1:length(DROP))
      {
        for(j in 1:length(MODELS))
        {
          # take last/best model of same covariance structure
          F1 <- DROP[[i]]$features
          F2 <- MODELS[[j]]$features
          if(length(F1)==length(F2) && all(F1==F2))
          {
            # copy over covariance parameters
            DROP[[i]] <- copy.parameters(DROP[[i]],get.mle(MODELS[[j]]))
            break # out of MODELS loop
          }
        } # end MODELS loop
      } # end DROP loop
    } # end refined drop adjustments

    GUESS <- c(DROP,REFINE)
    # add to TRYS before potential feature collapse
    TRYS.OLD <- TRYS # don't count new attempts for next select
    if(length(GUESS))
    {
      TRYS <<- c(TRYS,names(GUESS))
      TRYS <<- unique(TRYS)
    }

   if(length(GUESS)>1) # keep selecting (recursive)
    {
      # fit every model
      if(trace && FALSE) { message("* Fitting models ",paste(names(GUESS),collapse=", "),".") }
      # GUESS <- plapply(GUESS,function(g){M<-c(list(g),MODELS); attr(M,"attempted")<-TRYS.OLD; ctmm.select(data,M,verbose=TRUE,level=0,IC=IC,MSPE=MSPE,trace=trace,...)},cores=cores)
      GUESS <- plapply(GUESS,function(g){M<-c(list(g),MODELS); ctmm.iterate(data,M,verbose=TRUE,level=0,IC=IC,MSPE=MSPE,trace=trace,recurse=TRUE,TRYS=TRYS.OLD,...)},cores=cores)
      GUESS[sapply(GUESS,is.null)] <- NULL # delete collapsed repeats
      if(length(GUESS)) # concatenate list of lists
      {
        ATTEMPT <- do.call(c, lapply(GUESS,function(g){attr(g,"attempted")}) ) # do.call-lapply in case of NULL
        ATTEMPT[sapply(ATTEMPT,is.null)] <- NULL # delete NULL results... ? does this still happen ?
        TRYS <- c(TRYS,ATTEMPT)
        TRYS <- unique(TRYS)
        GUESS <- do.call(c,GUESS)
      }
    }
    else if(length(GUESS)==1) # fit last model
    {
      if(trace) { message("* Fitting model ",names(GUESS)) }
      GUESS[[1]] <- ctmm.fit(data,GUESS[[1]],trace=trace2,...)

      # cross-validate
      if(CV)
      {
        if(trace) { message("** Cross validating model ",names(GUESS)[1]) }
        GUESS[[1]][[IC]] <- do.call(IC,list(data=data,CTMM=GUESS[[1]],cores=cores,...))
      }
    }
    # corrected name after potential feature collapse
    names(GUESS) <- sapply(GUESS,name.ctmm)
    # add to TRYS after potential feature collapse
    if(length(GUESS))
    {
      TRYS <<- c(TRYS,names(GUESS))
      TRYS <<- unique(TRYS)
    }
    # TRYS can be longer than MODELS

    MODELS <<- c(GUESS,MODELS)

    # only sort newer models, not input candidate models, which we will not return
    N <- length(MODELS) - CAND
    if(N>1) { MODELS[1:N] <<- sort.ctmm(MODELS[1:N],IC=IC,MSPE=MSPE) }

    # what is the new best model?
    OLD <<- CTMM
    CTMM <<- MODELS[[1]]
    update.fix(CTMM) # update best non-zero error estimate
    # CTMM <<- min.ctmm(c(GUESS,list(CTMM)),IC=IC,MSPE=MSPE)
  } # end iterate()

  ########################
  # PHASE 0: pair down to essential features for 'compatibility'
  # all of the features we need to fit numerically
  FEATURES <- id.parameters(CTMM,UERE.FIT=UERE.FIT)$NAMES
  # consider only features unnecessary for Stein "compatibility"
  FEATURES <- FEATURES[!(FEATURES=="major")]
  FEATURES <- FEATURES[!grepl("error",FEATURES)]
  FEATURES <- FEATURES[!grepl("tau",FEATURES)]
  FEATURES <- FEATURES[!(FEATURES=="omega")]
  # if((CV || is.na(IC)) && CTMM$range && length(CTMM$tau)) { FEATURES <- c(FEATURES,"range") }

  # # clean up non-parametrics
  # if(!is.null(CTMM$timelink) && CTMM$timelink!="identity")
  # {
  #   FEATURES <- FEATURES[!grepl("timelink",FEATURES)]
  #   FEATURES <- c(FEATURES,"timelink")
  # }

  TARGET <- CTMM
  # start with the most basic "compatible" model, if not included in candidates
  GUESS <- simplify.ctmm(CTMM,FEATURES)
  NAME <- name.ctmm(GUESS)
  if(NAME %in% TRYS) # we have tried the paired down model
  { CTMM <- MODELS[[1]] }
  else # fit paired down model
  {
    if(trace) { message("* Fitting model ",NAME) }
    TRYS <- c(NAME,TRYS)
    CTMM <- ctmm.fit(data,GUESS,trace=trace2,...)
    if(CV)
    {
      if(trace) { message("** Cross validating model ",name.ctmm(GUESS)) }
      CTMM[[IC]] <- do.call(IC,list(data=data,CTMM=CTMM,cores=cores,...))
    }
    TRYS <- c(name.ctmm(CTMM),TRYS)
    TRYS <- unique(TRYS)
    MODELS <- c(list(CTMM),MODELS) # update initial guess to fit
    update.fix(CTMM)
  }
  names(MODELS) <- sapply(MODELS,name.ctmm)

  ##############################
  # PHASE 1: work our way up to complicated autocorrelation models
  OLD <- ctmm()
  while(!identical(CTMM,OLD))
  {
    MLE <- get.mle(CTMM)
    GUESS <- lapply(FEATURES,function(feat){complexify.ctmm(MLE,feat,TARGET)})

    # consider a bunch of new models and update best model without duplication
    iterate(GUESS)
  }

  #############################
  # PHASE 2: work our way down to simpler autocorrelation models & work our way up to more complex trend models
  # we only do phase 2 if(level), so we can call ctmm.select recursively to build up features again
  OLD <- ctmm()
  # CTMM <- min.ctmm(MODELS)
  while(level && !identical(CTMM,OLD))
  {
    GUESS <- list()
    MLE <- get.mle(CTMM)
    beta <- alpha.ctmm(CTMM,alpha)

    VARS <- dimnames(CTMM$COV)[[1]] # could be NULL

    # consider if smallest timescale is zero
    CI <- confint.ctmm(CTMM,alpha=beta)
    if(length(CTMM$tau)==2 && !is.na(IC)) # OUX -> OUx
    {
      if("tau velocity" %in% VARS) # OUF -> OU / IOU -> BM
      {
        Q <- CI["tau velocity",1]
        if(level==1 || is.na(Q) || (Q<=0))
        { GUESS <- c(GUESS,list(simplify.ctmm(MLE,'tau velocity'))) }
      }
      else if("omega" %nin% VARS) # OUf -> OU, IID
      {
        Q <- CI["tau",1]
        if(level==1 || is.na(Q) || (Q<=0))
        {
          GUESS <- c(GUESS,list(simplify.ctmm(MLE,'tau'))) # OUf -> IID
          GUESS <- c(GUESS,list(simplify.ctmm(MLE,'tau velocity'))) # OUf -> OU
        }
      }
      else # OUO -> OUf
      {
        Q <- 1/CI["tau period",3]
        if(level==1 || is.na(Q) || (Q<=0))
        { GUESS <- c(GUESS,list(simplify.ctmm(MLE,'omega'))) }
      }
    } # end # OUX -> OUx
    else if(CTMM$range && length(CTMM$tau)==1 && !is.na(IC)) # OU -> IID
    {
      Q <- CI["tau position",1]
      if(level==1 || is.na(Q) || (Q<=0))
      { GUESS <- c(GUESS,list(simplify.ctmm(MLE,'tau position'))) }
    } # end # OU -> IID

    # can autocorrelation timescales be distinguished?
    if(all(c('tau position','tau velocity') %in% VARS) || "omega" %in% VARS) # isn't omega redundant?
    {
      TEMP <- get.taus(CTMM,zeroes=TRUE)
      nu <- TEMP$f.nu[2] # frequency/difference
      J <- TEMP$J.nu.tau[2,] # Jacobian of nu WRT canonical parameters
      Q <- TEMP$tau.names
      Q <- c(J %*% CTMM$COV[Q,Q] %*% J) # variance of nu
      Q <- ci.tau(nu,Q,alpha=beta)[1]

      if(level==1 || is.na(Q) || Q<=0 || is.na(IC))
      { GUESS <- c(GUESS,list(simplify.ctmm(MLE,"diff.tau"))) }
    }
    else if("tau" %in% VARS) # try other side if boundary if chosen model is critically damped
    {
      # try overdamped
      TEMP <- TARGET
      if(length(TEMP$tau)<2 || !TEMP$tau[2] || TEMP$tau[2] >= MLE$tau[1]) { TEMP$tau <- MLE$tau * exp(c(1,-1)*sqrt(.Machine$double.eps)) }
      GUESS <- c(GUESS,list(complexify.ctmm(MLE,"tau velocity",TEMP)))

      # try underdamped
      TEMP <- TARGET
      if(!TEMP$omega) { TEMP$omega <- sqrt(.Machine$double.eps) }
      GUESS <- c(GUESS,list(complexify.ctmm(MLE,"omega",TEMP)))
    }
    else if("tau position" %in% VARS && level==1) # OU -> OUf (bimodal likelihood)
    { GUESS <- c(GUESS,list(simplify.ctmm(MLE,"diff.tau"))) }

    # consider if there is no circulation
    if("circle" %in% VARS)
    {
      Q <- CI["circle",3]
      if(level==1 || is.na(Q) || (Q==Inf) || is.na(IC))
      { GUESS <- c(GUESS,list(simplify.ctmm(MLE,"circle"))) }
    }

    # consider if eccentricity is zero --- log(major)==log(minor)
    if(!CTMM$isotropic)
    {
      Q <- c("major","minor")
      GRAD <- c(1/CTMM$sigma@par['major'],-1/CTMM$sigma@par['minor'])
      SD <- ifelse(all(Q %in% CTMM$features),sqrt(c(GRAD %*% CTMM$COV[Q,Q] %*% GRAD)),Inf) # variance could collapse early
      Q <- log(CTMM$sigma@par['major']/CTMM$sigma@par['minor'])
      Q <- stats::qnorm(beta/2,mean=Q,sd=SD)
      if(level==1 || is.na(Q) || Q<=0 || is.na(IC))
      { GUESS <- c(GUESS,list(simplify.ctmm(MLE,"minor"))) }
    }

    # is the animal even moving?
    if(!CTMM$sigma@par['major'] && any(grepl("error",VARS)))
    { GUESS <- c(GUESS,list(simplify.ctmm(MLE,"major"))) }

    # consider if we can relax range residence (non-likelihood comparison only)
    if(CTMM$range && length(CTMM$tau) && (is.na(IC) || CV))
    { GUESS <- c(GUESS,list(simplify.ctmm(MLE,"range"))) }
    else if(!CTMM$range && (is.na(IC) || CV) && TARGET$tau[1]<Inf)
    { GUESS <- c(GUESS,list(complexify.ctmm(MLE,"range",TARGET))) }

    # consider if time link could be more or less detailed
    if(!is.null(TARGET$timelink) && TARGET$timelink!="identity")
    { GUESS <- c(GUESS,list(timelink.complexify(MLE),timelink.simplify(MLE))) }

    # What was special about these cases?
    # consider if the mean could be more detailed or less detailed
    REFINE <- c(drift.complexify(MLE),drift.simplify(MLE))

    # consider a bunch of new models and update best model without duplication
    iterate(GUESS,REFINE,phase=2)
  }

  # remove any input candidate models
  N <- length(MODELS) - CAND
  MODELS <- MODELS[1:N]

  # return the best or return the full list of models
  if(verbose)
  {
    if(N>1) { MODELS <- sort.ctmm(MODELS,IC=IC,MSPE=MSPE) }
    NAMES <- sapply(MODELS,name.ctmm) -> names(MODELS)

    if(recurse) # store all model attempts in attribute
    { attr(MODELS,'attempted') <- TRYS }
    else # finishing stuff
    {
      # remove redundant models if outer most run
      if(N>1)
      {
        IN <- rep(TRUE,N)
        for(i in 2:N) { if(NAMES[i] %in% NAMES[1:(i-1)]) { IN[i] <- FALSE } }
        MODELS <- MODELS[IN]
      }

      # store model selection criterion
      attr(MODELS,"IC") <- IC
      attr(MODELS,"MSPE") <- MSPE
    }

    return(MODELS)
  } # end verbose
  else
  { return(CTMM) }
}


################
name.ctmm <- function(CTMM,whole=TRUE)
{
  if(is.null(CTMM)) { return(NULL) }

  FEATURES <- CTMM$features

  # base model
  tau <- CTMM$tau
  if(length(tau)==2)
  {
    if(tau[1]==Inf) { NAME <- "IOU" }
    else if(tau[1]>tau[2]) { NAME <- "OUF" }
    else if(CTMM$omega) { NAME <- "OU\u03A9" } # underdamped
    else { NAME <- "OUf" } # identical timescales
  }
  else if(length(tau)==1)
  { if(tau[1]<Inf) { NAME <- "OU" } else { NAME <- "BM" } }
  else if(length(tau)==0)
  {
    POS <- "sigma" %in% names(CTMM) && CTMM$sigma@par['major']<=0
    if(POS || "major" %in% FEATURES || "sigma" %nin% names(CTMM))
    { NAME <- "IID" }
    else
    { NAME <- "inactive" }
  }

  # isotropy
  isotropic <- CTMM$isotropic
  if(length(CTMM$axes)>1 && !is.null(isotropic))
  {
    if(all(!isotropic))
    { NAME <- c(NAME,"anisotropic") }
    else if(!isotropic[1])
    { NAME <- c(NAME,"anisotropic-\u03C3") }
    else if(length(isotropic)>1 && !isotropic[2])
    { NAME <- c(NAME,"anisotropic-\u03BC") }
    # else
    # { NAME <- c(NAME,"isotropic") }
  }

  # circulation
  if(CTMM$circle || "circle" %in% FEATURES)
  { NAME <- c(NAME,"circulation") }

  # error
  if(any(CTMM$error>0) || any(grepl("error",FEATURES)))
  { NAME <- c(NAME,"error") }

  # timelink
  if(!is.null(CTMM$timelink) && CTMM$timelink!="identity" && length(CTMM$timelink.par))
  {
    fn <- get(paste0(CTMM$timelink,".timelink.name"))
    NAME <- c(NAME,fn(CTMM))
  }

  # data link
  link <- get.link(CTMM)
  if(link$name!="identity") { NAME <- c(NAME,link$name) }

  # mean
  DNAME <- drift.name(CTMM)

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
sort.ctmm <- function(x,decreasing=FALSE,IC="AICc",MSPE="position",flatten=TRUE,INF=FALSE,...)
{
  CV <- c("LOOCV","HSCV")
  CV <- IC %in% CV

  if(is.na(MSPE))
  { ICS <- sapply(x,function(m){get.IC(m,IC)}) }
  else if(is.na(IC))
  { ICS <- sapply(x,function(m){get.MSPE(m,MSPE)}) }

  if(is.na(MSPE) || is.na(IC))
  {
    ICS <- nant(ICS,Inf)

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

  # BM/IOU log-likelihood is infinitely lower than OU/OUF log-likelihood
  RANGE <- sapply(y,function(Y){Y[[1]]$range})
  if(!is.na(IC) && !CV && any(RANGE) && any(!RANGE))
  {
    if(INF) { for(i in which(!RANGE)) { for(j in 1:length(y[[i]])) { y[[i]][[j]][[IC]] <- Inf } } }
    y <- c( y[RANGE] , y[!RANGE] )
  }

  if(flatten) { y <- do.call(c,y) }

  return(y)
}

############
min.ctmm <- function(x,IC="AICc",MSPE="position",...)
{ sort.ctmm(x,IC=IC,MSPE=MSPE,...)[[1]] }


########
summary.ctmm.list <- function(object, IC=NULL, MSPE=NULL, units=TRUE, ...)
{
  # use what ctmm.select used
  if(is.null(IC)) { IC <- attr(object,"IC") }
  if(is.null(MSPE)) { MSPE <- attr(object,"MSPE") }

  # defaults for manual lists
  if(is.null(IC)) { IC <- "AICc" }
  if(is.null(MSPE)) { MSPE <- "position" }

  CV <- c("LOOCV","HSCV")
  IC <- match.arg(IC,c("AICc","AIC","BIC",CV,NA))
  MSPE <- match.arg(MSPE,c("position","velocity",NA))
  CV <- IC %in% CV

  N <- length(object)
  object <- sort.ctmm(object,IC=IC,MSPE=MSPE,flatten=FALSE,INF=TRUE)
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
    MIN <- which.min(MSPES)
    MSPES <- sqrt(MSPES)
    if(MSPES[1]<Inf) { MSPES <- MSPES - MSPES[MIN] }
    # MIN <- min(c(MSPES[MSPES>0],Inf))
    AVE <- stats::median(MSPES[MSPES>0])
    if(MSPE=="position") { UNIT <- "length" } else { UNIT <- "speed" }
    UNIT <- unit(AVE,UNIT,concise=TRUE,SI=!units)
    MSPES <- MSPES/UNIT$scale
    CNAME <- paste0(CNAME," (",UNIT$name,")")
    MSPES <- cbind(MSPES)
    colnames(MSPES) <- CNAME
  }
  else { MSPES <- NULL }

  ICS <- cbind(ICS,MSPES)
  rownames(ICS) <- names(object)

  if(is.na(MSPE))
  {
    DOF <- sapply(object,DOF.mean)
    DOF <- cbind(DOF)
    colnames(DOF) <- "DOF[mean]"
  }
  else if(MSPE=="position")
  {
    DOF <- sapply(object,DOF.area)
    DOF <- cbind(DOF)
    colnames(DOF) <- "DOF[area]"
  }
  else if(MSPE=="velocity")
  {
    DOF <- sapply(object,DOF.speed)
    DOF <- cbind(DOF)
    colnames(DOF) <- "DOF[speed]"
  }

  METH <- sapply(object,function(m){m$method})
  if(FALSE) # only prints correctly in unicode locale (Windows R bug)
  {
    DOF <- data.frame(DOF,METH)
    colnames(DOF) <- c("DOF[mean]","method")
  }

  ICS <- cbind(ICS,DOF)

  return(ICS)
}

