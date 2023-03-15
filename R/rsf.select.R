rsf.select <- function(data,UD,R=list(),formula=NULL,verbose=FALSE,IC="AICc",trace=TRUE,...)
{
  CTMM <- UD@CTMM
  isotropic <- CTMM$isotropic
  SISO <- ifelse(isotropic,"isotropic","anisotropic")
  M <- list() # list of candidate models
  if(!isotropic)
  {
    if("ISO" %in% names(CTMM))
    { ISO <- CTMM$ISO }
    else
    {
      ISO <- simplify.ctmm(CTMM,'minor')
      if(trace) { message("Fitting isotropic autocorrelation model.") }
      ISO <- ctmm.fit(data,ISO,trace=max(trace-1,0))
    }
    UD@CTMM <- ISO
    UD$DOF.area <- DOF.area(ISO)
    M <- rsf.select(data,UD,R=R,formula=formula,verbose=verbose,IC=IC,trace=trace,...)
    return(M)

    UD@CTMM <- CTMM
  }

  if(length(R))
  {
    RVARS <- names(R)
    FORMULA <- !is.null(formula)

    if(!FORMULA)
    { TERMS <- RVARS }
    else
    {
      VARS <- all.vars(formula)
      DVARS <- VARS[ VARS %nin% RVARS ]

      for(D in DVARS) { data[[D]] <- as.numeric(data[[D]]) } # model.matrix will rename otherwise

      # this doesn't work with poly(), etc.
      # TERMS <- attr(stats::terms(formula),"term.labels")
      # dummy data for model parameter names
      DATA <- data.frame(data)[1:2,]
      DATA[RVARS] <- as.list(rep(0,length(RVARS)))
      TERMS <- colnames(stats::model.matrix(formula,data=DATA))
      TERMS <- TERMS[TERMS!="(Intercept)"]
    }
  }
  DIM <- length(TERMS)

  ON <- array(FALSE,DIM)
  names(ON) <- TERMS

  ########################
  # fit without covariates first
  NAME <- paste(SISO,"~ 0")
  if(trace) { message("Fitting RSF model ",NAME) }
  M[[NAME]] <- rsf.fit(data,UD,R=R,formula=~0,trace=max(trace-1,0),...)
  # re-use estimates for next fit
  UD@CTMM <- M[[NAME]]

  ################
  # build up phase
  while(any(!ON))
  {
    # which terms are presently off
    try <- which(!ON)
    TRYS <- array(ON,c(DIM,length(try)))
    TRYS <- t(TRYS) # [off->on,all]
    for(i in 1:length(try)) { TRYS[i,try[i]] <- TRUE }
    # models that we will try by turning one term on
    FORM <- sapply(1:length(try),function(i){paste("~",paste(TERMS[TRYS[i,]],collapse="+"))})
    NAMES <- paste(SISO,FORM)
    FORM <- lapply(FORM,stats::as.formula)
    # don't try the same model twice
    SUB <- NAMES %nin% names(M)
    NAMES <- NAMES[SUB]
    FORM <- FORM[SUB]

    # fit all
    NEW <- list()
    for(i in 1%:%length(NAMES))
    {
      if(trace) { message("Fitting RSF model ",NAMES[i]) }
      NEW[[NAMES[i]]] <- rsf.fit(data,UD,R=R,formula=FORM[[i]],trace=max(trace-1,0),...)
    }
    M <- c(M,NEW)

    # best RSF parameters so far (going up)
    ICS <- sapply(NEW,function(m){m[[IC]]})
    i <- which.min(ICS)
    beta <- NEW[[i]]$beta
    ON <- TERMS %in% names(beta)
    # copy for initial guess in fitting
    # UD@CTMM$beta <- beta
    UD@CTMM <- NEW[[i]]
  }

  ###################
  # pair down phase
  while(any(ON))
  {
    # which terms are presently off
    try <- which(ON)
    TRYS <- array(ON,c(DIM,length(try)))
    TRYS <- t(TRYS) # [off->on,all]
    for(i in 1:length(try)) { TRYS[i,try[i]] <- FALSE }
    # models that we will try by turning one term on
    FORM <- sapply(1:length(try),function(i){paste("~",paste(TERMS[TRYS[i,]],collapse="+"))})
    FORM[FORM=="~ "] <- "~ 0" # make a valid formula string
    NAMES <- paste(SISO,FORM)
    FORM <- lapply(FORM,stats::as.formula)
    # don't try the same model twice
    SUB <- NAMES %nin% names(M)
    NEW.NAMES <- NAMES[SUB]
    FORM <- FORM[SUB]

    # fit all
    NEW <- list()
    for(i in 1%:%length(NEW.NAMES))
    {
      if(trace) { message("Fitting RSF model ",NEW.NAMES[i]) }
      NEW[[NEW.NAMES[i]]] <- rsf.fit(data,UD,R=R,formula=FORM[[i]],trace=max(trace-1,0),...)
    }
    M <- c(M,NEW)

    # best RSF parameters so far (going down)
    NEW <- M[NAMES]
    ICS <- sapply(NEW,function(m){m[[IC]]})
    i <- which.min(ICS)
    beta <- NEW[[i]]$beta
    ON <- TERMS %in% names(beta)
    # copy for initial guess in fitting
    # UD@CTMM$beta <- beta
    UD@CTMM <- NEW[[i]]
  }

  ICS <- sapply(M,function(m){m[[IC]]})
  i <- sort(ICS,index.return=TRUE)$ix
  M <- M[i]

  if(!verbose) { M <- M[[1]] }

  return(M)
}
