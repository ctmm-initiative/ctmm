rsf.select <- function(data,UD,R=list(),formula=NULL,verbose=FALSE,IC="AICc",trace=TRUE,...)
{
  if(length(R))
  {
    RVARS <- names(R)

    if(is.null(formula))
    { TERMS <- RVARS }
    else
    {
      # this doesn't work with poly(), etc.
      # TERMS <- attr(stats::terms(formula),"term.labels")
      # dummy data for model parameter names
      DATA <- data.frame(data)[1:2,]
      for(r in names(R)) { DATA[[r]] <- 0 }
      TERMS <- colnames(stats::model.matrix(formula,data=DATA))
      TERMS <- TERMS[TERMS!="(Intercept)"]
    }
  }
  DIM <- length(TERMS)

  M <- list()
  ON <- array(FALSE,DIM)
  names(ON) <- TERMS

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
    NAMES <- sapply(1:length(try),function(i){paste("~",paste(TERMS[TRYS[i,]],collapse="+"))})
    # don't try the same model twice
    NAMES <- NAMES[NAMES %nin% names(M)]

    # fit all
    NEW <- list()
    for(i in 1%:%length(NAMES))
    {
      formula <- stats::as.formula(NAMES[i])
      if(trace) { message("Fitting RSF model ",NAMES[i]) }
      NEW[[NAMES[i]]] <- rsf.fit(data,UD,R=R,formula=formula,trace=trace-1,...)
    }
    M <- c(M,NEW)

    # best RSF parameters so far (going up)
    ICS <- sapply(NEW,function(m){m[[IC]]})
    i <- which.min(ICS)
    beta <- NEW[[i]]$beta
    ON <- TERMS %in% names(beta)
    # copy for initial guess in fitting
    UD@CTMM$beta <- beta
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
    NAMES <- sapply(1:length(try),function(i){paste("~",paste(TERMS[TRYS[i,]],collapse="+"))})
    NAMES[NAMES=="~ "] <- "~ 0" # make a valid formula string
    # don't try the same model twice
    NEW.NAMES <- NAMES[NAMES %nin% names(M)]

    # fit all
    NEW <- list()
    for(i in 1%:%length(NEW.NAMES))
    {
      if(NEW.NAMES[i]=="~ 0") # no terms left
      { formula <- ~ 0 }
      else
      { formula <- stats::as.formula(NEW.NAMES[i]) }
      if(trace) { message("Fitting RSF model ",NEW.NAMES[i]) }
      NEW[[NEW.NAMES[i]]] <- rsf.fit(data,UD,R=R,formula=formula,trace=trace-1,...)
    }
    M <- c(M,NEW)

    # best RSF parameters so far (going down)
    NEW <- M[NAMES]
    ICS <- sapply(NEW,function(m){m[[IC]]})
    i <- which.min(ICS)
    beta <- NEW[[i]]$beta
    ON <- TERMS %in% names(beta)
    # copy for initial guess in fitting
    UD@CTMM$beta <- beta
  }

  ICS <- sapply(M,function(m){m[[IC]]})
  i <- sort(ICS,index.return=TRUE)$ix
  M <- M[i]

  if(!verbose) { M <- M[[1]] }

  return(M)
}
