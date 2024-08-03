# match COV.mu structure to flattened mean c(mu) # [axes*m]
flatten.cov.mu <- function(COV.mu)
{
  if(length(dim(COV.mu))==4) # [axes,m,m,axes]
  {
    COV.mu <- aperm(COV.mu,c(2,1,3,4)) # [m,axes,m,axes]
    DIM <- dim(COV.mu)
    DIM <- prod(DIM[1:2])
    dim(COV.mu) <- c(DIM,DIM) # [axes*m,axes*m]
  }

  return(COV.mu)
}

##########
mean.mu <- function(x,debias=TRUE,weights=NULL,trace=FALSE,IC="AICc",...)
{
  if(is.null(weights))
  { weights <- rep(1,length(x)) }

  axes <- x[[1]]$axes
  AXES <- length(axes)
  N <- length(x)

  mean <- sapply(x,function(y){y$mean})
  mean <- unique(mean)

  # dimensions of mean functions
  m <- sapply(x,function(y){nrow(y$mu)})
  M <- max(m)

  # Gaussian-Gaussian in all cases
  MU <- array(0,c(N,AXES*M))
  SIGMA <- array(0,c(N,AXES*M,AXES*M))
  INF <- array(TRUE,c(N,AXES*M)) # missing data

  if(M==1)
  { CNAMES <- axes }
  else
  {
    NAMES <- which(m==M)[1] # index of most detailed mean
    NAMES <- rownames(x[[NAMES]]$mu)
    CNAMES <- c( outer(NAMES,axes,function(s1,s2){paste(s1,s2,sep="-")}) )
  }
  colnames(MU) <- CNAMES
  dimnames(SIGMA) <- list(NULL,CNAMES,CNAMES)

  for(i in 1:N)
  {
    # fill in for missing variances
    diag(SIGMA[i,,]) <- Inf

    # existing indices
    SUB <- logical(M*AXES)
    for(j in 1:AXES) { SUB[(j-1)*M+1:m[i]] <- TRUE }
    SUB <- which(SUB)

    MU[i,SUB] <- c(x[[i]]$mu)
    SIGMA[i,SUB,SUB] <- flatten.cov.mu(x[[i]]$COV.mu)
    INF[i,SUB] <- FALSE
  }
  DOF <- colSums(!INF) # amount of data per mode

  # list of candidate models
  MM <- list()

  # no variance in mean location
  S <- "Dirac-\u03B4(\u03BC)"
  if(trace) { message("Fitting location-mean model ",S) }
  MM[[S]] <- meta.normal(MU,SIGMA,VARS=FALSE,isotropic=TRUE,debias=debias,weights=weights,WARN=FALSE,...)
  GUESS <- MM[[S]]

  if(M==1) # stationary
  {
    DOF <- sum(DOF)
    if(DOF > AXES*M) # can support a variance
    {
      # symmetric mean distribution
      S <- "isotropic-\u03BC"
      if(trace) { message("Fitting location-mean model ",S) }
      MM[[S]] <- meta.normal(MU,SIGMA,isotropic=TRUE,debias=debias,weights=weights,GUESS=GUESS,WARN=FALSE,...)
      GUESS <- MM[[S]]

      # can support a covariance
      if(DOF > AXES*M + (AXES*M*(AXES*M+1))/2)
      {
        S <- "anisotropic-\u03BC"
        if(trace) { message("Fitting location-mean model ",S) }
        MM[[S]] <- meta.normal(MU,SIGMA,debias=debias,weights=weights,GUESS=GUESS,WARN=FALSE,...)
      }
    } # end if sum(DOF) > AXES*M
  } # end if M==1 stationary
  else # non-stationary
  {
    VAR <- DOF>1
    if(any(VAR))
    {
      # truncated block diagonal
      VARS <- diag(VAR)
      S <- "VAR[\u03BC]"
      if(trace) { message("Fitting location-mean model ",S) }
      MM[[S]] <- meta.normal(MU,SIGMA,VARS=VARS,debias=debias,weights=weights,GUESS=GUESS,WARN=FALSE,...)
      GUESS <- MM[[S]]

      dofs <- sort(unique(DOF),decreasing=TRUE)
      dofs <- dofs[dofs>1]
      # build up correlations
      for(d in dofs)
      {
        SUB <- DOF>=d
        VARS[SUB,SUB] <- TRUE
        if(sum(VARS) <= sum(DOF[diag(VARS)]))
        {
          # truncated unstructured
          VARS <- diag(VAR)
          S <- paste0("COV[\u03BC] ",sum(SUB),"/",length(SUB))
          if(trace) { message("Fitting location-mean model ",S) }
          MM[[S]] <- meta.normal(MU,SIGMA,VARS=VARS,debias=debias,weights=weights,GUESS=GUESS,WARN=FALSE,...)
          GUESS <- MM[[S]]
        }
        else # ran out of data
        { break }
      } # end correlation build up
    }
  } # end non-stationary

  ICS <- sapply(MM,function(m){m[[IC]]})
  names(ICS) <- names(MM)
  i <- which.min(ICS)
  i <- names(MM)[i]
  MM <- MM[[i]]
  MM$name <- i

  # report model selection
  if(trace)
  {
    i <- sort(ICS,index.return=TRUE)$ix
    ICS <- ICS[i] # sorted
    ICS <- ICS - ICS[1] # start at zero
    ICS <- data.frame(ICS)
    colnames(ICS) <- paste0("\u0394",IC)
    message("* Model selection for location-mean \u03BC distribution.")
    print(ICS)
  }

  # R$mu # population mean of mean locations
  # R$COV.mu # uncertainty in mean of means estimate
  names(MM)[ which(names(MM)=="sigma") ] <- "POV.mu" # dispersion of means
  names(MM)[ which(names(MM)=="COV.sigma") ] <- "COV.POV.mu" # uncertainty in dispersion of means

  if(M>1)
  {
    MM$mu <- array(MM$mu,c(M,AXES))
    dimnames(MM$mu) <- list(NAMES,axes)

    dim(MM$COV.mu) <- c(M,AXES,M,AXES)
    MM$COV.mu <- aperm(MM$COV.mu,c(2,1,3,4))
    dimnames(MM$COV.mu) <- list(axes,NAMES,NAMES,axes)

    dim(MM$POV.mu) <- c(M,AXES,M,AXES)
    MM$POV.mu <- aperm(MM$POV.mu,c(2,1,3,4))
    dimnames(MM$POV.mu) <- list(axes,NAMES,NAMES,axes)

    CNAMES <- outer(NAMES,axes,function(s1,s2){paste(s1,s2,sep="-")})
  }
  else
  {
    MM$mu <- rbind(MM$mu)
    colnames(MM$mu) <- axes

    dimnames(MM$COV.mu) <- list(axes,axes)
    dimnames(MM$POV.mu) <- list(axes,axes)

    CNAMES <- outer(axes,axes,function(s1,s2){paste(s1,s2,sep="-")})
  }
  CNAMES <- c(CNAMES)
  CNAMES <- outer(CNAMES,CNAMES,function(s1,s2){paste(s1,s2,sep="-")})
  CNAMES <- CNAMES[upper.tri(CNAMES,diag=TRUE)]
  dimnames(MM$COV.POV.mu) <- list(CNAMES,CNAMES)

  MM$mean <- mean

  return(MM)
}
