suitability <- function(data=NULL,CTMM=NULL,R=list(),level=0.95,grid=NULL,log=FALSE,...)
{
  if(is.null(grid)) { stop("Please provide a grid argument, such as a UD or raster object.") }

  if(class(data)[1]=="ctmm")
  {
    TEMP <- CTMM
    CTMM <- data
    data <- TEMP
    rm(TEMP)
  }

  grid <- format_grid(grid=grid,axes=CTMM$axes)
  x <- grid$r$x
  y <- grid$r$y

  R <- expand.factors(R,CTMM$formula,fixed=TRUE)$R

  proj <- CTMM@info$projection
  # calculate RASTERs on spatial grid
  R <- lapply(R,function(S){R.grid(r=grid$r,proj=proj,S)})

  R <- R.suit(R,CTMM,data=data,log=log,VAR=TRUE)
  VAR <- R$VAR
  R <- R$R

  alpha <- 1-level
  z <- stats::qnorm(1-alpha/2)

  if(!log)
  {
    # details in lognorm.ci()
    VAR <- log(1 + VAR/R^2)
    log.R <- log(R) - VAR/2
  }

  STD <- z * sqrt(VAR)
  rm(VAR)

  R <- array(R,c(dim(R),3))
  if(!log)
  {
    R[,,1] <- exp(log.R - STD)
    R[,,3] <- exp(log.R + STD)
    rm(log.R)
  }
  else
  {
    R[,,1] <- R[,,2] - STD
    R[,,3] <- R[,,2] + STD
  }

  proj <- sp::CRS(proj)
  R <- lapply(1:3,function(i){raster::raster(list(x=x,y=y,z=R[,,i]),crs=proj)})
  names(R) <- NAMES.CI
  R <- raster::brick(R[[1]],R[[2]],R[[3]])
  names(R) <- NAMES.CI

  return(R)
}


# evaluate habitat suitability raster(s)
# data.frame can only have one row used in this function
R.suit <- function(R,CTMM,data=NULL,log=FALSE,VAR=FALSE)
{
  DIM <- dim(R[[1]])
  beta <- CTMM$beta
  formula <- CTMM$formula

  OFFSET <- get.offset(CTMM$formula)

  RVARS <- names(R)
  VARS <- all.vars(formula)
  DVARS <- VARS[ VARS %nin% RVARS ]

  R <- lapply(R,c)
  R <- data.frame(R)

  if(length(OFFSET))
  {
    OFFSET <- stats::model.frame(formula,R)
    OFFSET <- stats::model.offset(OFFSET)
  }

  for(D in DVARS) { R[[D]] <- as.numeric(data[[D]])[1] } # model.matrix will rename otherwise; only use first data row

  # working subset that model.matrix will return (NA will be skipped)
  SUB <- apply(R,1,function(x){!any(is.na(x))})
  # model.matrix - NA have been skipped
  R <- stats::model.matrix(formula,data=R)

  TERMS <- colnames(R)
  TERMS <- TERMS[TERMS!="(Intercept)"]
  R <- R[,TERMS,drop=FALSE]

  if(VAR)
  {
    COV <- CTMM$COV[TERMS,TERMS,drop=FALSE]
    COV <- sapply(1:nrow(R),function(i){R[i,] %*% COV %*% R[i,]})

    if("POV" %in% names(CTMM))
    {
      POV <- CTMM$POV[TERMS,TERMS]
      POV <- sapply(1:nrow(R),function(i){R[i,] %*% POV %*% R[i,]})

      CTERMS <- cross.terms(TERMS)
      COV.POV <- CTMM$COV.POV[CTERMS,CTERMS]
      r <- quad2lin(R,diag=TRUE)
      COV.POV <- sapply(1:nrow(R),function(i){r[i,] %*% COV.POV %*% r[i,]})
    }
    else
    {
      POV <- 0
      COV.POV <- 0
    }
  }
  R <- c(R %*% beta[TERMS])

  if(!log)
  {
    if(VAR)
    {
      R <- exp_log(R,VAR.est=COV,VAR=POV,VAR.VAR=COV.POV)
      COV <- R$VAR
      R <- R$mu
    }
    else
    { R <- exp(R) }

    if(length(OFFSET))
    {
      R <- OFFSET * R
      if(VAR) { COV <- OFFSET * COV }
    }
  }
  else if(length(OFFSET))
  {
    R <- R + log(OFFSET)
    if(VAR) { COV <- OFFSET * COV }
  }

  # NA default
  if(log)
  { FULL <- array(-Inf,DIM) }
  else
  { FULL <- array(0,DIM) }

  # copy over non-NA values
  FULL[SUB] <- R
  R <- FULL

  if(VAR)
  {
    # NA default
    FULL <- array(0,DIM)

    FULL[SUB] <- COV
    COV <- FULL

    return(list(R=R,VAR=COV))
  }
  else
  { return(R) }
}
