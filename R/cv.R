# cross-validation model-selection routines
###########################################

# conditional log-likelihood of !IN | IN
cv.like <- function(data,CTMM,IN,method,...)
{
  axes <- CTMM$axes
  CTMM <- ctmm.fit(data[IN,],CTMM,method=method,COV=FALSE,...)

  ## detrend mean, so it isn't modified in ctmm.loglike
  drift <- drift.mean(CTMM,data$t) %*% CTMM$mu
  data[,axes] <- get.telemetry(data,axes=axes) - drift
  CTMM$mean <- "zero" # set to mean-zero process

  # conditional log-likelihood of OUT=!IN | IN
  LIKE <- ctmm.loglike(data,CTMM,REML=FALSE,profile=FALSE) - ctmm.loglike(data[IN,],CTMM,REML=FALSE,profile=FALSE)

  return(LIKE)
}


# leave-one-out cross validated likelihood for model selection
LOOCV <- function(data,CTMM,cores=1,method=CTMM$method,...)
{
  cores <- resolveCores(cores,fast=FALSE)
  n <- nrow(data)

  cv.like.i <- function(i) { cv.like(data,CTMM,IN=-i,method=method,...) }

  LIKE <- plapply(1:n,cv.like.i,cores=cores,fast=FALSE)
  LIKE <- unlist(LIKE)
  LIKE <- sum(LIKE)

  return(-2*LIKE)
}


# half-sample cross validation
HSCV <- function(data,CTMM,cores=1,method=CTMM$method,...)
{
  cores <- resolveCores(cores,fast=FALSE)

  t <- data$t
  t.mid <- t[1] + (last(t)-t[1])/2

  IN <- which(t<=t.mid)
  OUT <- which(t>t.mid)
  LIKE <- cv.like(data,CTMM,IN,method=method,...) + cv.like(data,CTMM,OUT,method=method,...)

  if(t.mid %in% t) # average splits
  {
    IN <- which(t<t.mid)
    OUT <- which(t>=t.mid)
    LIKE <- LIKE + cv.like(data,CTMM,IN,method=method,...) + cv.like(data,CTMM,OUT,method=method,...)
    LIKE <- LIKE/2
  }

  return(-2*LIKE)
}
