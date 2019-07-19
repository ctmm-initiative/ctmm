# cross-validation model-selection routines

# leave-one-out cross validated likelihood for model selection
LOOCV <- function(data,CTMM,cores=1,method=CTMM$method,...)
{
  cores <- resolveCores(cores,fast=FALSE)
  axes <- CTMM$axes
  AXES <- length(axes)
  n <- nrow(data)

  cv.like <- function(i)
  {
    # error matrix on the ith location
    ERROR1 <- get.error(data[i,],CTMM,DIM=AXES)
    dim(ERROR1) <- c(AXES,AXES)

    # predict ith location using the same class model
    CTMM <- ctmm.fit(data[-i,],CTMM,method=method,COV=FALSE,...)
    PRED <- predict(data[-i,],CTMM,t=data$t[i])

    # uncertainty matrix on the prediction
    CTMM$error <- TRUE
    ERROR2 <- get.error(PRED,CTMM,DIM=AXES)
    dim(ERROR2) <- c(AXES,AXES)

    # total uncertainty
    ERROR <- ERROR1 + ERROR2

    # standardized distance^2
    D <- get.telemetry(PRED,axes) - get.telemetry(data[i,],axes)
    D <- c(D)
    D <- D %*% PDsolve(ERROR) %*% D
    D <- c(D)

    return( -log(det(2*pi*ERROR))/2 - D/2 )
  }

  LIKE <- plapply(1:n,cv.like,cores=cores,fast=FALSE)
  LIKE <- unlist(LIKE)
  LIKE <- sum(LIKE)

  return(-2*LIKE)
}
