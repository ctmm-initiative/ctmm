AICc.list <- function(object,...)
{
  LIST <- sapply(object,function(M) { AIC <- AICc(M,...) ; c(AIC$mean,AIC$error) } )
  LIST <- t(LIST)

  colnames(LIST) <- c("mean","error")

  ORDER <- sort(LIST[,1],method="quick",index.return=TRUE)$ix
  LIST <- LIST[ORDER,]

  return(LIST)
}

# double-bootstrap AICc estimate for one model
AICc.ctmm <- function(object,data,n=100,fast=FALSE,...)
{
  CTMM <- object
  # toss location information for simulations
  data[,CTMM$axes] <- NULL

  LIKE <- numeric(n)
  for(i in 1:n)
  {
    DATA <- simulate(CTMM,data=data)

    FIT <- emulate(CTMM,data=data,fast=fast,...)

    LIKE[i] <- ctmm.loglike(DATA,FIT,profile=FALSE)
  }

  # now these are AICc targets
  LIKE <- -2*LIKE

  MEAN <- mean(LIKE)
  SE <- sqrt(stats::var(LIKE)/n)

  RESULT <- list(mean=MEAN,error=SE)
  return(RESULT)
}
