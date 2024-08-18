# relative encounter rates (trajectory based)
encounter.ecdf <- function(data,CTMM,debias=TRUE,res.time=1,...)
{
  DIFF <- difference(data,CTMM,uniform=TRUE,res.time=1,...)
  # prediction information
  R2 <- DIFF$x^2+DIFF$y^2
  VAR <- 2*DIFF$VAR.xy
  # worst/null prediction
  VAR.0 <- sum((CTMM[[1]]$mu-CTMM[[2]]$mu)^2) + var.covm(CTMM[[1]]$sigma) + var.covm(CTMM[[2]]$sigma)
  # added information
  DOF <- VAR.0/VAR - 1
  DOF <- clamp(DOF,0,Inf)
  # uncertainty of added information
  VAR <- VAR.0/DOF
  # natural weights
  w <- VAR.0/(VAR.0+VAR)
  # chi^2 DOF

}


# relative encounter rates (UD based)
encounter <- function(data,UD,debias=FALSE,level=0.95,method="ECDF",normalize=FALSE,self=TRUE,...)
{
  units <- FALSE
  method <- match.arg(method,c("DF","ECDF"))

  if(method=="ECDF")
  {
    CTMM <- lapply(UD,function(ud){ud@CTMM})
    return(encounter.ecdf(data,CTMM,debias=debias,...))
  }

  if(class(data[[1]])=="UD")
  {
    UD <- data
    rm(data)
  }
  object <- UD
  rm(UD)

  R <- overlap(object,debias=debias,level=level,method="Rate",...)
  R$CI <- nant(R$CI,0)
  R$DOF <- nant(R$DOF,0)

  if(normalize)
  {
    # calculate mean self-rate
    s <- diag(R$CI[,,'est'])
    dof <- diag(R$DOF)
    s <- meta.chisq(s,dof)$CI["mean","est"]

    R$CI <- R$CI/s
  }

  # fix diagonals # self encounter rate
  if(self)
  {
    diag(R$CI[,,1]) <- diag(R$CI[,,2]) <- diag(R$CI[,,3]) <- 1
    diag(R$DOF) <- Inf
  }

  R$CI <- R$CI * pi # standardized encounter probability (1-meter)

  return(R)
}
