# relative encounter rates (trajectory based)
encounter.ecdf <- function(data,UD,level=0.95,debias=TRUE,res.time=1,r=NULL,...)
{
  if(class(UD[[1]])[1]=="UD")
  { CTMM <- lapply(UD,function(ud){ud@CTMM}) }
  else
  { CTMM <- UD }

  t1 <- c(data[[1]]$t[1],data[[2]]$t[1])
  t2 <- c(last(data[[1]]$t),last(data[[2]]$t))
  OVER <- ( min(t2)-max(t1) ) / ( max(t2)-min(t1) )

  if(OVER<=0) { stop("No overlapping times.") }

  DIFF <- difference(data,CTMM,uniform=TRUE,res.time=1,...)
  # prediction information
  R2 <- DIFF$x^2+DIFF$y^2
  VAR <- 2*DIFF$VAR.xy
  # worst/null prediction
  VAR0 <- sum((CTMM[[1]]$mu-CTMM[[2]]$mu)^2) + var.covm(CTMM[[1]]$sigma) + var.covm(CTMM[[2]]$sigma)
  # added information
  DOF <- VAR0/VAR - 1
  DOF <- clamp(DOF,0,Inf)
  # uncertainty of added information
  VAR <- VAR0/DOF
  # natural weights
  w <- VAR0/(VAR0+VAR)
  w <- w/sum(w)

  if(is.null(r)) # default grid (roughly 1% increments)
  {
    dr <- sqrt(stats::median(R2))/50
    r <- (1:100)*dr
  }
  # cumulative probability
  P <- array(0,length(r))

  # smoothed estimate
  R <- sqrt(R2)
  for(i in 1:length(R))
  {
    j <- r>=R[i]
    P[j] <- P[j] + w[i]
  }

  if(debias)
  {
    # varying estimate
    P2 <- array(0,length(r))
    VAR <- 2*DIFF$VAR.xy
    DOF <- 2*VAR0/VAR
    r2 <- r^2

    for(i in 1:length(R2))
    {
      X2 <- r2/R2[i] # reduced X^2
      X2 <- X2 * DOF[i] # X^2
      P2 <- P2 + w[i] * stats::pchisq(X2,DOF[i])
    }

    BIAS <- P2/P
    BIAS <- nant(BIAS,1)
    P <- P/BIAS
  }

  # PDF-based calculation
  DOF <- encounter(data,CTMM,debias=FALSE,method="PDF")$DOF[1,2]
  DOF <- DOF * OVER # correction for partial temporal overlap

  R <- data.frame(r=r,P=P,P=P,P=P)
  names(R) <- c('r',NAMES.CI)
  alpha <- 1-level
  # binomial CIs
  for(i in 1:length(P)) { R[i,NAMES.CI] <-  beta.ci(P[i],2*P[i]^2/DOF,level=level) }

  # class(R) <- "ECDF"

  return(R)
}


# relative encounter rates (UD based)
encounter <- function(data,UD,method="ECDF",debias=TRUE,level=0.95,r=NULL,res.time=1,normalize=FALSE,self=TRUE,...)
{
  units <- FALSE
  method <- match.arg(method,c("PDF","ECDF"))

  if(method=="ECDF")
  { return(encounter.ecdf(data,UD,debias=debias,r=r,res.time=res.time,...)) }

  if(class(data[[1]])[1]=="UD")
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
