overlap <- function(object1,object2,...) UseMethod("overlap") #S3 generic

overlap.ctmm <- function(object1,object2,level=0.95,...)
{
  CTMM1 <- object1
  CTMM2 <- object2
  
  # ML parameters
  par <- NULL
  
  par <- c(par,CTMM1$mu[1,])
  par <- c(par,CTMM1$sigma@par)

  par <- c(par,CTMM2$mu[1,])
  par <- c(par,CTMM2$sigma@par)
  
  D <- function(p)
  {
    CTMM1$mu[1,] <- p[1:2]
    CTMM1$sigma <- covm(p[3:5])
    
    CTMM2$mu[1,] <- p[6:7]
    CTMM2$sigma <- covm(p[8:10])
    
    return(BhattacharyyaD(CTMM1,CTMM2))
  }
  
  # propagate uncertainty
  grad <- numDeriv::grad(D,par)
  VAR <- 0
  
  VAR <- VAR + grad[1:2] %*% CTMM1$COV.mu %*% grad[1:2]
  if(CTMM1$isotropic)
  { VAR <- VAR + grad[3] * CTMM1$COV[1,1] * grad[3] }
  else
  { VAR <- VAR + grad[3:5] %*% CTMM1$COV[1:3,1:3] %*% grad[3:5] }
    
  VAR <- VAR + grad[6:7] %*% CTMM1$COV.mu %*% grad[6:7]
  if(CTMM1$isotropic)
  { VAR <- VAR + grad[8] * CTMM1$COV[1,1] * grad[8] }
  else
  { VAR <- VAR + grad[8:10] %*% CTMM1$COV[1:3,1:3] %*% grad[8:10] }
  
  # this quantity is roughly chi-square
  MLE <- BhattacharyyaD(CTMM1,CTMM2)
  CI <- chisq.ci(MLE,COV=VAR,alpha=1-level)
  
  # transform from (square) distance to overlap measure
  CI <- exp(-rev(CI))
  names(CI) <- c("low","ML","high")
  
  if(level) { return(CI) }
  else { return(list(MLE=MLE,VAR=VAR)) }
}

# square distance between stationary Gaussian distributions
BhattacharyyaD <- function(CTMM1,CTMM2)
{
  sigma <- (CTMM1$sigma + CTMM2$sigma)/2
  mu <- CTMM1$mu[1,] - CTMM2$mu[1,]
  
  D <- as.numeric(mu %*% solve(sigma) %*% mu)/8 + log(det(sigma)/sqrt(det(CTMM1$sigma)*det(CTMM2$sigma)))/2

  return(D)
}

#overlap density function
overlay <- function(data1,data2,CTMM1,CTMM2,debias=TRUE,error=0.001,res=10,grid=NULL)
{
  HP1 <- akde.bandwidth(data=data1,CTMM=CTMM1,verbose=TRUE)
  HP2 <- akde.bandwidth(data=data2,CTMM=CTMM2,verbose=TRUE)
  
  # concatenate bandwidth arrays
  H1 <- prepare.H(data1,HP1$H)
  H2 <- prepare.H(data2,HP2$H)
  
  n1 <- length(data1$x)
  n2 <- length(data2$x)
  
  H1 <- array(H1,c(n1,2*2))
  H2 <- array(H2,c(n2,2*2))

  H <- rbind(H1,H2)
  H <- array(H,c(n1+n2,2,2))
  H1 <- array(H1,c(n1,2,2))
  H2 <- array(H2,c(n2,2,2))
  
  # concatenate data
  data <- rbind(data1,data2)
  
  # highest resolution
  dr <- sapply( 1:2 , function(i) { sqrt(min(HP1$H[i,i],HP2$H[i,i])) } )/res
  
  # construct joint grid
  grid <- kde.grid(data,H,alpha=error,dr=dr)
  DX <- grid$DX
  DY <- grid$DY
  
  # calcualte individual density esitmates
  grid$DX <- utils::head(DX,n=n1)
  grid$DY <- utils::head(DY,n=n1)
  if(debias) { debias <- HP1$bias }
  UD1 <- kde(data1,H=H1,grid=grid,alpha=error,bias=debias)

  grid$DX <- utils::tail(DX,n=n2)
  grid$DY <- utils::tail(DY,n=n2)
  if(debias) { debias <- HP2$bias }
  UD2 <- kde(data2,H=H2,grid=grid,alpha=error,bias=debias)

  PDF <- sqrt(UD1$PDF * UD2$PDF)
  dA <- prod(dr)

  # overlap point estimate
  OVER <- sum(PDF)*dA
  # calculate Gaussian overlap distance^2 variance
  CI <- overlap(CTMM1,CTMM2,level=FALSE)

  # normalize the denisty function for plotting
  PDF <- PDF/OVER
  
  # create CDF for plotting
  CDF <- pdf2cdf(PDF*dA)
  
  # choose largest bandwidth for grid.....
  H <- diag(c(max(HP1$H[1,1],HP2$H[1,1]),max(HP1$H[2,2],HP2$H[2,2])),2)
  
  UD <- list(PDF=PDF,CDF=CDF,x=grid$x,y=grid$y,dA=dA,H=H)
  UD$D <- -log(OVER)
  UD$COV.D <- CI$VAR
  
  UD <- new.UD(UD,info=mean.info(list(data1,data2)))
  return(UD)
}

# overlap quantity
overlap.telemetry <- function(object1,object2,CTMM1,CTMM2,level=0.95,...)
{
  data1 <- object1
  data2 <- object2
  
  OVER = overlay(data1,data2,CTMM1,CTMM2,...)
  
  # calculate new distance^2 with KDE point estimate
  CI <- chisq.ci(OVER$D,COV=OVER$COV.D,alpha=1-level)
  
  # transform from (square) distance to overlap measure
  CI <- exp(-rev(CI))
  names(CI) <- c("low","ML","high")
  
  return(CI)
}