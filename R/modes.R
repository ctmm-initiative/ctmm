# statistical mode
modes.numeric <- function(x,method="exact")
{
  x <- sort(x)
  unique.x <- unique(x)

  if(length(unique.x)<length(x))
  {
    freq <- tabulate(match(x, unique.x))
    MAX <- max(freq)
    IND <- which(freq==MAX)

    # well defined sample mode
    if(length(IND)==1) { return(x[IND]) }
    # else minimize log-distance
  }
  else
  { IND <- 1:length(x) }

  # log distance
  D <- log(abs(outer(x,x,'-')))
  # no self distance
  diag(D) <- 0
  # no ties
  if(length(IND)<length(x))
  {
    D[IND,IND] <- 0
    D[!IND,!IND] <- Inf
  }
  # total distance to other locations
  TD <- rowSums(D)
  # mode index
  IND <- which.min(TD)

  x[IND]
}


# assuming stationary model
modes.ctmm <- function(object,...)
{
  mu <- object$mu[1,]
  COV.mu <- object$COV.mu
  if(length(dim(COV.mu))==4) { COV.mu <- COV.mu[,1,1,] }

  axes <- object$axes
  data <- data.frame(t(mu))
  names(data) <- axes

  for(i in 1:length(axes))
  {
    for(j in i:length(axes))
    {
      NAME <- paste("COV.",axes[i],".",axes[j],sep="") # consistent with imported ARGOS error ellipse notation
      data[,NAME] <- COV.mu[i,j]
    }
  }

  data <- new.telemetry(data,info=object@info)

  return(data)
}


# calculate second-order, 2D gradient and Hessian from 3x3 grid
DiffGrid <- function(SUB,dx,dy,dx2=dx^2,dy2=dy^2,dr2=dx2+dy2,dxdydr2=dx*dy/dr2)
{
  # calculate gradient
  GRAD <- numeric(2)
  GRAD[1] <- mean( diff(SUB[2,1:2]) , diff(SUB[2,2:3]) )/dx
  GRAD[2] <- mean( diff(SUB[1:2,2]) , diff(SUB[2:3,2]) )/dy

  # calculate Hessian
  HESS <- diag(2)
  HESS[1,1] <- diff(diff(SUB[2,]))/dx2 # (+x,+x)
  HESS[2,2] <- diff(diff(SUB[,2]))/dy2 # (+y,+y)
  HESS[1,2] <- diff(diff(c(SUB[1,1],SUB[2,2],SUB[3,3])))/dr2 # (+x+y,+x+y)
  HESS[2,1] <- diff(diff(c(SUB[1,3],SUB[2,2],SUB[3,1])))/dr2 # (+x-y,+x-y)
  HESS[1,2] <- HESS[2,1] <- dxdydr2 * (HESS[1,2]-HESS[2,1]) # (+x,+y)

  return(list(GRAD=GRAD,HESS=HESS))
}


#
modes.UD <- function(object,...)
{
  CTMM <- object@CTMM
  axes <- CTMM$axes
  COV <- CTMM$sigma
  COV.mu <- CTMM$COV.mu
  if(length(dim(COV.mu))==4) { COV.mu <- COV.mu[,1,1,] }
  COV.mu <- covm(COV.mu,axes=axes,isotropic=CTMM$isotropic)
  # SCALE %*% COV %*% t(SCALE) -> COV.mu
  # SCALE == DOF^(-1/2)
  SCALE <- mpow.covm(COV.mu,1/2) %*% mpow.covm(COV,-1/2)

  PDF <- object$PDF
  CDF <- object$CDF
  dx <- object$dr['x']
  dy <- object$dr['y']
  dx2 <- dx^2
  dy2 <- dy^2
  dr2 <- dx^2+dy^2
  dxdydr2 <- dx*dy/dr2

  data <- data.frame(t(rep(0,length(axes))))
  names(data) <- axes
  m <- 0

  for(i in 2:(nrow(PDF)-1))
  {
    for(j in 2:(ncol(PDF)-1))
    {
      SUB <- PDF[i+(-1):1,j+(-1):1]
      TEST <- all( SUB[5] > SUB[-5] ) # SUB[5] == P[i,j]
      if(TEST) # local maximum
      {
        m <- m + 1

        SUB <- log(SUB)
        STUFF <- DiffGrid(SUB,dx,dy,dx2,dy2,dr2,dxdydr2)
        GRAD <- STUFF$GRAD
        HESS <- STUFF$HESS
        COV <- -PDsolve(HESS)

        # locate mode
        dr <- c(COV %*% GRAD)
        # in case of bad Hessian
        dr[1] <- clamp(dr[1],-dx/2,dx/2)
        dr[2] <- clamp(dr[2],-dy/2,dy/2)
        r <- c(object$r$x[i],object$r$y[j]) + dr
        # store mode
        data[m,axes] <- r

        # GRF mode uncertainty ellipse
        COV <- SCALE %*% COV %*% t(SCALE)

        # CDF/N heuristic correction
        SUB <- CDF[i+(-1):1,j+(-1):1]
        STUFF <- DiffGrid(SUB,dx,dy,dx2,dy2,dr2,dxdydr2)
        GRAD <- STUFF$GRAD
        HESS <- STUFF$HESS
        cdf <- CDF[i,j] + c(GRAD %*% dr) + c(dr %*% HESS %*% dr)/2
        cdf <- clamp(cdf,0,1)
        COV <- COV/(1-cdf) # DOES THIS MAKE SENSE?

        # store uncertainty ellipse
        for(I in 1:length(axes))
        {
          for(J in I:length(axes))
          {
            NAME <- paste("COV.",axes[I],".",axes[J],sep="") # consistent with imported ARGOS error ellipse notation
            data[m,NAME] <- COV[I,J]
          }
        } # end uncertainty ellipse
      } # end if local minimum
    } # end col loop
  } # end row loop

  data <- new.telemetry(data,info=CTMM@info)
  return(data)
}


# rewrite to work from CDF, via least steep descent

#
ridges.UD <- function(object,...)
{
  CTMM <- object@CTMM
  axes <- CTMM$axes
  COV <- CTMM$sigma
  COV.mu <- CTMM$COV.mu
  if(length(dim(COV.mu))==4) { COV.mu <- COV.mu[,1,1,] }
  COV.mu <- covm(COV.mu,axes=axes,isotropic=CTMM$isotropic)
  # SCALE %*% COV %*% t(SCALE) -> COV.mu
  # SCALE == DOF^(-1/2)
  SHRINK <- mpow.covm(COV.mu,1/2) %*% mpow.covm(COV,-1/2)
  GROW <- mpow.covm(COV.mu,-1/2) %*% mpow.covm(COV,1/2)

  PDF <- object$PDF # log first in ridge code? NO!
  CDF <- object$CDF
  dx <- object$dr['x']
  dy <- object$dr['y']
  dx2 <- dx^2
  dy2 <- dy^2
  dr2 <- dx^2+dy^2
  dxdydr2 <- dx*dy/dr2

  data <- data.frame(t(rep(0,length(axes))))
  names(data) <- axes
  m <- 0

  for(i in 2:(nrow(PDF)-1))
  {
    for(j in 2:(ncol(PDF)-1))
    {
      SUB <- PDF[i+(-1):1,j+(-1):1]
      if(all(SUB>0)) # non-zero density (under logarithm)
      {
        DIFF <- DiffGrid(SUB,dx,dy,dx2,dy2,dr2,dxdydr2)
        GRAD <- DIFF$GRAD # 1/p D p # under log
        HESS <- DIFF$HESS # 1/p DD p - 1/p^2 Dp Dp # under log
        #HESS <- HESS + (GRAD %o% GRAD) # logarithm correction
        EIGEN <- eigen(HESS)
        if(last(EIGEN$values) < -.Machine$double.eps) # some negative curvature
        {
          VAR <- -1/EIGEN$values[2] # -1/ridge curvature > 0
          DIR <- EIGEN$vectors[,2] # tangent to ridge line
          dr <- VAR*c(DIR%*%GRAD) * DIR # move from center of cell up to the ridge line
          # this is the closest ridge-line point to the center of the cell

          # local ridge-line equation
          # r(t) = r + dr + t*EIGEN$vectors[,1]
          # grad(t) = (GRAD + t*HESS) * EIGEN$vectors[,1]
          # grad(t) %*% EIGEN$vectors[,2] == 0

          # closest rige-line point is actually in cell
          TEST <- all(abs(dr)<c(dx,dy)/2)
          if(TEST)
          {
            m <- m + 1

            # middling ridge point
            r <- c(object$r$x[i],object$r$y[j]) + dr
            # store ridge point
            data[m,axes] <- r

            # gradient at middling ridge point
            # t <- which.max(abs(DIR))
            # t <- DR[t]/DIR[t]
            # grad <- (grad + t*hess) * DIR
            # # store as velocity
            # data[m,paste0("v",axes)] <- grad

            grad <- GRAD + c(HESS%*%dr) # gradient on ridge line
            DIR <- EIGEN$vectors[,1] # direction along ridge line
            data[m,paste0("v",axes)] <- sign(c(grad%*%DIR)) * DIR

            # # GRF adjusted Hessian
            # iCOV <- GROW %*% HESS %*% t(GROW)
            #
            # # CDF/N heuristic correction
            # SUB <- CDF[i+(-1):1,j+(-1):1]
            # STUFF <- DiffGrid(SUB,dx,dy,dx2,dy2,dr2,dxdydr2)
            # GRAD <- STUFF$GRAD
            # HESS <- STUFF$HESS
            # dr <- dr + DR
            # cdf <- CDF[i,j] + c(GRAD %*% dr) + c(dr %*% HESS %*% dr)/2
            # cdf <- clamp(cdf,0,1)
            #
            # Fstat <- EIGEN$values[1]/EIGEN$values[2]
            #
            # VAR <- VAR / (1-cdf) # DOES THIS MAKE SENSE?

            # ridge-only covariance
            # DIR <- STUFF$vectors[,2]
            # COV <- VAR * (DIR %o% DIR)
            # GRF mode uncertainty ellipse
            # COV <- SCALE %*% COV %*% t(SCALE)
            # VAR <- eigen(COV,only.values=TRUE)$values[1]

            # data[m,"VAR.xy"] <- VAR
          } # end if ridge-line passes through cell
        } # end if some negative curvature
      } # end if non-zero density
    } # end col loop
  } # end row loop

  # temporary fix
  data$t <- (1:nrow(data)) * sqrt(dr2)

  data <- new.telemetry(data,info=CTMM@info)
  return(data)
}
