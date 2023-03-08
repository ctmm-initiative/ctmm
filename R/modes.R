# statistical mode
modes.numeric <- function(x,method="exact",na.rm=TRUE)
{
  if(na.rm) { x <- x[!is.na(x)] }

  x <- sort(x)
  unique.x <- unique(x)

  if(length(unique.x)<length(x))
  {
    freq <- tabulate(match(x, unique.x))
    MAX <- max(freq)
    IND <- freq==MAX
    MODES <- unique.x[IND]

    # well defined sample mode
    if(length(MODES)==1) { return(MODES) }

    # else minimize log-distance
    IND <- x %in% MODES
  }
  else
  { IND <- rep(TRUE,length(x)) }

  # log distance
  D <- log(abs(outer(x,x,'-')))
  # no self distance
  diag(D) <- 0
  # no ties
  if(sum(IND)<length(x))
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
  GRAD[1] <- mean( c( diff(SUB[1:2,2]) , diff(SUB[2:3,2]) ) )/dx
  GRAD[2] <- mean( c( diff(SUB[2,1:2]) , diff(SUB[2,2:3]) ) )/dy

  # calculate Hessian
  HESS <- diag(2)
  HESS[1,1] <- diff(diff(SUB[,2]))/dx2 # (+x,+x)
  HESS[2,2] <- diff(diff(SUB[2,]))/dy2 # (+y,+y)
  HESS[1,2] <- diff(diff(c(SUB[1,3],SUB[2,2],SUB[3,1])))/dr2 # (+x+y,+x+y)
  HESS[2,1] <- diff(diff(c(SUB[1,1],SUB[2,2],SUB[3,3])))/dr2 # (+x-y,+x-y)
  HESS[1,2] <- HESS[2,1] <- dxdydr2 * (HESS[2,1]-HESS[1,2]) # (+x,+y)

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
  PDF <- log(PDF) # log transform for numerics
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

        STUFF <- DiffGrid(SUB,dx,dy,dx2,dy2,dr2,dxdydr2)
        GRAD <- STUFF$GRAD
        HESS <- STUFF$HESS
        HESS <- HESS + (GRAD %o% GRAD) # logarithm correction
        COV <- -PDsolve(HESS)

        # locate mode
        dr <- c(COV %*% GRAD)

        # is mode well estimated as within pixel?
        if(dr[1] > -dx/2 && dr[1] < dx/2 && dr[2] > -dy/2 && dr[2] < dy/2)
        {
          r <- c(object$r$x[i],object$r$y[j]) + dr
          # store mode
          data[m,axes] <- r

          # GRF mode uncertainty ellipse
          COV <- SCALE %*% COV %*% t(SCALE)

          # store uncertainty ellipse
          for(I in 1:length(axes))
          {
            for(J in I:length(axes))
            {
              NAME <- paste("COV.",axes[I],".",axes[J],sep="") # consistent with imported ARGOS error ellipse notation
              data[m,NAME] <- COV[I,J]
            }
          } # end uncertainty ellipse
        } # end mode recording
      } # end if local minimum
    } # end col loop
  } # end row loop

  data <- new.telemetry(data,info=CTMM@info)
  FIX <- data
  uere(FIX) <- 1
  data@UERE <- FIX@UERE
  return(data)
}


# rewrite to work from CDF, via least steep descent

#
ridges.UD <- function(object,...)
{
  # if(all(c('grad','hess') %nin% names(object))) { stop("Please run akde() with grad=TRUE") }

  PDF <- object$PDF
  PDF <- log(PDF) # log transform
  # COST <- array(Inf,dim(PDF))
  POINT <- array(0,dim(PDF))
  dx <- object$dr['x']
  dy <- object$dr['y']
  dx2 <- dx^2
  dy2 <- dy^2
  dr2 <- dx^2+dy^2
  dxdydr2 <- dx*dy/dr2

  ROWS <- 2:(nrow(PDF)-1)
  # ROWS <- 1:nrow(PDF)
  for(i in ROWS)
  {
    COLS <- 2:(ncol(PDF)-1)
    # COLS <- 1:ncol(PDF)
    for(j in COLS)
    {
      SUB <- PDF[i+(-1):1,j+(-1):1]
      TEST <- all(SUB>-Inf)
      # TEST <- PDF[i,j]>0
      if(TEST) # non-zero density (under logarithm)
      {
        # GRAD <- object$grad[i,j,]
        # HESS <- object$hess[i,j,,]

        DIFF <- DiffGrid(SUB,dx,dy,dx2,dy2,dr2,dxdydr2)
        GRAD <- DIFF$GRAD # 1/p D.p # from log
        HESS <- DIFF$HESS # 1/p D.D.p - 1/p^2 D.p D.p # from log
        HESS <- HESS + (GRAD %o% GRAD) # logarithm correction # now 1/p D.D.p

        EIGEN <- eigen(HESS)
        if(EIGEN$values[2]<0) # most negative eigenvalue is negative # a ridge exists somewhere
        {
          # COST[i,j] <- (GRAD %*% EIGEN$vectors[,2])^2 / EIGEN$values[2]
          # if(EIGEN$values[2]>0) { COST[i,j] <- COST[i,j] + EIGEN$values[1]/EIGEN$values[2] }

          # ridge constraint
          # (GRAD + HESS %*% DL ) %*% EIGEN$vectors[,1] == 0

          # 1. Lagrangian (minimal distance)
          # L == 1/2 * DL %*% DL - lambda * (GRAD + HESS %*% DL ) %*% EIGEN$vectors[,2]
          # 0 == DL - lambda * HESS %*% EIGEN$vectors[,2]
          # 0 == DL - lambda * EIGEN$values[2] * EIGEN$vectors[,2]
          # DL == lambda * EIGEN$values[2] * EIGEN$vectors[,2]

          # 2. climbing directly up to the ridge
          # DL == dl * EIGEN$vectors[,2]
          # GRAD %*% EIGEN$vectors[,2] + dl*(EIGEN$vectors[,2] %*% HESS %*% EIGEN$vectors[,2]) == 0
          # GRAD %*% EIGEN$vectors[,2] + EIGEN$values[2]*dl == 0

          dr <- -c(GRAD %*% EIGEN$vectors[,2])/EIGEN$values[2]
          dr <- dr * EIGEN$vectors[,2]

          # pixel displacement
          dij <- dr/c(dx,dy)
          ip <- i+dij[1]
          jp <- j+dij[2]

          # is point in an adjacent pixel or this pixel
          if(max(abs(dij))<=2)
          {
            # point is in this pixel
            if(max(abs(dij))<=1)
            {
              POINT <- anti.alias(ip,jp,POINT,1)

              # what pixels are next?
              u <- EIGEN$vectors[,1] # bi-direction of ridge
              u <- u/c(dx,dy) # now in units of pixels
              u <- u/sqrt(sum(u^2)) # one pixel length

              POINT <- anti.alias(ip+u[1],jp+u[2],POINT,1/2)
              POINT <- anti.alias(ip-u[1],jp-u[2],POINT,1/2)
            } # point is in this pixel
            else # point is in adjacent pixel
            { POINT <- anti.alias(ip,jp,POINT,1/8) }
          } # point is in nearby pixel
        } # end ridge exists
      } # end if non-zero density
    } # end col loop
  } # end row loop

  # CTMM <- object@CTMM
  # DOF <- DOF.mean(CTMM)
  # COST <- object$PDF * COST # cancel 1/p term
  # COST <- DOF* COST # now like grad^2/std.err^2 ~ F

  POINT <- POINT/(1+2/2+6/8)
  POINT[POINT>1] <- 1
  return(POINT)
}

anti.alias <- function(i,j,M,dM=1)
{
  I <- unique( c(floor(i),ceiling(i)) )
  J <- unique( c(floor(j),ceiling(j)) )

  I <- I[0<I & I<nrow(M)]
  J <- J[0<J & J<ncol(M)]

  dM <- array(dM,c(length(I),length(J)))

  if(!length(dM)) { return(M) }

  if(length(I)>1)
  {
    w <- abs(i-I[1])
    dM[1,] <- w*dM[1,]
    dM[2,] <- (1-w)*dM[2,]
  }

  if(length(J)>1)
  {
    w <- abs(j-J[1])
    dM[,1] <- w*dM[,1]
    dM[,2] <- (1-w)*dM[,2]
  }

  M[I,J] <- M[I,J] + dM
  return(M)
}
