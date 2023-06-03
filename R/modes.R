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
        COV <- PDsolve(-HESS)

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


