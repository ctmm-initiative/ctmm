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


#
ridges.UD <- function(object,level.UD=0.95,precision=1/8,...)
{
  tol <- .Machine$double.eps^precision
  VOTE.MIN <- 2
  LENGTH.MIN <- 2

  PDF <- object$PDF
  PDF <- log(PDF) # log transform for more accurate derivative estimation

  dx <- object$dr['x']
  dy <- object$dr['y']
  dx2 <- dx^2
  dy2 <- dy^2
  dr2 <- dx^2+dy^2
  dxdydr2 <- dx*dy/dr2

  DIM <- dim(PDF)
  VOTE <- array(0,DIM)

  X <- object$r$x
  Y <- object$r$y
  dR <- c(dx,dy)

  if(all(c('grad','hess') %in% names(object)))
  {
    GRAD <- object$grad
    HESS <- object$hess
  }
  else
  {
    # gradient information (under log)
    GRAD <- array(0,c(DIM,2)) # akima won't tolerate -Inf nor NA
    HESS <- array(0,c(DIM,2,2))

    ROWS <- 2:(nrow(PDF)-1)
    for(i in ROWS)
    {
      COLS <- 2:(ncol(PDF)-1)
      for(j in COLS)
      {
        SUB <- PDF[i+(-1):1,j+(-1):1]
        TEST <- all(SUB>-Inf)
        if(TEST) # non-zero density (under logarithm)
        {
          DIFF <- DiffGrid(SUB,dx,dy,dx2,dy2,dr2,dxdydr2)
          # these are better estimates than from splines
          GRAD[i,j,] <- DIFF$GRAD # 1/p D.p # from log
          HESS[i,j,,] <- DIFF$HESS # 1/p D.D.p - 1/p^2 D.p D.p # from log
          HESS[i,j,,] <- HESS[i,j,,] + (GRAD[i,j,] %o% GRAD[i,j,]) # logarithm correction # now 1/p D.D.p
          # last part is unnecessary for modes, but for ridges this looks important
        } # end if non-zero density
      } # end col loop
    } # end row loop
  }

  dim(HESS) <- c(DIM,4)
  HESS <- HESS[,,c(1,2,4)] # xx, xy, yy

  IND <- sort(object$CDF,decreasing=FALSE,index.return=TRUE)
  CDF <- IND$x
  i <- 0
  while(CDF[i+1]<=level.UD) { i <- i + 1 }
  rm(CDF)
  IND <- IND$ix[1:i]
  IND <- arrayInd(IND,DIM) # [sort,(row,col)]

  # ridge coordinates
  RIDGE <- matrix(0,nrow(IND),2)
  # is a mode
  MODE <- array(FALSE,nrow(IND))
  # CURV <- array(0,nrow(IND))

  pixels <- function(dr,DIR)
  {
    SIDE <- rbind(right=dr,left=dr,top=dr,bottom=dr)
    colnames(SIDE) <- c('x','y')

    # right solution
    dir <- DIR/sign(DIR[1])
    dir <- nant(dir,0)
    t <- ( +dx/2 - dr[1] ) / dir[1]
    SIDE['right',] <- dr + t*dir

    # left solution
    dir <- -dir
    t <- ( -dx/2 - dr[1] ) / dir[1]
    SIDE['left',] <- dr + t*dir

    # top solution
    dir <- DIR/sign(DIR[2])
    dir <- nant(dir,0)
    t <- ( +dy/2 - dr[2] ) / dir[2]
    SIDE['top',] <- dr + t*dir

    # bottom solution
    dir <- -dir
    t <- ( -dy/2 - dr[2] ) / dir[2]
    SIDE['bottom',] <- dr + t*dir

    SIDE <- t(t(SIDE)/dR*2)

    # first intersections
    PIX <- matrix(FALSE,3,3)

    if(all(abs(dr)<=dR/2))
    { PIX[2,2] <- TRUE }
    else
    { return(PIX) }

    if(abs(SIDE['left','y'])<1) { PIX[1,2] <- TRUE }
    else if(SIDE['left','y']==+1) { PIX[1,3] <- TRUE }
    else if(SIDE['left','y']==-1) { PIX[1,1] <- TRUE }

    if(abs(SIDE['right','y'])<1) { PIX[3,2] <- TRUE }
    else if(SIDE['right','y']==+1) { PIX[3,1] <- TRUE }
    else if(SIDE['right','y']==-1) { PIX[3,3] <- TRUE }

    if(abs(SIDE['bottom','x'])<1) { PIX[2,1] <- TRUE }
    if(abs(SIDE['top','x'])<1) { PIX[2,3] <- TRUE }

    # does cross the middle pixel
    if(!PIX[2,2] && sum(apply(PIX,1,any))>1 && sum(apply(PIX,2,any))) { PIX[2,2] <- TRUE }

    return(PIX)
  }

  pb <- utils::txtProgressBar(style=3)
  for(i in 1:nrow(IND))
  {
    if(any(IND[i,]==1) || any(IND[i,]==DIM)) { next }

    SUB.x <- IND[i,1]+(-1):1
    SUB.y <- IND[i,2]+(-1):1
    SUB <- PDF[SUB.x,SUB.y]

    M.TEST <- all(SUB[5] > SUB[-5]) # SUB[5] == P[i,j] # mode
    R.TEST <- min(SUB[-5]) < SUB[5] && SUB[5] <= max(SUB[-5]) # possible ridge

    if(!M.TEST && !R.TEST) { next }

    G <- GRAD[IND[i,1],IND[i,2],]
    H <- matrix(HESS[IND[i,1],IND[i,2],c(1,2,2,3)],2,2)
    EIGEN <- eigen(H)
    DIR <- EIGEN$vectors[,1]

    r <- c(X[IND[i,1]],Y[IND[i,2]])

    # first step without interpolation
    if(M.TEST)
    {
      MODE[i] <- TRUE
      dr <- c(PDsolve(-H) %*% G)
    }
    else
    {
      if(EIGEN$values[2]>=0) { next }

      # closest ridge point
      dr <- -c(G %*% EIGEN$vectors[,2])/EIGEN$values[2]
      dr <- dr * EIGEN$vectors[,2]
    }
    if(any(abs(dr)>1.5*dR)) { next } # no ridge in vicinity

    PIX <- pixels(dr,DIR)
    VOTE[SUB.x,SUB.y] <- VOTE[SUB.x,SUB.y] + PIX
    # VOTE[IND[i,,drop=FALSE]] <- 1

    # store results
    r <- r + dr
    RIDGE[i,] <- r

    G <- G + H %*% dr
    # CURV[i] <- -EIGEN$values[2]/max(-EIGEN$values[1],0)

    utils::setTxtProgressBar(pb,i/nrow(IND))
  } # end loop over IND

  colnames(RIDGE) <- c('x','y')
  VOTE <- VOTE[IND]
  rm(IND)
  VOTE <- (VOTE>=VOTE.MIN | MODE)
  MODE <- MODE[VOTE]
  RIDGE <- RIDGE[VOTE,,drop=FALSE]
  rm(VOTE)

  ## connect ridge points from top to bottom ##
  LINE <- list()
  for(i in 1:nrow(RIDGE))
  {
    if(MODE[i]>0) # start new ridge-line from mode
    { LINE[[length(LINE)+1]] <- RIDGE[i,,drop=FALSE] }
    else # point should connect to existing ridge-line
    {
      r <- RIDGE[i,c('x','y')]
      DIST <- sapply(LINE,function(L){sum((r-tail(L,1))^2)})
      j <- which.min(DIST) # ridge-line to connect to

      # greedy search up tail end
      n <- nrow(LINE[[j]])
      DIST <- DIST[j]
      for(k in n:1)
      {
        if(k==n) { next }
        D.NEXT <- sum((r-LINE[[j]][k,])^2)
        if(DIST<=D.NEXT)
        {
          k <- k+1 # last point was best point
          break
        }
        DIST <- D.NEXT
      }

      DIST <- abs(r-LINE[[j]][k,])
      if(any(DIST>2*dR)) # far from other points
      { LINE[[length(LINE)+1]] <- RIDGE[i,,drop=FALSE] }
      if(k==n) # connect to tail end
      { LINE[[j]] <- rbind(LINE[[j]],RIDGE[i,]) }
      else # branch off and make new ridge
      { LINE[[length(LINE)+1]] <- rbind(LINE[[j]][k,],RIDGE[i,]) }
    } # end connect point to existing ridge
  } # end connect ridge points
  rm(RIDGE)

  ## connect end points into loops
  n <- length(LINE)
  LOOP <- list()

  if(n>1) # connect loops
  {
    DIST2 <- array(Inf,c(n,n,2))
    for(i in 1%:%(n-1))
    {
      for(j in (i+1)%:%n)
      { DIST2[i,j,] <- DIST2[j,i,] <- abs( tail(LINE[[i]],1) - tail(LINE[[j]],1) ) }
    }
    DIST <- DIST2[,,1]^2+DIST2[,,2]^2

    while(dim(DIST)[1]>1)
    {
      IND <- which.min(DIST)
      IND <- arrayInd(IND,dim(DIST))

      if(all(DIST2[IND[1],IND[2],]<=2*dR))
      {
        m <- nrow(LINE[[IND[1]]])
        LOOP[[1+length(LOOP)]] <- rbind( LINE[[IND[2]]] , LINE[[IND[1]]][m:1,] )
        LINE <- LINE[-c(IND)]
        DIST2 <- DIST2[-c(IND),-c(IND),,drop=FALSE]
        DIST <- DIST[-c(IND),-c(IND),drop=FALSE]
      }
      else
      { break }
    }
  } # end multiple ridges

  LINE <- c(LINE,LOOP)
  rm(LOOP)

  # LENGTH <- sapply(LINE,nrow)
  # LINE <- LINE[LENGTH>=LENGTH.MIN]

  LINE <- lapply(LINE,function(L){new.telemetry(data.frame(cbind(L,t=1:nrow(L))))})

  close(pb)
  return(LINE)
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
