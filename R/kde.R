# akde object generator
# list of kde objects with info slots
new.UD <- methods::setClass("UD", representation("list",info="list"))


# Slow lag counter
tab.lag.DOF <- function(data,fast=NULL,dt=NULL,w=NULL)
{
  t <- data$t
  # intelligently select algorithm
  n <- length(t)

  # uniform weights
  # if(is.null(w)) { w <- rep(1,n)/n }
  # finish this....
  
  if(is.null(fast))
  {
    if(n<100) { fast <- FALSE }
    else { fast <- TRUE }
  }
  
  # calculate lag,n(lag) vectors
  if(fast)
  {
    lag.DOF <- variogram.fast(data,dt=dt,CI="IID",axes="t") # count all lag pairs
    lag.DOF$SVF <- NULL
    
    # lag==0 should not be doubled
    lag.DOF$DOF[1] <- lag.DOF$DOF[1]/2 
  }
  else
  {
    lag <- outer(t,t,FUN="-")
    lag <- abs(lag) # may have to change far in the future
    lag <- c(lag) # collapse to 1D array
    
    # can we think of a faster sorter?
    r <- lag # to be the remainder
    lag <- NULL
    DOF <- NULL
    while(length(r)>0)
    {
      n <- length(r)
      lag <- c(lag,r[1])
      r <- r[r!=r[1]]
      DOF <- c(DOF,n-length(r))
    }
    
    lag.DOF <- list(lag=lag,DOF=DOF)
  }
  
  return(lag.DOF)
}


##################################
# Bandwidth optimizer
#lag.DOF is an unsupported option for end users
akde.bandwidth <- function(data,CTMM,fast=NULL,dt=NULL,verbose=FALSE)
{
  
  lag.DOF <- tab.lag.DOF(data,fast=fast,dt=dt)
  # Extract lag data
  DOF <- lag.DOF$DOF
  lag <- lag.DOF$lag
  
  sigma <- methods::getDataPart(CTMM$sigma)
  # standardized SVF
  CTMM$sigma <- covm(diag(1,2))
  svf <- svf.func(CTMM,moment=FALSE)$svf
  g <- Vectorize(svf)
  
  n <- length(data$t)
  
  # Mean Integrated Square Error modulo a constant
  MISE <- function(h)
  {
    if(h<=0) {Inf}
    else { (1/n^2)*sum(DOF/(2*g(lag)+2*h^2)) - 2/(2+h^2) + 1/2 }
  }
  
  h <- 1/n^(1/6) # User Silverman's rule of thumb to place lower bound
  h <- stats::optimize(f=MISE,interval=c(h/2,2))$minimum
  
  H <- h^2
  
  DOF.H <- ( 1/(2*H)^2 - 1/(2+2*H)^2 ) / ( 1/(2+H)^2 - 1/(2+2*H)^2 )
  
  H <- H*sigma
  
  rownames(H) <- c("x","y")
  colnames(H) <- c("x","y")

  CTMM$sigma <- covm(sigma)
  bias <- akde.bias(CTMM,H=H,lag=lag,DOF=DOF)
  
  if(verbose)
  { return(list(H=H,DOF.H=DOF.H,bias=bias,DOF.area=DOF.area(CTMM))) }
  else
  { return(H) }
}

# bias of Gaussian Reference function AKDE
# generalize to non-stationary mean
akde.bias <- function(CTMM,H,lag,DOF)
{
  sigma <- methods::getDataPart(CTMM$sigma)
  
  # weighted correlation
  ACF <- Vectorize( svf.func(CTMM,moment=FALSE)$ACF )
  COV <- sum(ACF(lag)*DOF)/DOF[1]^2
  
  # area inflation factor      
  bias <- sqrt( det( (1-COV)*sigma + H )/det(sigma) )
  
  return(bias)
}

###################
# homerange wrapper function
homerange <- function(data,CTMM,method="AKDE",...)
{
  akde(data,CTMM,...)
}


#######################################
# wrap the kde function for our telemetry data format and CIs.
akde <- function(data,CTMM,debias=TRUE,smooth=TRUE,error=0.001,res=10,grid=NULL,...)
{
  # smooth out errors
  if(CTMM$error && smooth) { data <- predict(CTMM,data=data,t=data$t) }
  
  # calculate optimal bandwidth and some other information
  KDE <- akde.bandwidth(data=data,CTMM=CTMM,verbose=TRUE,...)
  if(debias) { debias <- KDE$bias }

  # absolute resolution
  dr <- sqrt(c(KDE$H[1,1],KDE$H[2,2]))/res

  KDE <- c(KDE,kde(data,KDE$H,bias=debias,alpha=error,dr=dr,grid=grid))

  KDE <- new.UD(KDE,info=attr(data,"info"))
  
  return(KDE)
}


# if a single H matrix is given, make it into an array of H matrices
prepare.H <- function(data,H)
{
  n <- length(data$x)

  # one matrix given
  if(length(dim(H))==2)
  {
    H <- array(H,c(2,2,n))
    H <- aperm(H,c(3,1,2))
  }
  
  return(H)
}


############################################
# construct a grid for the density function
kde.grid <- function(data,H,alpha=0.001,res=1,dr=NULL)
{
  x <- data$x
  y <- data$y
  
  H <- prepare.H(data,H)
  
  # how far to extend range from data as to ensure alpha significance in total probability
  z <- sqrt(-2*log(alpha))
  DX <- z * apply(H,1,function(h){sqrt(h[1,1])})
  DY <- z * apply(H,1,function(h){sqrt(h[2,2])})
  
  # now to find the necessary extent of our grid
  rx <- c( min(x - DX) , max(x + DX) )
  ry <- c( min(y - DY) , max(y + DY) )
  
  # grid center
  mu.x <- mean(rx)
  mu.y <- mean(ry)
  
  if(is.null(dr))
  {
    dr <- NULL
    
    # grid resolution 
    dr[1] <- diff(rx)/res
    dr[2] <- diff(ry)/res
    
    # grid locations
    res <- res/2+1 # half resolution
    X <- mu.x + (-res):(res)*dr[1]
    Y <- mu.y + (-res):(res)*dr[2]
  }
  else
  {
    res <- diff(rx)/dr[1]
    res <- ceiling(res/2+1)
    X <- mu.x + (-res):(res)*dr[1]
    
    res <- diff(ry)/dr[2]
    res <- ceiling(res/2+1)
    Y <- mu.y + (-res):(res)*dr[2]
  }
  
  grid <- list(x=X,y=Y,dr=dr,DX=DX,DY=DY,rx=rx,ry=ry)
  return(grid)
}


##################################
# construct my own kde objects
# was using ks-package but it has some bugs
# alpha is the error goal in my total probability
kde <- function(data,H,bias=FALSE,W=rep(1,length(data$x)),alpha=0.001,res=NULL,dr=NULL,grid=NULL)
{
  n <- length(data$x)
  x <- data$x
  y <- data$y

  # normalize weights
  W <- W/sum(W)
  
  # format bandwidth matrix
  H <- prepare.H(data,H)
  
  if(is.null(grid)) { grid <- kde.grid(data,H=H,alpha=alpha,res=res,dr=dr) }
  X <- grid$x
  Y <- grid$y
  # generalize this for future grid option use
  DX <- grid$DX
  DY <- grid$DY
  dx <- grid$dr[1]
  dy <- grid$dr[2]
  
  cdf <- array(0,c(length(X),length(Y)))
  for(i in 1:n)
  {
    # sub-grid row indices
    r1 <- floor((x[i]-DX[i]-X[1])/dx) + 1
    r2 <- ceiling((x[i]+DX[i]-X[1])/dx) + 1
    
    # sub-grid column indices
    c1 <- floor((y[i]-DY[i]-Y[1])/dy) + 1
    c2 <- ceiling((y[i]+DY[i]-Y[1])/dy) + 1
    
    cdf[r1:r2,c1:c2] <- cdf[r1:r2,c1:c2] + W[i]*pnorm2(X[r1:r2]-x[i],Y[c1:c2]-y[i],H[i,,],dx,dy,alpha)
  }

  dA <- dx*dy
  if(!bias) { pdf <- cdf/dA }

  # cdf: cell probability -> probability included in contour
  cdf <- pdf2cdf(cdf,finish=FALSE)
  DIM <- cdf$DIM
  IND <- cdf$IND
  cdf <- cdf$cdf
  
  # areas are biased to be estimated as debias*area
  if(bias)
  {
    # counting area by dA
    AREA <- 1:length(cdf)

    # evaluate the debiased cdf on the original area grid
    cdf <- stats::approx(x=AREA/bias,y=cdf,xout=AREA,yleft=0,yright=1)$y
    
    # recalculate pdf
    pdf <- diff(c(0,cdf))/dA
    pdf[IND] <- pdf
    pdf <- array(pdf,DIM)
  }

  cdf[IND] <- cdf # back in spatial order
  cdf <- array(cdf,DIM) # back in table form
  
  result <- list(PDF=pdf,CDF=cdf,x=X,y=Y,dA=dA)
  
  return(result)
}

########################
# cdf: cell probability -> probability included in contour
pdf2cdf <- function(cdf,finish=TRUE)
{
  #cdf <- pdf * dA
  DIM <- dim(cdf)
  cdf <- c(cdf) # flatten table
  cdf <- sort(cdf,decreasing=TRUE,method="quick",index.return=TRUE)
  IND <- cdf[[2]] # sorted indices
  cdf <- cdf[[1]]
  cdf <- cumsum(cdf)

  if(finish)
  {
    cdf[IND] <- cdf # back in spatial order
    cdf <- array(cdf,DIM) # back in table form
    return(cdf)
  }
  else
  { return(list(cdf=cdf,IND=IND,DIM=DIM)) }
}

#######################
# robust bi-variate CDF (mean zero assumed)
#######################
pnorm2 <- function(X,Y,sigma,dx=stats::mean(diff(X)),dy=stats::mean(diff(Y)),alpha=0.001)
{
  cdf <- array(0,c(length(X),length(Y)))
  
  # eigensystem of kernel covariance
  v <- eigen(sigma)
  s <- v$values
  
  # effective degree of degeneracy at best resolution
  ZERO <- sum(s <= 0 | (min(dx,dy)/2)^2/s > -2*log(alpha))

  # correlation
  S <- sqrt(sigma[1,1]*sigma[2,2])
  if(S>0)
  {
    rho <- sigma[1,2]/S
    # prevent some tiny numerical errors just in case
    rho <- clamp(rho,min=-1,max=1)
  }
  else { rho <- 0 }
  
  # main switch
  if(ZERO==0 && abs(rho)<1) # no degeneracy
  {
    # relative grid size (worst case)
    z <- sqrt((dx^2+dy^2)/s[2])
    
    if(z^3/12<=alpha) # midpoint integration
    {
      cdf <- (dx*dy) * Gauss(X,Y,sigma)
    }
    else if(z^5/2880<=alpha) # Simpson integration
    {
      W <- c(1,4,1)
      cdf <- NewtonCotes(X,Y,sigma,W,dx,dy)
    }
    else if(z^7/1935360<=alpha) # Boole integration
    {
      W <- c(7,32,12,32,7)
      cdf <- NewtonCotes(X,Y,sigma,W,dx,dy)
    }
    else # exact calculation
    {
      # offset to corners
      x <- c(X-dx/2,last(X)+dx/2)
      y <- c(Y-dy/2,last(Y)+dy/2)
      
      # standardized all locations
      x <- (x)/sqrt(sigma[1,1])
      y <- (y)/sqrt(sigma[2,2])
      
      # dimensions
      n.x <- length(x)
      n.y <- length(y)
      
      # corner grid of cell probabilities
      CDF <- outer(x,y,function(x,y){pbivnorm::pbivnorm(x,y,rho=rho)})
      
      # integrate over cell and add
      cdf <- CDF[-1,-1] - CDF[-n.x,-1] - CDF[-1,-n.y] + CDF[-n.x,-n.y]
      
      # pbivnorm is very fragile
      # if(any(is.nan(cdf))) { stop(" dx=",dx," dy=",dy," sigma=",sigma," x=",x," y=",y," rho=",rho," CDF=",cdf)}
    }
  }
  else if(ZERO==1 || abs(rho)==1) # line degeneracy
  {
    # max variance
    s <- s[1]
    # unit vector of max variance
    v <- v$vectors[,1]
    
    # crossings along X grid
    x.cell <- c()
    y.cross <- c()
    m.y <- v[2]/v[1]
    if(abs(m.y)<Inf)
    {
      x.cell <- c(X-dx/2,last(X)+dx/2)
      y.cross <- (x.cell)*m.y
    }
    
    # crossings along Y grid
    y.cell <- c()
    x.cross <- c()
    m.x <- v[1]/v[2]
    if(abs(m.x)<Inf)
    {
      y.cell <- c(Y-dy/2,last(Y)+dy/2)
      x.cross <- (y.cell)*m.x
    }
    
    # all crossings
    x.cross <- c(x.cell,x.cross)
    y.cross <- c(y.cross,y.cell)
    
    # standardized location along line
    z.cross <- ((x.cross)*v[1] + (y.cross)*v[2])/sqrt(s)
    z.cross <- sort(z.cross,method="quick")
    z.cross <- unique(z.cross)
    
    for(i in 1:(length(z.cross)-1))
    {
      # what cell is this line segment in?
      z.mid <- mean(z.cross[i:(i+1)])
      
      x <- sqrt(s)*v[1]*z.mid
      y <- sqrt(s)*v[2]*z.mid
      
      r <- abs(x-X) <= dx/2
      c <- abs(y-Y) <= dy/2
      
      cdf[r,c] <- (stats::pnorm(z.cross[i+1])-stats::pnorm(z.cross[i]))/(sum(r)*sum(c))
    }
  }
  else if(ZERO==2) # point degeneracy
  {
    # the closest point(s)
    r <- abs(X) <= dx/2
    c <- abs(Y) <= dy/2
    
    # increment the closest point(s)
    cdf[r,c] <- cdf[r,c] + 1/(sum(r)*sum(c))
  }
  else stop("something is wrong in matrix: sigma == ",sigma)
  
  return(cdf)
}

#################
# Newton-Cotes integrators
NewtonCotes <- function(X,Y,sigma,W,dx=mean(diff(X)),dy=mean(diff(Y)))
{
  W <- W/sum(W)
  n <- length(W)
  m <- n-1
  
  # refined grid to have Simpson's rule in between X,Y points
  # changed from to= to length.out= to avoid roundoff error
  x <- seq(from=X[1]-dx/2,by=dx/m,length.out=length(X)*m+1)
  y <- seq(from=Y[1]-dy/2,by=dy/m,length.out=length(Y)*m+1)
  
  # weight arrays
  w.x <- dx * array(W[-n],length(x))
  w.y <- dy * array(W[-n],length(y))
  
  # weight table
  W <- (w.x %o% w.y)

  cdf <- W * Gauss(x,y,sigma)

  # coarsen grid
  # index order is (x,y)
  cdf <- vapply(1:length(X)-1,function(i){colSums(cdf[1:n+m*i,])},rep(0,length(y)))
  # index order is (y,X)
  cdf <- vapply(1:length(Y)-1,function(i){colSums(cdf[1:n+m*i,])},rep(0,length(X)))
  # index order is (X,Y)
  
  return(cdf)
}

#####################
# gaussian pdf
Gauss <- function(X,Y,sigma=NULL,sigma.inv=solve(sigma),sigma.GM=sqrt(det(sigma)))
{
  cdf <- outer(X^2*sigma.inv[1,1],Y^2*sigma.inv[2,2],"+")/2
  cdf <- cdf + (X %o% Y)*sigma.inv[1,2]
  cdf <- exp(-cdf)/(2*pi*sigma.GM)
  return(cdf)
}


#####################
# AKDE CIs
CI.UD <- function(object,level.UD=0.95,level=0.95,P=FALSE)
{
  # point estimate
  if(is.na(level.UD))
  {
    # calculate mean area
    area <- sort(object$PDF,decreasing=TRUE,method="quick")
    area <- sum(area * (1:length(area))) * object$dA
  }
  else
  { area <- sum(object$CDF <= level.UD) * object$dA }
  
  # chi square approximation of uncertainty
  area <- chisq.ci(area,DOF=2*object$DOF.area,alpha=1-level)
  
  if(!P) { return(area) }
  
  # probabilities associated with these areas
  P <- round(area / object$dA)
  P <- sort(object$CDF,method="quick")[P]

  # recorrect point estimate level
  if(!is.na(level.UD)) { P[2] <- level.UD }
  
  return(P)
}

#######################
# summarize details of akde object
summary.UD <- function(object,level.UD=0.95,level=0.95,...)
{
  area <- CI.UD(object,level.UD,level)
  
  # pretty units
  unit.info <- unit(area[2],"area")
  name <- unit.info$name
  scale <- unit.info$scale
  
  area <- array(area/scale,c(1,3))
  rownames(area) <- paste("area (",name,")",sep="")

  colnames(area) <- c("low","ML","high")

  SUM <- list()
  
  SUM$DOF <- c(object$DOF.area,object$DOF.H)
  names(SUM$DOF) <- c("area","bandwidth")
  
  SUM$CI <- area
  
  return(SUM)
}
#methods::setMethod("summary",signature(object="UD"), function(object,...) summary.UD(object,...))


################################
# create a raster of the ML akde
raster.UD <- function(x,DF="CDF",...)
{
  UD <- x
  
  dx <- UD$x[2]-UD$x[1]
  dy <- UD$y[2]-UD$y[1]
  
  xmn <- UD$x[1]-dx/2
  xmx <- last(UD$x)+dx/2
  
  ymn <- UD$y[1]-dy/2
  ymx <- last(UD$y)+dy/2
  
  Raster <- raster::raster(t(UD[[DF]][,dim(UD[[DF]])[2]:1]),xmn=xmn,xmx=xmx,ymn=ymn,ymx=ymx,crs=attr(UD,"info")$projection)
  
  return(Raster)
}
methods::setMethod("raster",signature(x="UD"), function(x,DF="CDF",...) raster.UD(x,DF=DF,...))


################
# Is contour A inside contour B
inside <- function(A,B)
{
  result <- mode(sp::point.in.polygon(A$x,A$y,B$x,B$y))
  if(1<=result && result<=2) { return(1) } else { return(0) }
}


##############
SpatialPolygonsDataFrame.UD <- function(object,level.UD=0.95,level=0.95,...)
{
  UD <- object
  
  P <- CI.UD(UD,level.UD,level,P=TRUE)
  NAMES <- c("low","ML","high")
  
  ID <- paste(UD@info$identity," ",round(100*level.UD),"% ",NAMES,sep="")

  polygons <- list()
  for(i in 1:length(P))
  {
    CL <- grDevices::contourLines(x=UD$x,y=UD$y,z=UD$CDF,levels=P[i])
    
    # create contour heirarchy matrix (half of it)
    H <- array(0,c(1,1)*length(CL))    
    for(row in 1:length(CL))
    {
      for(col in row:length(CL))
      {
        H[row,col] <- inside(CL[[row]],CL[[col]]) 
      }
    }
    
    # number of contours that this contour is inside
    I <- rowSums(H)
    
    # if I is odd, then you are a hole inside a positive area
    hole <- is.odd(I)
    
    # polygon
    polygons[[i]] <- list()
    for(j in 1:length(CL))
    {
      polygons[[i]][[j]] <- sp::Polygon( cbind( CL[[j]]$x , CL[[j]]$y ) , hole=hole[j] )
    }

    # polygonS
    polygons[[i]] <- sp::Polygons(polygons[[i]],ID=ID[i])
  }
  names(polygons) <- ID

    # spatial polygons
  polygons <- sp::SpatialPolygons(polygons, proj4string=sp::CRS(attr(UD,"info")$projection))

  # spatial polygons data frame  
  data <- data.frame(name=rev(ID))
  rownames(data) <- rev(ID)
  polygons <- sp::SpatialPolygonsDataFrame(polygons,data)
  
  return(polygons)
}
#methods::setMethod("SpatialPolygonsDataFrame",signature(Sr="UD"), function(Sr,level.UD=0.95,level=0.95) SpatialPolygonsDataFrame.UD(Sr,level.UD=level.UD,level=level))


################
writeShapefile.UD <- function(object, folder, file=UD@info$identity, level.UD=0.95 ,level=0.95,  ...)
{
  UD <- object
  
  SP <- SpatialPolygonsDataFrame.UD(UD,level.UD=level.UD,level=level)
  
  rgdal::writeOGR(SP, dsn=folder, layer=file, driver="ESRI Shapefile",...)
}

