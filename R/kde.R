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
akde.bandwidth <- function(data,CTMM,VMM=NULL,fast=NULL,dt=NULL,verbose=FALSE)
{
  
  lag.DOF <- tab.lag.DOF(data,fast=fast,dt=dt)
  # Extract lag data
  DOF <- lag.DOF$DOF
  lag <- lag.DOF$lag
  
  # is CTMM a list of H,V models or an H model?
  # generic function to delineate an H&V list
  
  sigma <- methods::getDataPart(CTMM$sigma)
  # standardized SVF
  CTMM$sigma <- covm(diag(1,2))
  svf <- svf.func(CTMM,moment=FALSE)$svf
  g <- Vectorize(svf)
  
  n <- length(data$t)
  
  if(is.null(VMM))
  {
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
    
    rownames(H) <- CTMM$axes
    colnames(H) <- CTMM$axes
    
    if(verbose)
    {
      CTMM$sigma <- covm(sigma,axes=CTMM$axes,isotropic=CTMM$isotropic)
      bias <- akde.bias(CTMM,H=H,lag=lag,DOF=DOF)
      
      return(list(H=H,DOF.H=DOF.H,bias=bias,DOF.area=DOF.area(CTMM)))
    }
  }  
  else
  {
    sigmaz <- methods::getDataPart(VMM$sigma)
    # standardized SVF
    VMM$sigma <- covm(1,axes=VMM$axes,isotropic=VMM$isotropic)
    svfz <- svf.func(VMM,moment=FALSE)$svf
    gz <- Vectorize(svfz)
    
    # Mean Integrated Square Error modulo a constant
    MISE <- function(h)
    {
      if(any(h==0)) {Inf}
      else { (1/n^2)*sum(DOF/(2*g(lag)+2*h[1]^2)/sqrt(2*gz(lag)+2*h[2]^2)) - 2/(2+h[1]^2)/sqrt(2+h[2]^2) + 1/2^(3/2) }
    }
    
    h <- 4/5/n^(1/7) # User Silverman's rule of thumb to place lower bound
    h <- stats::optim(par=c(h,h),fn=MISE,control=list(maxit=.Machine$integer.max))$par
    
    H <- h^2
    
    # horizontal DOF
    DOF.H <- ( 1/(2*H[1])^2/sqrt(2*H[2]) - 1/(2+2*H[1])^2/sqrt(2+2*H[2]) ) / ( 1/(2+H[1])^2/sqrt(2+H[2]) - 1/(2+2*H[1])^2/sqrt(2+2*H[2]) )
    # vertical DOF
    DOF.H[2] <- ( 1/(2*H[1])/(2*H[2])^(3/2) - 1/(2+2*H[1])/(2+2*H[2])^(3/2) ) / ( 1/(2+H[1])/sqrt(2+H[2])^(3/2) - 1/(2+2*H[1])/sqrt(2+2*H[2])^(3/2) )
    # Above are bad
    DOF.H <- c(NA,NA)
    
    VH <- H[2]*sigmaz
    HH <- H[1]*sigma
    
    H <- matrix(0,3,3)
    H[1:2,1:2] <- HH
    H[3,3] <- VH
    
    axes <- c(CTMM$axes,VMM$axes)
    
    rownames(H) <- axes
    colnames(H) <- axes
    
    if(verbose)
    {
      CTMM$sigma <- covm(sigma,axes=CTMM$axes,isotropic=CTMM$isotropic)
      VMM$sigma <- covm(sigmaz,axes=VMM$axes,isotropic=TRUE)
      
      bias <- akde.bias(CTMM,H=HH,lag=lag,DOF=DOF)
      bias[3] <- akde.bias(VMM,H=VH,lag=lag,DOF=DOF)
      names(bias) <- axes
      
      DOF.area <- c(DOF.area(CTMM),DOF.area(VMM))
      
      axes <- c("Horizontal","Vertical")
      names(DOF.area) <- axes
      names(DOF.H) <- axes
      
      return(list(H=H,DOF.H=DOF.H,bias=bias,DOF.area=DOF.area))
    }
    
  }
    
  return(H) 
}

# bias of Gaussian Reference function AKDE
# generalize to non-stationary mean
akde.bias <- function(CTMM,H,lag,DOF)
{
  sigma <- methods::getDataPart(CTMM$sigma)
  
  # weighted correlation
  ACF <- Vectorize( svf.func(CTMM,moment=FALSE)$ACF )
  COV <- sum(ACF(lag)*DOF)/DOF[1]^2
  
  # variance inflation factor
  bias <- ( det(cbind((1-COV)*sigma + H))/det(cbind(sigma)) )^(1/length(CTMM$axes))
  # remove cbind if/when I can get det.numeric working
  
  # name dimensions of bias
  axes <- CTMM$axes
  bias <- rep(bias,length(axes))
  names(bias) <- axes
    
  return(bias)
}


###################
# homerange wrapper function
homerange <- function(data,CTMM,method="AKDE",...)
{ akde(data,CTMM,...) }


#######################################
# wrap the kde function for our telemetry data format and CIs.
akde <- function(data,CTMM,VMM=NULL,debias=TRUE,smooth=TRUE,error=0.001,res=10,grid=NULL,...)
{
  axes <- CTMM$axes
  
  # smooth out errors
  z <- NULL
  if(!is.null(VMM))
  {
    axis <- VMM$axes
    if(VMM$error && smooth) { z <- predict(VMM,data=data,t=data$t)[,axis] }
    axes <- c(axes,axis)
  }
  if(CTMM$error && smooth) { data <- predict(CTMM,data=data,t=data$t) }
  # copy back
  if(!is.null(VMM)) { data[,axis] <- z }
  
  # calculate optimal bandwidth and some other information
  KDE <- akde.bandwidth(data=data,CTMM=CTMM,VMM=VMM,verbose=TRUE,...)
  if(debias) { debias <- KDE$bias }

  # absolute resolution
  dr <- sqrt(diag(KDE$H))/res

  KDE <- c(KDE,kde(data,KDE$H,axes=axes,bias=debias,alpha=error,dr=dr,grid=grid))

  KDE <- new.UD(KDE,info=attr(data,"info"))
  
  return(KDE)
}


# if a single H matrix is given, make it into an array of H matrices
prepare.H <- function(H,n)
{
  # one matrix given
  if(length(dim(H))==2)
  {
    d <- nrow(H)
    H <- array(H,c(d,d,n))
    H <- aperm(H,c(3,1,2))
  }
  
  return(H)
}


############################################
# construct a grid for the density function
kde.grid <- function(data,H,axes=c("x","y"),alpha=0.001,res=1,dr=NULL)
{
  R <- get.telemetry(data,axes) # (times,dim)
  n <- nrow(R) # (times)
  
  H <- prepare.H(H,n) # (times,dim,dim)
  
  # how far to extend range from data as to ensure alpha significance in total probability
  z <- sqrt(-2*log(alpha))
  dH <- z * apply(H,1,function(h){sqrt(diag(h))}) # (dim,times)
  dH <- t(dH) # (times,dim)

  # now to find the necessary extent of our grid
  EXT <- rbind( apply(R-dH,2,min) , apply(R+dH,2,max) ) # (ext,dim)
  dEXT <- EXT[2,]-EXT[1,]
  
  # grid center
  mu <- apply(EXT,2,mean)
  
  if(is.null(dr)) { dr <- dEXT/res }

  res <- dEXT/dr

  # half resolution
  res <- ceiling(res/2+1)
  
  # grid locations
  R <- lapply(1:length(res),function(i){ (-res[i]:res[i])*dr[i] + mu[i] } ) # (grid,dim)
  names(R) <- axes
  
  grid <- list(R=R,dH=dH,dr=dr,EXT=EXT,dEXT=dEXT)
  return(grid)
}


##################################
# construct my own kde objects
# was using ks-package but it has some bugs
# alpha is the error goal in my total probability
kde <- function(data,H,axes=c("x","y"),bias=FALSE,W=rep(1,length(data$x)),alpha=0.001,res=NULL,dr=NULL,grid=NULL)
{
  r <- get.telemetry(data,axes)
  n <- nrow(r)
  
  # normalize weights
  W <- W/sum(W)
  
  # format bandwidth matrix
  H <- prepare.H(H,n)
  
  if(is.null(grid)) { grid <- kde.grid(data,H=H,axes=axes,alpha=alpha,res=res,dr=dr) }
  
  R <- grid$R
  # generalize this for future grid option use
  dH <- grid$dH
  dr <- grid$dr

  R0 <- sapply(R,first)
  # corner origin to minimize arithmetic later
  #r <- t(t(r) - R0)
  #R <- lapply(1:length(R0),function(i){ R[[i]] - R0[i] })
  
  cdf <- array(0,sapply(R,length))
  for(i in 1:n)
  {
    # sub-grid lower/upper bound indices
    i1 <- floor((r[i,]-dH[i,]-R0)/dr) + 1
    i2 <- ceiling((r[i,]+dH[i,]-R0)/dr) + 1
    
    SUB <- lapply(1:length(i1),function(d){ i1[d]:i2[d] })
    
    # I can't figure out how to do this in one line
    if(length(SUB)==1)
    { cdf[SUB[[1]]] <- cdf[SUB[[1]]] + W[i]*pnorm1(R[[1]][SUB[[1]]]-r[i,1],H[i,,],dr,alpha) }
    else if(length(SUB)==2)
    { cdf[SUB[[1]],SUB[[2]]] <- cdf[SUB[[1]],SUB[[2]]] + W[i]*pnorm2(R[[1]][SUB[[1]]]-r[i,1],R[[2]][SUB[[2]]]-r[i,2],H[i,,],dr,alpha) }
    else if(length(SUB)==3)
    { cdf[SUB[[1]],SUB[[2]],SUB[[3]]] <- cdf[SUB[[1]],SUB[[2]],SUB[[3]]] + W[i]*pnorm3(R[[1]][SUB[[1]]]-r[i,1],R[[2]][SUB[[2]]]-r[i,2],R[[3]][SUB[[3]]]-r[i,3],H[i,,],dr,alpha) }
  }

  dV <- prod(dr)
  if(!sum(bias)) { pdf <- cdf/dV }

  # cdf: cell probability -> probability included in contour
  cdf <- pdf2cdf(cdf,finish=FALSE)
  DIM <- cdf$DIM
  IND <- cdf$IND
  cdf <- cdf$cdf
  
  # just using minimum bias for now... not ideal in 3D
  vbias <- min(bias) # bias in variance along least biased axis
  vbias <- sqrt(vbias)^length(dr) # convert to volume bias
  # areas are biased to be estimated as debias*area
  if(vbias) 
  {
    # counting volume by dV
    VOL <- 1:length(cdf)

    # evaluate the debiased cdf on the original area grid
    cdf <- stats::approx(x=VOL/vbias,y=cdf,xout=VOL,yleft=0,yright=1)$y
    
    # recalculate pdf
    pdf <- diff(c(0,cdf))/dV
    pdf[IND] <- pdf
    pdf <- array(pdf,DIM)
  }
  # residual biases
  # bias <- bias/min(bias)
  
  cdf[IND] <- cdf # back in spatial order
  cdf <- array(cdf,DIM) # back in table form
  
  result <- list(PDF=pdf,CDF=cdf,r=R,dr=dr)
  
  return(result)
}

########################
# cdf: cell probability -> probability included in contour
pdf2cdf <- function(cdf,finish=TRUE)
{
  #cdf <- pdf * dV
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
pnorm2 <- function(X,Y,sigma,dr,alpha=0.001)
{
  dx <- dr[1]
  dy <- dr[2]
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


# This function is not ready for Kriging
pnorm3 <- function(X,Y,Z,sigma,dr,alpha=0.001)
{
  cdf <- prod(dr) * Gauss3(X,Y,Z,sigma)
  
  return(cdf)
}

# UNFINISHED
pnorm1 <- function(X,sigma,dr,alpha=0.001) { 0 }

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
# gaussian pdf
# assumes uncorrelated z-axis
Gauss3 <- function(X,Y,Z,sigma=NULL,sigma.inv=solve(sigma[1:2,1:2]),sigma.GM=sqrt(det(sigma[1:2,1:2])))
{
  cdf <- Gauss(X,Y,sigma=sigma[1:2,1:2],sigma.inv=sigma.inv,sigma.GM=sigma.GM)
  cdf <- cdf %o% (exp(-Z^2/(2*sigma[3,3]))/sqrt(2*pi*sigma[3,3]))
  
  return(cdf)
}


#####################
# AKDE CIs
CI.UD <- function(object,level.UD=0.95,level=0.95,P=FALSE)
{
  dV <- prod(object$dr)
  
  # point estimate
  if(is.na(level.UD))
  {
    # calculate mean area
    area <- sort(object$PDF,decreasing=TRUE,method="quick")
    area <- sum(area * (1:length(area))) * dV
  }
  else
  { area <- sum(object$CDF <= level.UD) * dV }
  
  # chi square approximation of uncertainty
  area <- chisq.ci(area,DOF=2*object$DOF.area,alpha=1-level)
  
  if(!P) { return(area) }
  
  # probabilities associated with these areas
  P <- round(area / dV)
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
  
  dx <- UD$dr[1]
  dy <- UD$dr[2]
  
  xmn <- UD$r$x[1]-dx/2
  xmx <- last(UD$r$x)+dx/2
  
  ymn <- UD$r$y[1]-dy/2
  ymx <- last(UD$r$y)+dy/2
  
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
    CL <- grDevices::contourLines(UD$r,z=UD$CDF,levels=P[i])
    
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

