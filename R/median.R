# convert lat, long, alt data.frame to x,y,z (not projected) matrix
ellipsoid2cartesian <- function(data)
{
  # assume Movebank data.frame
  lon <- data$longitude * (2*pi/360)
  lat <- data$latitude * (2*pi/360)
  s <- data$z # distance from surface
  if(is.null(s)) { s <- 0 }

  # convert to x,y,z coordinates (3D)
  z <- (DATA.EARTH$R.PL + s) * sin(lat)
  r <- (DATA.EARTH$R.EQ + s) * cos(lat)

  x <- r * cos(lon)
  y <- r * sin(lon)

  data <- cbind(x,y,z)
  return(data)
}


# convert back
cartesian2ellipsoid <- function(mu)
{
  mu <- rbind(mu)

  x <- mu[,1]
  y <- mu[,2]
  z <- mu[,3]

  lon <- atan2(y,x)

  r <- sqrt(x^2 + y^2)

  r <- r/DATA.EARTH$R.EQ
  z <- z/DATA.EARTH$R.PL

  lat <- atan2(z,r)

  mu <- cbind(longitude=lon,latitude=lat)  * (360/(2*pi))

  return(mu)
}


# median of dataset
median.telemetry <- function(x,na.rm=FALSE,...)
{
  if(class(x)[1]=="list")
  {
    if(length(x)>1)
    {
      id <- mean_info(x)$identity
      # is there a common projection to preserve
      proj <- sapply(x,projection.telemetry)
      proj <- unlist(proj) # need for multiple NULLs

      if(length(proj)<length(x))
      { proj <- NULL }
      else
      {
        proj <- unique(proj)
        if(length(proj)>1) { proj <- NULL }
      }

      x <- lapply(x,function(d){ median.telemetry(d,...) }) # this will only ever return long-lat & x-y columns
      x <- do.call(rbind,x)
      x <- as.data.frame(x)
      x <- new.telemetry(x,info=list(identity=id,projection=proj),UERE=new.UERE())
    }
    else
    { x <- x[[1]] }
    # now flow through to per-individual analysis below
  }

  id <- paste0("median of ",attr(x,'info')$identity)
  # tz <- attr(x,'info')$timezone

  # t <- stats::median(x$t)
  # timestamp <- as.character(as.POSIXct(t,tz=tz,origin="1970/01/01"))

  if(all(c("longitude","latitude") %in% names(x)))
  {
    proj <- projection(x)

    mu <- median_longlat(x,...)
    x <- data.frame(longitude=mu[,"longitude"],latitude=mu[,"latitude"])
    x <- new.telemetry(x,info=list(identity=id,projection=proj),UERE=new.UERE())

    projection(x) <- proj # NULL is handled
  }
  else # fallback for turtle data (k>1 not supported)
  {
    x <- cbind(x$x,x$y)
    mu <- apply(x,2,stats::median)
    mu <- c(Gmedian::Gmedian(x,init=mu,...))

    x <- data.frame(x=mu[1],y=mu[2])
    x <- new.telemetry(x,info=list(identity=id),UERE=new.UERE())
  }

  return(x)
}


# ellipsoidal median (k-cluster)
median_longlat <- function(data,k=1,...)
{
  data <- ellipsoid2cartesian(data)

  # what is the centroid of the data in 3D
  mu <- apply(data,2,stats::median)
  STAT <- Gmedian::GmedianCov(data,init=mu,scores=0,nstart=10,...)
  mu <- c(STAT$median)
  COV <- STAT$covmedian

  if(k==1 || nrow(data)==1)
  {
    mu <- cartesian2ellipsoid(mu)
    return(mu)
  }

  # k==2 below
  if(nrow(data)==1)
  { mu <- rbind(data,data) }
  else if(nrow(data)==2)
  { mu <- data }
  else
  {
    # find the longest axis of variation
    COV <- eigen(COV)$vectors[,1]
    # distances along long axis
    COV <- t(t(data)-mu) %*% COV
    # centroid along long axis
    mu <- stats::median(COV)
    # detrended distances along long axis
    COV <- c(COV) - mu

    # positive half
    SUB <- data[COV>=0,,drop=FALSE]
    if(nrow(SUB)==1) { mu1 <- SUB }
    else
    {
      mu1 <- apply(SUB,2,stats::median)
      mu1 <- c(Gmedian::Gmedian(SUB,init=mu1))
    }

    # negative half
    SUB <- data[COV<=0,,drop=FALSE]
    if(nrow(SUB)==1) { mu2 <- SUB }
    else
    {
      mu2 <- apply(SUB,2,stats::median)
      mu2 <- c(Gmedian::Gmedian(SUB,init=mu2))
    }

    # 2-mode cluster
    mu <- rbind(mu1,mu2)
    if(nrow(mu)>3) { mu <- Gmedian::kGmedian(data,ncenters=mu,nstart=10,...)$centers }

    # make sure separation is not reduced by clustering
    if(sum((mu1-mu2)^2)>sum((mu[1,]-mu[2,])^2))
    {
      mu1 -> mu[1,]
      mu2 -> mu[2,]
    }
  }
  mu <- cartesian2ellipsoid(mu)

  # order from west to east
  if(mu[1,1] > mu[2,1])
  {
    mu[1,] -> mu2
    mu[2,] -> mu1

    mu <- rbind(mu1,mu2)
  }

  colnames(mu) <- c("longitude","latitude")
  return(mu)
}
