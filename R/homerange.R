###################
# homerange wrapper function
homerange <- function(data=NULL,CTMM=NULL,method="AKDE",...)
{
  method <- match.arg(method,c("AKDE","AGDE"))

  if(class(data)[1]=="ctmm")
  {
    TEMP <- CTMM
    CTMM <- data
    data <- TEMP
  }

  # these are the only possibilities so far
  if(is.null(data))
  {
    # if(method=="AKDE") {  warning('No data provided. Method changed to "AGDE".') }
    method <- "AGDE"
  }
  # else
  # { method <- "AKDE" }

  if(method=="AKDE")
  { RETURN <- akde(data,CTMM,...) }
  else if(method=="AGDE")
  { RETURN <- agde(CTMM,...) }
  return(RETURN)
}


agde <- function(data=NULL,CTMM=NULL,R=list(),variable="utilization",error=0.001,res=100,grid=NULL,...)
{
  if(class(data)[1]=="ctmm")
  {
    TEMP <- CTMM
    CTMM <- data
    data <- TEMP
  }

  axes <- CTMM$axes
  grid <- format_grid(grid,axes=axes)

  level.UD <- CTMM$level.UD
  if(is.null(level.UD))
  {
    data <- data.frame(CTMM$mu[1,,drop=FALSE])
    H <- methods::getDataPart(CTMM$sigma)
    dr <- sqrt(diag(H))/res
    dim(H) <- c(1,2,2)
    EXT <- extent(CTMM,level=1-error)[,axes] # Gaussian extent (includes uncertainty)
    grid <- kde.grid(data,H=H,axes=axes,alpha=error,res=res,dr=dr,grid=grid,EXT.min=EXT)

    pdf <- kde(data,H=H,RASTER=R,axes=axes,CTMM=CTMM,res=res,alpha=error,grid=grid)
    pdf <- new.UD(pdf,info=attr(CTMM,"info"),type='range',variable="utilization",CTMM=CTMM)

    pdf$DOF.area <- array( DOF.area(CTMM) , 2)
    names(pdf$DOF.area) <- CTMM$axes
  }
  else # traditional RSF with fixed boundary
  {
    proj <- projection(CTMM)

    level.UD <- level.UD@Polygons[[1]]@coords

    x <- range(level.UD[,1])
    y <- range(level.UD[,2])
    EXT <- cbind(x,y)
    rownames(EXT) <- c('min','max')
    data <- data.frame(x=mean(x),y=mean(y))
    x <- diff(x)
    y <- diff(y)
    dr <- c(x=x/res,y=y/res)
    H <- diag(c(x,y))/qmvnorm(1-error,length(axes)) # cancels out in kde.grid()
    H <- H^2
    dim(H) <- c(1,2,2)

    grid <- kde.grid(data,H=H,axes=axes,alpha=error,res=res,dr=dr,grid=grid,EXT.min=EXT)
    r <- grid$r
    dr <- grid$dr
    dV <- prod(dr)

    DIM <- c(length(r$x),length(r$y))
    xy <- array(0,c(DIM,2))
    xy[,,1] <- r$x
    xy <- aperm(xy,c(2,1,3))
    xy[,,2] <- r$y
    xy <- aperm(xy,c(2,1,3))
    dim(xy) <- c(prod(DIM),2)

    if(length(R))
    {
      R <- expand.factors(R,CTMM$formula,fixed=TRUE)

      # calculate RASTERs on spatial grid
      R <- lapply(R,function(RAS){R.grid(r,proj=proj,RAS)})
      # this needs to be moved up for multiple individuals?
    }

    # not finished with this part
    STATIONARY <- is.stationary(CTMM,R)

    # suitability raster
    if(length(R))
    { R <- R.suit(R,CTMM) }

    PMF <- sp::point.in.polygon(xy[,1],xy[,2],level.UD[,1],level.UD[,2])
    dim(PMF) <- DIM

    if(length(R)) { PMF <- PMF*R }
    PMF <- PMF/sum(PMF)
    CDF <- pmf2cdf(PMF)

    pdf <- list(PDF=PMF/dV,CDF=CDF,r=r,dr=dr,axes=axes)
    pdf <- new.UD(pdf,info=attr(CTMM,"info"),type='range',variable="utilization",CTMM=CTMM)
  }

  return(pdf)
}

