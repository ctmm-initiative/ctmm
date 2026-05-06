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
  data <- listify(data)
  CTMM <- listify(CTMM)

  if(class(data[[1]])[1]=="ctmm")
  {
    TEMP <- CTMM
    CTMM <- data
    data <- TEMP
    rm(TEMP)
  }
  NAMES <- names(data)

  n <- length(CTMM)
  axes <- CTMM[[1]]$axes
  grid <- format_grid(grid,axes=axes)
  level.UD <- sapply(CTMM,function(M){M$level.UD}) # conventional RSF polygon

  H <- array(0,c(n,2,2))
  pdf <- list()

  if(is.null(level.UD[[1]])) # AGDE or AGDE RSF
  {
    # bandwidth loop
    for(i in 1:n)
    {
      data[[i]] <- data.frame(CTMM[[i]]$mu[1,,drop=FALSE])
      H[i,,] <- methods::getDataPart(CTMM[[i]]$sigma)
    }
    data <- do.call(rbind,data)

    ## grid construction
    EXT <- extent(CTMM,level=1-error)[,axes] # Gaussian extent (includes uncertainty)
    grid <- kde.grid(data,H=H,axes=axes,alpha=error,res=res,grid=grid,EXT.min=EXT)
    dr <- grid$dr

    # kde loop
    for(i in 1:n)
    {
      pdf[[i]] <- kde(data[i,],H=H[i,,],RASTER=R,axes=axes,CTMM=CTMM[[i]],res=res,alpha=error,grid=grid)
      pdf[[i]] <- new.UD(pdf[[i]],info=attr(CTMM[[i]],"info"),type='range',variable="utilization",CTMM=CTMM[[i]])

      pdf[[i]]$DOF.area <- array( DOF.area(CTMM[[i]]) , 2)
      names(pdf[[i]]$DOF.area) <- CTMM[[i]]$axes
    }
  }
  else
  {
    proj <- projection(CTMM)
    level.UD <- lapply(level.UD,function(L){L@Polygons[[1]]@coords})
    level.UD <- do.call(rbind,level.UD)

    # bandwidth loop
    for(i in 1:n)
    {
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
    }

    grid <- kde.grid(data,H=H,axes=axes,alpha=error,res=res,dr=dr,grid=grid,EXT.min=EXT)
    r <- grid$r
    dr <- grid$dr
    dV <- prod(dr)

    DIM <- c(length(r$x),length(r$y))

    # kde loop
    for(i in 1:n)
    {
      xy <- array(0,c(DIM,2))
      xy[,,1] <- r$x
      xy <- aperm(xy,c(2,1,3))
      xy[,,2] <- r$y
      xy <- aperm(xy,c(2,1,3))
      dim(xy) <- c(prod(DIM),2)

      if(length(R))
      {
        Ri <- expand.factors(R,CTMM[[i]]$formula,fixed=TRUE)$R

        # calculate RASTERs on spatial grid
        Ri <- lapply(Ri,function(RAS){R.grid(r,proj=proj,RAS)})
        # this needs to be moved up for multiple individuals?
      }

      # not finished with this part
      STATIONARY <- is.stationary(CTMM[[i]],Ri)

      # suitability raster
      if(length(Ri))
      { Ri <- R.suit(Ri,CTMM[[i]]) }

      PMF <- sp::point.in.polygon(xy[,1],xy[,2],level.UD[,1],level.UD[,2])
      dim(PMF) <- DIM

      if(length(Ri)) { PMF <- PMF*Ri }
      PMF <- PMF/sum(PMF)
      CDF <- pmf2cdf(PMF)

      pdf[[i]] <- list(PDF=PMF/dV,CDF=CDF,r=r,dr=dr,axes=axes)
      pdf[[i]] <- new.UD(pdf,info=attr(CTMM[[i]],"info"),type='range',variable="utilization",CTMM=CTMM[[i]])
    }
  }

  if(length(pdf)==1)
  { return(pdf[[1]]) }
  else
  {
    names(pdf) <- NAMES
    return(pdf)
  }
}

