# plot overhead 3D UD
plot3d <- function(data=NULL,UD,level=0.95,level.UD=0.95,xlim=NULL,ylim=NULL,ext=NULL,...)
{
  UD <- listify(UD)
  z <- NULL; rm(z) # R check bug
  Z <- dZ <- zext <- MODE <- list()
  if(is.null(ext)) { ext <- list() }
  for(i in 1:length(UD))
  {
    Z[[i]] <- UD[[i]]$r$z
    dZ[[i]] <- UD[[i]]$dr['z']

    # fix 2D extent
    if(!length(ext) && (is.null(xlim) || is.null(ylim)))
    {
      FLAT <- UD[[i]]
      FLAT$CDF <- apply(FLAT$PDF,1:2,sum) * prod(FLAT$dr) # marginal PMF
      FLAT$CDF <- pmf2cdf(FLAT$CDF)
      ext[[i]] <- extent(FLAT,level=level,level.UD=level.UD)
      rm(FLAT)
    }

    # fix 3D extent
    zext[[i]] <- UD[[i]]$CDF <= level.UD
    zext[[i]] <- apply(zext[[i]],3,any)
    zext[[i]] <- Z[[i]][zext[[i]]]
    zext[[i]] <- c(zext[[i]][1],last(zext[[i]]))

    # z mode to start at
    MODE[[i]] <- apply(UD[[i]]$PDF,3,sum)
    MODE[[i]] <- which.max(MODE[[i]])
    MODE[[i]] <- Z[[i]][MODE[[i]]]
  }

  EXT <- extent(ext)
  zext <- range(unlist(zext))
  MODE <- mean(unlist(MODE))

  plot.fn <- function(z)
  {
    BAD <- rep(FALSE,length(UD))
    for(i in 1:length(UD))
    {
      # if(z<ext[[i]][1] || ext[[i]][2]<z)
      # { BAD[i] <- TRUE }
      # else
      # {
      j <- round((z-Z[[i]][1])/dZ[[i]]) + 1
      UD[[i]]$PDF <- UD[[i]]$PDF[,,j]
      UD[[i]]$CDF <- UD[[i]]$CDF[,,j]
      UD[[i]]$H <- UD[[i]]$H[1:2,1:2]
      UD[[i]]$h <- UD[[i]]$h[1:2]
      UD[[i]]$bias <- UD[[i]]$bias[1:2]
      UD[[i]]$DOF.area <- UD[[i]]$DOF.area[1:2]
      UD[[i]]$axes <- UD[[i]]$axes[1:2]
      UD[[i]]$dr <- UD[[i]]$dr[1:2]
      # }
    } # end for
    # UD <- UD[!BAD] # other arguments would be modified too...
    if(is.null(data))
    { plot.UD(UD,level=level,level.UD=level.UD,xlim=xlim,ylim=ylim,ext=EXT,...) }
    else
    { plot.telemetry(data,UD=UD,level=level,level.UD=level.UD,xlim=xlim,ylim=ylim,ext=EXT,...) }
  }
  manipulate::manipulate(plot.fn(z),z=manipulate::slider(min=zext[1],max=zext[2],initial=MODE,label="z (m)"))
}

# calculate vertex mesh of 3D UD
# UNFINISHED
vertex <- function(UD,level.UD=0.95)
{
  VERTEX <- array(0,c(0,3))

  DIM <- dim(UD$CDF)
  IN <- UD$CDF < level.UD

  S <- -1:1
  C <- 3*3+5 # center index of cube

  pb <- utils::txtProgressBar(min=1,max=DIM[1]-1,initial=1,style=3)
  for(i in 2:(DIM[1]-1))
  {
    for(j in 2:(DIM[2]-1))
    {
      for(k in 2:(DIM[3]-1))
      {
        SIN <- IN[i+S,j+S,k+S]
        # is there a surface running through pixel (i,j,k) ?
        if(any(SIN[C] != SIN[-C]))
        {
          # fit a solvable regression function to cube data
          x <- y <- z <- array(-1:1,c(3,3,3))
          y <- aperm(x,c(3,1,2))
          z <- aperm(x,c(2,3,1))

          M <- matrix(0,3^3-1,6)
          M[,1] <- c(x)[-C]
          M[,2] <- c(y)[-C]
          M[,3] <- c(z)[-C]
          M[,4] <- c(x*y)[-C]
          M[,5] <- c(y*z)[-C]
          M[,6] <- c(z*x)[-C]

          SUB <- UD$CDF[i+S,j+S,k+S]
          V <- c(SUB)[-C] - SUB[C]

          V <- t(M) %*% V
          M <- t(M) %*% M
          M <- PDsolve(M)
          V <- M %*% V # coefficients

          # minimize distance from center while constrained to surface
          # L == x^2+y^2+z^2 + lambda*(level - V*(x,y,z,xy,yz,zx) )
          # minimized by
          # 2*x == lambda*( V['x'] + V['xy']*y + V['zx']*z)
          # 2*y == lambda*( V['y'] + V['yz']*z + V['xy']*x)
          # 2*z == lambda*( V['z'] + V['zx']*x + V['yz']*y)
          #
          # (2/lambda*I - VV)*(x,y,z) = V['x','y','z']
          VV <- matrix(0,3,3)
          VV[1,2] <- VV[2,1] <- V[4]
          VV[2,3] <- VV[3,2] <- V[5]
          VV[3,1] <- VV[1,3] <- V[6]

          level.fn <- function(lambda,return.X=FALSE)
          {
            if(lambda<=0) { return(Inf) }
            M <- 2/lambda*diag(3) - VV
            X <- V[1:3]

            X <- t(M) %*% X
            M <- t(M) %*% M
            M <- PDsolve(M)
            X <- M %*% X # this is now (x,y,z) for lambda

            if(return.X) { return(X) }

            level <- X[1]*V[1] + X[2]*V[2] + X[3]*V[3] + X[1]*X[2]*V[4] + X[2]*X[3]*V[5] + X[3]*X[1]*V[6]
            return((level-level.UD)^2)
          }

          LEVEL <- max(2/sum(abs(V)),1) # initial guess (order of magnitude)
          LEVEL <- optimizer(LEVEL,level.fn,lower=0)$par
          LEVEL <- level.fn(LEVEL,return.X=TRUE)

          if(all(abs(LEVEL)<1/2)) # record if within pixel
          {
            LEVEL <- LEVEL + c(i,j,k)
            LEVEL <- t(LEVEL)
            VERTEX <- rbind(VERTEX,LEVEL)
          }

          # NOT SURE IF THIS IS WORTH DOING WITHOUT GRADIENTS EXPLICITLY CALCULATED
          # TODO normal vector (gradient)
          # TODO normal vector (gradient)
          # TODO normal vector (gradient)
          # TODO normal vector (gradient)
        } # end vertex solver
      } # end z loop
    } # end y loop
    utils::setTxtProgressBar(pb,i)
  } # end x loop

  colnames(VERTEX) <- c('x','y','z')
  # convert from indices to locations
  VERTEX[,'x'] <- UD$r$x[1] + (VERTEX[,'x']-1)*UD$dr['x']
  VERTEX[,'y'] <- UD$r$y[1] + (VERTEX[,'y']-1)*UD$dr['y']
  VERTEX[,'z'] <- UD$r$z[1] + (VERTEX[,'z']-1)*UD$dr['z']

  close(pb)
  return(VERTEX)
}

