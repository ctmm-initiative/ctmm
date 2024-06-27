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
mesh3d.UD <- function(UD,level.UD=0.95)
{
  DIM <- dim(UD$CDF)
  IN <- UD$CDF < level.UD
  isVERT <- array(FALSE,DIM)
  VERTEX <- array(0,c(0,3)) # matrix of vertex (x,y,z) per index
  NORM <- array(0,c(0,3)) # matrix of normal vectors per index
  TRI <- array(integer(1),c(0,3)) # matrix of 3 vertex indices per triangle

  S <- -1:1
  C <- 3*3+5 # center index of cube

  SHELL <- array(integer(1),c(0,3))
  for(i in S) { for(j in S) { for(k in S) { SHELL <- rbind(SHELL,c(i,j,k)) } } }
  SHELL <- SHELL[-14,] # (0,0,0)

  # table intersection
  row.intersect <- function(X,Y)
  {
    Z <- rbind(X,Y)
    # sort by all columns
    Z <- Z[order(Z[,1],Z[,2],Z[,3]),]
    # do rows change
    I <- apply(Z,2,diff) # [n-1,3]
    I <- apply(abs(I),1,sum) # [n-1]
    I <- c(TRUE,I) # [n]
    # return repeated rows
    Z[!I,]
  }

  # table unique
  row.unique <- function(Z)
  {
    # sort by all columns
    Z <- Z[order(Z[,1],Z[,2],Z[,3]),]
    # do rows change
    I <- apply(Z,2,diff) # [n-1,3]
    I <- apply(abs(I),1,sum) # [n-1]
    I <- c(TRUE,I) # [n]
    # return unique rows  (different from previous)
    Z[I,]
  }

  # level.fn <- function(adbmal,return.X=FALSE)
  # {
  #   lambda <- 1/adbmal
  #   M <- 2/lambda*diag(3) - VV
  #   X <- V[1:3]
  #
  #   X <- t(M) %*% X
  #   M <- t(M) %*% M
  #   M <- pd.solve(M)
  #   X <- M %*% X # this is now (x,y,z) for lambda
  #
  #   if(return.X) { return(X) }
  #
  #   level <- SUB[C] + sum(X[1:3]*V[1:3]) + X[1]*X[2]*V[4] + X[2]*X[3]*V[5] + X[3]*X[1]*V[6] + sum(X[1:3]^2*V[7:9])
  #   return((level-level.UD)^2)
  # }

  # vertex loop
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
          # fit a solvable quadratic regression function to cube data
          x <- y <- z <- array(-1:1,c(3,3,3))
          y <- aperm(x,c(3,1,2))
          z <- aperm(x,c(2,3,1))

          # M <- matrix(0,3^3-1,9)
          M <- matrix(0,3^3-1,3)
          # linear regression terms
          M[,1] <- c(x)[-C]
          M[,2] <- c(y)[-C]
          M[,3] <- c(z)[-C]
          # # terms also found in tri-linear splines
          # M[,4] <- c(x*y)[-C]
          # M[,5] <- c(y*z)[-C]
          # M[,6] <- c(z*x)[-C]
          # # quadratic regression terms
          # M[,7] <- c(x*x)[-C]
          # M[,8] <- c(y*y)[-C]
          # M[,9] <- c(z*z)[-C]

          SUB <- UD$CDF[i+S,j+S,k+S]
          V <- c(SUB)[-C] - SUB[C]

          V <- t(M) %*% V
          M <- t(M) %*% M
          M <- pd.solve(M)
          V <- M %*% V # coefficients

          # minimize distance from center while constrained to surface
          # L == x^2+y^2+z^2 + lambda*(level - V*(x,y,z,xy,yz,zx,xx,yy,zz) )
          # minimized by
          # 2*x == lambda*( V['x'] + V['xy']*y + V['zx']*z + 2*V['xx'] )
          # 2*y == lambda*( V['y'] + V['yz']*z + V['xy']*x + 2*V['yy'] )
          # 2*z == lambda*( V['z'] + V['zx']*x + V['yz']*y + 2*V['zz'] )
          #
          # (2/lambda*I - VV)*(x,y,z) = V['x','y','z']
          # VV <- matrix(0,3,3)
          # VV[1,2] <- VV[2,1] <- V[4]
          # VV[2,3] <- VV[3,2] <- V[5]
          # VV[3,1] <- VV[1,3] <- V[6]
          # diag(VV) <- 2*V[7:9]

          # linear regression solution lambda/2*V[1:3]^2 == level.UD-SUB[C]
          lambda <- 2*(level.UD-SUB[C])/sum(V[1:3]^2) # initial guess (linear regression solution)
          # adbmal <- 1/lambda # doesn't diverge between -Inf and +Inf
          # adbmal <- optimizer(adbmal,level.fn)$par
          # lambda <- 1/adbmal
          # X <- level.fn(adbmal,return.X=TRUE)
          X <- lambda/2*V[1:3]

          if(all(abs(X)<1/2)) # record if within pixel
          {
            # store vertex coordinate
            Xc <- t(X + c(i,j,k))
            VERTEX <- rbind(VERTEX,Xc)
            isVERT[i,j,k] <- length(VERTEX) # store index for triangulation

            # store normal vector
            # Xn <- V[1:3] + c(X[2]*V[4]+X[3]*V[6],X[1]*V[4]+X[3]*V[5],X[2]*V[5]+X[1]*V[6]) + 2*X[1:3]*V[7:9]
            Xn <- V[1:3]
            Xn <- +lambda/2*t(Xn) # gradient points outwards to higher coverage
            # no point in normalizing this until we include (x,y,z) scale information
            NORM <- rbind(NORM,Xn)
          }
        } # end vertex solver
      } # end z loop
    } # end y loop
    utils::setTxtProgressBar(pb,i)
  } # end x loop

  # Nearest Neighbor triangles
  utils::setTxtProgressBar(pb,0)
  for(i in 2:(DIM[1]-1))
  {
    for(j in 2:(DIM[2]-1))
    {
      for(k in 2:(DIM[3]-1))
      {
        i1 <- isVERT[i,j,k] # first vertex of triangle
        if(i1)
        {
          V2 <- t(c(i,j,k) + t(SHELL))
          for(l in 1%:%nrow(V2))
          {
            i2 <- isVERT[V2[l,1],V2[l,2],V2[l,3]] # second vertex of triangle
            if(i2)
            {
              V3 <- t(V2[l,] + t(SHELL))
              # must also be a neighbor of v1
              V3 <- row.intersect(V3,V2)
              for(m in 1%:%nrow(V3))
              {
                i3 <- isVERT[V3[m,1],V3[m,2],V3[m,3]]
                if(i3){ TRI <- rbind(TRI,c(i1,i2,i3)) }
              }
            }
          }
        }
      } # end z loop
    } # end y loop
    utils::setTxtProgressBar(pb,i)
  } # end x loop

  # remove duplicate NN triangles
  TRI <- row.unique(TRI)

  # convert from fractional indices to locations
  colnames(VERTEX) <- c('x','y','z')
  VERTEX[,'x'] <- UD$r$x[1] + (VERTEX[,'x']-1)*UD$dr['x']
  VERTEX[,'y'] <- UD$r$y[1] + (VERTEX[,'y']-1)*UD$dr['y']
  VERTEX[,'z'] <- UD$r$z[1] + (VERTEX[,'z']-1)*UD$dr['z']

  # convert from index gradient to spatial gradient
  colnames(NORM) <- c('x','y','z')
  NORM[,'x'] <- NORM[,'x']*UD$dr['x']
  NORM[,'y'] <- NORM[,'y']*UD$dr['y']
  NORM[,'z'] <- NORM[,'z']*UD$dr['z']

#  MESH <- rgl::mesh3d(VERTEX,triangles=TRI)
#  return(MESH)

  close(pb)
  return(VERTEX)
}

