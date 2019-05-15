# statistical mode
Mode <- function(x)
{
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
  mean(ux)
}


# assuming stationary model
modes.ctmm <- function(object,...)
{
  mu <- object$mu[1,]
  COV.mu <- object$COV.mu
  if(length(dim(COV.mu)==4)) { COV.mu <- COV.mu[,1,1,] }

  COV <- methods::getDataPart(object$sigma)
  density <- 1/sqrt(det(2*pi*COV))

  MODES <- list(list())
  MODES[[1]]$density <- density
  MODES[[1]]$mode <- mu
  MODES[[1]]$COV <- COV.mu

  return(MODES)
}


#
modes.UD <- function(object,...)
{
  CTMM <- object@CTMM
  COV.mu <- object$COV.mu
  if(length(dim(COV.mu)==4)) { COV.mu <- COV.mu[,1,1,] }
  # best/worst output DOF
  DOF.max <- DOF.mean(CTMM)
  DOF.min <- 0
  # best/worst input densities
  H <- object$H
  PDF.min <- 1/sqrt(det(2*pi*H))/object$DOF.H
  PDF.max <- 1/sqrt(det.covm(scale.covm(CTMM$sigma,2*pi)))

  PDF <- object$PDF
  dx <- object$dr$x
  dy <- object$dr$y

  MODES <- list()
  for(i in 2:(nrow(PDF)-1))
  {
    for(j in 2:(ncol(PDF)-1))
    {
      if(PDF[i,j]==max(PDF[i+(-1):1,j+(-1):1])) # local maximum
      {
        MODES[[length(MODES)+1]] <- list()

        SUB <- log(PDF[i+(-1):1,j+(-1):1])

        # calculate gradient
        GRAD <- numeric(2)
        GRAD[1] <- mean( diff(SUB[2,1:2]) , diff(SUB[2,2:3]) )/dx
        GRAD[2] <- mean( diff(SUB[1:2,2]) , diff(SUB[2:3,2]) )/dy

        # calculate Hessian
        HESS <- diag(2)
        HESS[1,1] <- diff(diff(SUB[2,]))/dx^2 # (+x,+x)
        HESS[2,2] <- diff(diff(SUB[,2]))/dy^2 # (+y,+y)
        HESS[1,2] <- diff(diff(c(SUB[1,1],SUB[2,2],SUB[3,3])))/(dx^2+dy^2) # (+x+y,+x+y)
        HESS[2,1] <- diff(diff(c(SUB[1,3],SUB[2,2],SUB[3,1])))/(dx^2+dy^2) # (+x-y,+x-y)
        HESS[1,2] <- HESS[2,1] <- (HESS[1,2]-HESS[2,1])/2 # (+x,+y)
        COV <- -PDsolve(HESS)

        # locate mode
        r0 <- c(object$r$x[i],object$r$y[j])
        mode <- r0 + COV %*% GRAD
        # in case of bad Hessian
        mode[1] <- clamp(mode[1],object$r$x[i-1],object$r$x[i+1])
        mode[2] <- clamp(mode[2],object$r$x[j-1],object$r$x[j+1])
        names(mode) <- c('x','y')
        MODES[[length(MODES)]]$mode <- mode

        # estimate mode density
        dr <- mode - r0
        MODES[[length(MODES)]]$density <- PDF[i,j] + c(GRAD %*% dr) + (1/2)*c(dr %*% HESS %*% dr)

        # calculate weight from density
        MODES[[length(MODES)]]$weight <- MODES[[length(MODES)]]$density * (2*pi/sqrt(det(HESS)))

        # calculate COV.mode





      }
    }
  }

  # sort modes

  return(MODES)
}
