# x[i,*] <- value
# where * includes the remaining indices
Assign <- function(x,i,value,index=1)
{
  #
}

arrayify <- function(x)
{
  if(is.null(dim(x))) { x <- array(x,length(x)) }
  return(x)
}

# first-index bind
fbind <- function(x,y)
{
  x <- arrayify(x)
  y <- arrayify(y)

  DIMx <- dim(x)
  dim(x) <- c(DIMx[1],prod(DIMx[-1]))
  DIMy <- dim(y)
  dim(y) <- c(DIMy[1],prod(DIMy[-1]))
  x <- rbind(x,y)
  dim(x) <- c(DIMx[1]+DIMy[1],DIMx[-1])
  return(x)
}

# last-index bind
lbind <- function(x,y)
{
  x <- arrayify(x)
  y <- arrayify(y)

  DIMx <- dim(x)
  n <- length(DIMx)
  dim(x) <- c(prod(DIMx[-n]),DIMx[n])
  DIMy <- dim(y)
  dim(y) <- c(prod(DIMy[-n]),DIMy[n])
  x <- cbind(x,y)
  dim(x) <- c(DIMx[-n],DIMx[n]+DIMy[n])
  return(x)
}

# tensor (inner) product
'%.%' <- function(x,y)
{
  x <- arrayify(x)
  y <- arrayify(y)

  DIMx <- dim(x)
  DIMy <- dim(y)

  dx <- length(DIMx)
  dy <- 1

  dim(x) <- c(prod(DIMx[-dx]),DIMx[dx])
  dim(y) <- c(DIMy[dy],prod(DIMy[-dy]))
  xy <- x %*% y

  DIM <- c(DIMx[-dx],DIMy[-dy])
  if(length(DIM)>1) { dim(xy) <- DIM }

  return(xy)
}

# contract w.r.t. two indices
'%..%' <- function(x,y)
{
  x <- arrayify(x)
  y <- arrayify(y)

  DIMx <- dim(x)
  DIMy <- dim(y)

  dx <- length(DIMx) - (1:0)
  dy <- 1:2

  dim(x) <- c(prod(DIMx[-dx]),prod(DIMx[dx]))
  dim(y) <- c(prod(DIMy[dy]),prod(DIMy[-dy]))
  xy <- x %*% y

  DIM <- c(DIMx[-dx],DIMy[-dy])
  if(length(DIM)>1) { dim(xy) <- DIM }

  return(xy)
}
