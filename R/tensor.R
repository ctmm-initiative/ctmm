# tensor (inner) product
'%.%' <- function(x,y)
{
  DIMx <- dim(x)
  DIMy <- dim(y)

  dx <- length(DIMx)
  dy <- 1

  x <- array(x,c(prod(DIMx[-dx]),DIMx[dx]))
  y <- array(y,c(DIMy[dy],prod(DIMy[-dy])))
  xy <- x %*% y

  DIM <- c(DIMx[-dx],DIMy[-dy])
  if(length(DIM)>2) { xy <- array(xy,DIM) }

  return(xy)
}

# contract wrt two indices
'%..%' <- function(x,y)
{
  DIMx <- dim(x)
  DIMy <- dim(y)

  dx <- length(DIMx) - (0:1)
  dy <- 1:2

  x <- array(x,c(prod(DIMx[-dx]),prod(DIMx[dx])))
  y <- array(y,c(prod(DIMy[dy]),prod(DIMy[-dy])))
  xy <- x %*% y

  DIM <- c(DIMx[-dx],DIMy[-dy])
  if(length(DIM)>2) { xy <- array(xy,DIM) }

  return(xy)
}
