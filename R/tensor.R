# x[i,*] <- value
# where * includes the remaining indices
Assign <- function(x,i,value,index=1)
{
  #
}


# tensor (inner) product
'%.%' <- function(x,y)
{
  DIMx <- dim(x)
  DIMy <- dim(y)

  dx <- length(DIMx)
  dy <- 1

  dim(x) <- c(prod(DIMx[-dx]),DIMx[dx])
  dim(y) <- c(DIMy[dy],prod(DIMy[-dy]))
  xy <- x %*% y

  DIM <- c(DIMx[-dx],DIMy[-dy])
  if(length(DIM)>2) { dim(xy) <- DIM }

  return(xy)
}

# contract wrt two indices
'%..%' <- function(x,y)
{
  DIMx <- dim(x)
  DIMy <- dim(y)

  dx <- length(DIMx) - (1:0)
  dy <- 1:2

  dim(x) <- c(prod(DIMx[-dx]),prod(DIMx[dx]))
  dim(y) <- c(prod(DIMy[dy]),prod(DIMy[-dy]))
  xy <- x %*% y

  DIM <- c(DIMx[-dx],DIMy[-dy])
  if(length(DIM)>2) { dim(xy) <- DIM }

  return(xy)
}
