# mag.numeric <- function(x)
# {
#   x <- cbind(x)
#   x <- norm(x,type="F")
#   return(x)
# }
#
# mag.complex <- function(x)
# {
#   x <- cbind(x)
#   x <- Adj(x) %*% x
#   x <- diag(x,nrow(x))
#   x <- sum(x)
#   x <- sqrt(x)
#   return(x)
# }

# magnitude (for residuals)
mag.telemetry <- function(x,axes=c("x","y"),...)
{
  x <- get.telemetry(x,axes=axes)
  x <- rowSums(x^2)
  x <- sqrt(x)
  return(x)
}

