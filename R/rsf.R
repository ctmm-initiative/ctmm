rsf.fit <- function(data,UD,R,...)
{
  CTMM <- UD@CTMM
  W <- UD$weights * UD$DOF.area

  R <- listify(R)
  r <- NULL
  for(i in 1:length(R))
  {
    r <- cbind(r,...)
  }
}
