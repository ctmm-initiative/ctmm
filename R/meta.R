mean.ctmm <- function(x,summary="Wishart",prior="Inverse-Wishart",method="exact",...)
{
  summary <- match.arg(summary,c("Wishart","exp-log-chisq","log-normal"))
  prior <- match.arg(prior,c("Inverse-Wishart","log-normal"))
  method <- match.arg(method,c("exact","Laplace","MCMC"))

  AXES <- length(x[[1]]$axes)

  # analyticlly solvable
  if(summary=="Wishart" && prior=="Inverse-Wishart")
  {
    SIGMA <- list()
    DOF <- list()

    for(i in 1:length(x))
    {
      # extract covariance matrix SIGMA and covariance of covariance matrix COV
      SIGMA[[i]] <- x[[i]]$sigma # Wishart matrix / n
      if(x[[i]]$isotropic) # chi^2
      { PAR <- 'major' }
      else # approximate Wishart DOF (exact if Wishart)
      { PAR <- c('major','minor') }
      EST <- SIGMA[[i]]@par[PAR]
      DOF[[i]] <- (2/AXES) * c(EST %*% PDsolve(x[[i]]$COV[PAR,PAR]) %*% EST) # average multiple DOFs if not Wishart
    }


  }
  else if(summary=="log-normal" && prior=="log-normal")
  {
    #
  }

}
