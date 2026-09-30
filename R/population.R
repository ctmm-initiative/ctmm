# effective individual of a population (RSF sans missing variance)
pop2ind <- function(data,UD,smooth=TRUE,...)
{
  CTMM <- UD@CTMM
  isotropic <- CTMM$isotropic
  axes <- CTMM$axes

  data <- listify(data)
  n <- length(data)
  if(n==1)
  { FIT <- list(CTMM) }
  else
  { FIT <- CTMM$CTMM }

  for(i in 1:n)
  {
    # smooth the data, but don't drop
    if(smooth && any(FIT[[i]]$error>0) && is.bad(attr(data[[i]],"info")$smoothed))
    {
      data[[i]][,c(axes,GEO)] <- predict(data[[i]],CTMM=FIT[[i]],t=data[[i]]$t,complete=TRUE)[,c(axes,GEO)]
      data[[i]]@info$smoothed <- TRUE
    }

    if(!FIT[[i]]$isotropic)
    {
      message("Use isotropic=TRUE before rsf.fit")

      if("ISO" %in% names(FIT[[i]]))
      { ISO <- FIT[[i]]$ISO }
      else
      {
        ISO <- simplify.ctmm(FIT[[i]],'minor')
        if(trace) { message("Fitting isotropic autocorrelation model.") }
        ISO <- ctmm.fit(data[[i]],ISO,trace=max(trace-1,0))
      }
      FIT[[i]] <- ISO

      if(n==1)
      {
        CTMM <- ISO
        UD@CTMM <- ISO
        UD$DOF.area <- DOF.area(ISO)
      }
      else
      {
        CTMM$CTMM[[i]] <- ISO
        if(i==n)
        {
          CTMM <- mean(CTMM$CTMM)
          UD@CTMM <- CTMM
          UD$DOF.area <- DOF.area(CTMM)
        }
      } # end pop ISO fix
    } # end ISO fix
  } # end data smoothing and iso fix

  # extract weights and structure data
  if(n==1)
  { data <- data[[1]] }
  else
  {
    UD$weights <- unlist(UD$w.list)

    # data <- lapply(data,function(d){d[,c('t',axes)]})
    info <- mean_info(data)
    data <- do.call(rbind,data)
    data <- new.telemetry(data,info=info)
  }

  R <- list(data=data,UD=UD)
  return(R)
}
