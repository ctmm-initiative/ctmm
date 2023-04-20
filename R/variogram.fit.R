####################################
# guess variogram model parameters #
####################################
variogram.guess <- function(variogram,CTMM=ctmm())
{
  # make sure formatting is correct
  CTMM$axes <- attr(variogram,"info")$axes
  CTMM <- ctmm.prepare(variogram,CTMM,precompute=FALSE,tau=FALSE)

  # guess at some missing parameters
  n <- length(variogram$lag)

  if(n>1) # more than lag-0
  {
    # variance estimate
    sigma <- mean(variogram$SVF[2:n])

    # peak curvature estimate
    # should occur at short lags
    v2 <- 2*max((variogram$SVF/variogram$lag^2)[2:n])

    # free frequency
    Omega2 <- v2/sigma
    Omega <- sqrt(Omega2)

    # peak diffusion rate estimate
    # should occur at intermediate lags
    # diffusion parameters
    D <- (variogram$SVF/variogram$lag)[2:n]
    # index of max diffusion
    tauD <- which.max(D)
    # max diffusion
    D <- D[tauD]
    # time lag of max diffusion
    tauD <- variogram$lag[tauD]

    # average f-rate
    f <- -log(D/(sigma*Omega))/tauD

    tau <- c(sigma/D,D/v2)
    names(tau) <- c("position","velocity")

    if(length(CTMM$tau)==0)
    { CTMM$tau <- tau }
    else if(length(CTMM$tau)==1)
    { CTMM$tau[2] <- min(CTMM$tau[1],tau[2]) }

    if(!CTMM$range) { sigma <- D ; CTMM$tau[1] <- Inf }

    # preserve orientation and eccentricity if available/necessary
    if(is.null(CTMM$sigma))
    { CTMM$sigma <- sigma }
    # else
    # {
    #   CTMM$sigma <- CTMM$sigma@par
    #   CTMM$sigma[1] <- sigma / cosh(CTMM$sigma[2]/2)
    # }
  }

  # don't overwrite or lose ctmm parameters not considered here
  model <- as.list(CTMM) # getDataPart deletes names()
  model$info <- attr(variogram,"info")
  model <- do.call("ctmm",model)

  return(model)
}


######################################################################
# visual fit of the variogram
# manipulate abstraction layer
######################################################################
variogram.fit <- function(variogram,CTMM=ctmm(),name="GUESS",fraction=0.5,interactive=TRUE,...)
{
  CLASS <- class(variogram)[1]
  if(CLASS %nin% c('variogram','data.frame')) { stop("variogram is class ",CLASS) }

  # R CHECK CRAN BUG BYPASS
  z <- NULL
  sigma <- NULL
  tau1 <- NULL
  tau2 <- NULL
  omega.period <- NULL
  circle.period <- NULL
  error <- NULL
  rm(z,sigma,tau1,tau2,omega.period,circle.period,error)

  if(interactive && !manipulate::isAvailable()) { interactive <- FALSE }
  envir <- .GlobalEnv

  # fill in missing parameters non-destructively
  CTMM <- variogram.guess(variogram,CTMM)
  if(!interactive) { return(CTMM) }

  STUFF <- variogram.fit.backend(variogram,CTMM=CTMM,fraction=fraction)
  DF <- STUFF$DF

  manlist <- lapply(1:nrow(DF),function(r){ manipulate::slider(min=DF$min[r],max=DF$max[r],initial=DF$initial[r],label=DF$label[r],step=DF$step[r]) })
  names(manlist) <- rownames(DF)

  # put error checkbox option here
  CHECK <- DF$min==0 & DF$max==1 & DF$step==1
  if(any(CHECK))
  {
    CHECK <- which(CHECK)
    for(r in CHECK) { manlist[[r]] <- manipulate::checkbox(initial=as.logical(DF$initial[r]),label=DF$label[r]) }
  }

  # R CHECK CRAN BUG BYPASS
  store <- NULL ; rm(store)
  # storage button
  manlist <- c(manlist, list(store = manipulate::button(paste("Save to",name))))

  manipulate::manipulate(
    {
      slider_names <- rownames(DF)
      # zoom slider is not needed here
      slider_names <- slider_names[!(slider_names %in% "z")]
      # must use as.name to avoid immediate evaluation
      slider_values <- lapply(slider_names, as.name)
      names(slider_values) <- slider_names
      CTMM <- do.call(STUFF$storer, slider_values)

      fraction <- STUFF$fraction(z)

      if(store) { assign(name,CTMM,envir=envir) }

      plot.variogram(variogram,CTMM=CTMM,fraction=fraction,...)
    },
    manlist
  )

}


# data.frame DF: columns of label, initial, min, max, step
# row names of data.frame = variable names
# storer() function sets ctmm object
# fraction() sets plot fraction
variogram.fit.backend <- function(variogram,CTMM=ctmm(),fraction=0.5,b=4)
{
  # R CHECK CRAN BUG BYPASS
  z <- NULL
  sigma <- NULL
  tau1 <- NULL
  tau2 <- NULL
  omega.period <- NULL
  circle.period <- NULL
  error <- NULL
  rm(z,sigma,tau1,tau2,omega.period,circle.period,error)

  RES <- 1000

  # base error parameters
  TYPE <- DOP.match(CTMM$axes)
  UERE <- attr(variogram,"UERE")$UERE[,TYPE]
  names(UERE) <- rownames(attr(variogram,"UERE")$UERE)

  m <- 2 # slider length relative to point guestimate
  n <- length(variogram$lag)

  # parameters for logarithmic slider
  # b <- 4
  min.step <- 10*variogram$lag[2]/variogram$lag[n]
  #min.step <- max(min.step,fraction)

  # manipulation controls
  z=1+log(fraction,b)
  MIN <- 1+log(min.step,b)
  MIN <- min(z,MIN)
  DF <- data.frame(min=MIN,max=1,initial=z,label="zoom",step=log(min.step,b)/RES,stringsAsFactors=FALSE)
  NAMES <- c("z")

  K <- length(CTMM$tau)
  if(K==0) { CTMM$tau <- c(0,0) }
  else if(K==1) { CTMM$tau[2] <- 0 }

  range <- CTMM$range
  sigma <- mean(diag(CTMM$sigma))
  if(range)
  {
    sigma.unit <- unit(sigma,"area",concise=TRUE)
    sigma <- sigma / sigma.unit$scale
    label <- paste("\u03C3 variance (",sigma.unit$name,")",sep="")
  }
  else
  {
    sigma.unit <- unit(sigma,"diffusion",concise=TRUE)
    sigma <- sigma / sigma.unit$scale
    label <- paste("\u03C3 diffusion (",sigma.unit$name,")",sep="")
  }
  DF <- rbind(DF,data.frame(min=0,max=m*sigma,initial=sigma,label=label,step=sigma/RES,stringsAsFactors=FALSE))
  NAMES <- c(NAMES,"sigma")

  tau <- CTMM$tau
  tau1.unit <- unit(tau[1],"time",2,concise=TRUE)
  tau2.unit <- unit(tau[2],"time",2,concise=TRUE)
  tau[1] <- tau[1] / tau1.unit$scale
  tau[2] <- tau[2] / tau2.unit$scale

  omega <- CTMM$omega

  if(!omega) # position and velocity decay
  {
    label <- paste("\u03C4[position] (",tau1.unit$name,")",sep="")
    DF <- rbind(DF,data.frame(min=0,max=m*tau[1],initial=tau[1],label=label,step=tau[1]/RES,stringsAsFactors=FALSE))
    NAMES <- c(NAMES,"tau1")

    label <- paste("\u03C4[velocity] (",tau2.unit$name,")",sep="")
    DF <- rbind(DF,data.frame(min=0,max=m*tau[2],initial=tau[2],label=label,step=tau[2]/RES,stringsAsFactors=FALSE))
    NAMES <- c(NAMES,"tau2")
  }
  else # decay and oscillation period
  {
    label <- paste("\u03C4[decay] (",tau1.unit$name,")",sep="")
    DF <- rbind(DF,data.frame(min=0,max=m*tau[1],initial=tau[1],label=label,step=tau[1]/RES,stringsAsFactors=FALSE))
    NAMES <- c(NAMES,"tau1")

    omega.period <- 2*pi/omega
    omega.unit <- unit(omega.period,"time",concise=TRUE)
    omega.period <- omega.period / omega.unit$scale
    label <- paste("\u03C4[period] (",omega.unit$name,")",sep="")
    DF <- rbind(DF,data.frame(min=0,max=m*omega.period,initial=omega.period,label=label,step=omega.period/RES,stringsAsFactors=FALSE))
    NAMES <- c(NAMES,"omega.period")
  }

  # circulation
  circle <- CTMM$circle
  if(circle)
  {
    circle.period <- 2*pi/circle
    circle.unit <- unit(circle.period,"time",concise=TRUE)
    circle.period <- circle.period / circle.unit$scale
    label <- paste("circulation period (",circle.unit$name,")",sep="")
    c1 <- min(0,m*circle.period)
    c2 <- max(0,m*circle.period)
    DF <- rbind(DF,data.frame(min=c1,max=c2,initial=circle.period,label=label,step=abs(circle.period)/RES,stringsAsFactors=FALSE))
    NAMES <- c(NAMES,"circle.period")
  }

  ## Depreciated error slider
  # e2 <- max(100,2*error[i])
  # label <- paste("error",names(UERE.DOF)[i],"(m)")
  # DF <- rbind(DF,data.frame(min=0,max=e2,initial=as.numeric(error[i]),label=label,step=e2/RES/2,stringsAsFactors=FALSE))

  # ERROR
  error <- any(CTMM$error>0)
  label <- paste("error (logical)")
  DF <- rbind(DF,data.frame(min=0,max=1,initial=sign(error),label=label,step=1,stringsAsFactors=FALSE))
  NAMES <- c(NAMES,"error")

  rownames(DF) <- NAMES

  # discard rows we don't want
  if(!range)
  {
    NAMES <- rownames(DF)
    DELETE <- which(NAMES=="tau1")
    DF <- DF[-DELETE,]
    tau1 <- Inf
  }

  if(K<1 && FALSE) # need to update
  {
    NAMES <- rownames(DF)
    DELETE <- which(NAMES=="tau1")
    DF <- DF[-DELETE,]
    tau1 <- 0
  }

  if(K<2 && FALSE) # need to update
  {
    NAMES <- rownames(DF)
    DELETE <- which(NAMES=="tau2")
    DF <- DF[-DELETE,]
    tau2 <- 0
  }

  # non-destructive parameter overwrite
  storer <- function(sigma=sigma,tau1=tau[1],tau2=tau[2],omega.period=omega.period,circle.period=circle.period,error=error)
  {
    # store trace, but preserve angle & eccentricity
    if(length(CTMM$axes)==2)
    { CTMM$sigma <- scale.covm(CTMM$sigma,(sigma*sigma.unit$scale)/mean(diag(CTMM$sigma))) }
    else
    { CTMM$sigma <- sigma*sigma.unit$scale }

    if(!CTMM$omega)
    { CTMM$tau <- c(tau1*tau1.unit$scale, tau2*tau2.unit$scale) }
    else
    {
      CTMM$tau <- c(1,1)*tau1*tau1.unit$scale
      CTMM$omega <- 2*pi/(omega.period * omega.unit$scale)
    }

    if(circle) { CTMM$circle <- 2*pi/(circle.period * circle.unit$scale) }

    CTMM$error <- error * UERE

    CTMM <- as.list(CTMM)
    CTMM$info <- attr(variogram,"info")
    CTMM <- do.call("ctmm",CTMM)
    if(any(CTMM$tau>0)) { CTMM$tau <- CTMM$tau[CTMM$tau>0] }
    else { CTMM$tau <- NULL }

    return(CTMM)
  }

  fraction <- function(z=z) { b^(z-1) }

  return(list(DF=DF,storer=storer,fraction=fraction))
}
