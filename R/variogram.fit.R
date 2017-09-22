####################################
# guess variogram model parameters #
####################################
variogram.guess <- function(variogram,CTMM=ctmm())
{
  # guess at some missing parameters
  n <- length(variogram$lag)

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

  CPF <- CTMM$CPF
  if(CPF) # frequency, rate esitmate
  {
    omega2 <- Omega2 - f^2
    if(f>0 && omega2>0)
    { tau <- c(2*pi/sqrt(omega2),1/f) }
    else # bad backup estimate
    {
      tau <- sqrt(2)/Omega
      tau <- c(2*pi*tau , tau)
    }
  }
  else # position, velocity timescale estimate
  { tau <- c(sigma/D,D/v2)}

  if(!CTMM$range) { sigma <- D ; tau[1] <- Inf }

  if(length(CTMM$tau)==0) { CTMM$tau <- tau }
  #else if(length(CTMM$tau)==1) { CTMM$tau[2] <- tau[2] }

  # preserve orientation and eccentricity if available/necessary
  if(is.null(CTMM$sigma))
  { CTMM$sigma <- sigma }
  # else
  # {
  #   CTMM$sigma <- CTMM$sigma@par
  #   CTMM$sigma[1] <- sigma / cosh(CTMM$sigma[2]/2)
  # }

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
  # R CHECK CRAN BUG BYPASS
  z <- NULL
  sigma <- NULL
  tau1 <- NULL
  tau2 <- NULL
  circle.period <- NULL
  error <- NULL
  rm(z,sigma,tau1,tau2,circle.period,error)

  if(interactive && !manipulate::isAvailable()) { interactive <- FALSE }
  envir <- .GlobalEnv

  # fill in missing parameters non-destructively
  CTMM <- variogram.guess(variogram,CTMM)
  if(!interactive) { return(CTMM) }

  STUFF <- variogram.fit.backend(variogram,CTMM=CTMM,fraction=fraction)
  DF <- STUFF$DF

  manlist <- lapply(1:nrow(DF),function(r){ manipulate::slider(min=DF$min[r],max=DF$max[r],initial=DF$initial[r],label=DF$label[r],step=DF$step[r]) })
  names(manlist) <- rownames(DF)

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
  circle.period <- NULL
  rm(z,sigma,tau1,tau2,circle.period)

  RES <- 1000

  error <- CTMM$error

  m <- 2 # slider length relative to point guestimate
  n <- length(variogram$lag)

  # parameters for logarithmic slider
  # b <- 4
  min.step <- 10*variogram$lag[2]/variogram$lag[n]
  #min.step <- max(min.step,fraction)

  # manipulation controls
  z=1+log(fraction,b)
  DF <- data.frame(min=1+log(min.step,b),max=1,initial=z,label="zoom",step=log(min.step,b)/RES,stringsAsFactors=FALSE)
  NAMES <- c("z")

  K <- length(CTMM$tau)
  if(K==1) { CTMM$tau[2] <- 0 }

  range <- CTMM$range
  sigma <- mean(diag(CTMM$sigma))
  if(range)
  {
    sigma.unit <- unit(sigma,"area",concise=TRUE)
    sigma <- sigma / sigma.unit$scale
    label <- paste("sigma variance (",sigma.unit$name,")",sep="")
  }
  else
  {
    sigma.unit <- unit(sigma,"diffusion",concise=TRUE)
    sigma <- sigma / sigma.unit$scale
    label <- paste("sigma diffusion (",sigma.unit$name,")",sep="")
  }
  DF <- rbind(DF,data.frame(min=0,max=m*sigma,initial=sigma,label=label,step=sigma/RES,stringsAsFactors=FALSE))
  NAMES <- c(NAMES,"sigma")

  tau <- CTMM$tau
  tau1.unit <- unit(tau[1],"time",2,concise=TRUE)
  tau2.unit <- unit(tau[2],"time",2,concise=TRUE)
  tau[1] <- tau[1] / tau1.unit$scale
  tau[2] <- tau[2] / tau2.unit$scale

  label <- paste("tau position (",tau1.unit$name,")",sep="")
  DF <- rbind(DF,data.frame(min=0,max=m*tau[1],initial=tau[1],label=label,step=tau[1]/RES,stringsAsFactors=FALSE))
  NAMES <- c(NAMES,"tau1")

  label <- paste("tau velocity (",tau2.unit$name,")",sep="")
  DF <- rbind(DF,data.frame(min=0,max=m*tau[2],initial=tau[2],label=label,step=tau[2]/RES,stringsAsFactors=FALSE))
  NAMES <- c(NAMES,"tau2")

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

  # error
  e2 <- max(100,2*error)
  DF <- rbind(DF,data.frame(min=0,max=e2,initial=as.numeric(error),label="error (m)",step=e2/RES/2,stringsAsFactors=FALSE))
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

  if(K==1)
  {
    NAMES <- rownames(DF)
    DELETE <- which(NAMES=="tau2")
    DF <- DF[-DELETE,]
    tau2 <- 0
  }

  # non-destructive parameter overwrite
  storer <- function(sigma=sigma,tau1=tau[1],tau2=tau[2],circle.period=circle.period,error=error)
  {
    # store trace, but preserve angle & eccentricity
    if(length(CTMM$axes)==2)
    {
      CTMM$sigma <- CTMM$sigma@par
      CTMM$sigma[1] <- sigma * sigma.unit$scale / cosh(CTMM$sigma[2]/2)
    }
    else
    { CTMM$sigma <- sigma }

    CTMM$tau <- c(tau1*tau1.unit$scale, tau2*tau2.unit$scale)
    if(circle) { CTMM$circle <- 2*pi/(circle.period * circle.unit$scale) }
    CTMM$error <- error

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