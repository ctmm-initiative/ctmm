# solve a quadratic programming problem constrained to the probability simplex
# this method uses the KKT relations to identify what appears to be the current feasible/free region
# and then within that region we apply preconditioned conjugate gradient to the equality constrained optimization problem
PQP.solve <- function(G,FLOOR=NULL,p=NULL,lag=NULL,error=.Machine$double.eps,PC="circulant",IG=NULL,MARKOV=NULL)
{
  silent <- TRUE

  if(length(dim(G))==2) { fast <- FALSE } else { fast <- TRUE }
  # are we doing direct matrix calculations or FFT tricks?
  if(!fast)
  {
    n <- dim(G)[1]

    # O(n^2) matrix-vector multiplication
    G.VEC <- function(V,inverse=FALSE) { c(G %*% V) }

    # stuff needed to apply preconditioners
    if(PC=="IID")
    {
      ###############################
      #IID PRECONDITIONER: diagonal + constant matrix
      # pre-conditioner eigen-values
      PC0 <- (G[1,1] + (n-1)*G[1,n])
      PC1 <- (G[1,1] - G[1,n])
      PC.UPDATE <- function()
      {
        PC0 <<- (G[1,1] + (m-1)*G[1,n])
        PC1 <<- (G[1,1] - G[1,n])
      }
    }
    else if(PC=="banded")
    {
      ######################
      # Tri-banded + Constant matrix
      CONS <- G[1,n]
      DIAG <- G[1,1] - CONS
      BAND <- sapply(1:(n-1),function(i) { G[i,i+1] }) - CONS
      PC.UPDATE <- function()
      {
        INDEX <- (1:n)[FREE]
        BAND <<- sapply(1:(m-1), function(i) { G[INDEX[i],INDEX[i+1]] } ) - CONS
      }
    }
    else if(PC=="circulant")
    {
      ###############################
      # PSD CIRCULANT PRECONDITIONER
      # T. Chan's PSD preconditioner for symmetric, almost-Toeplitz matrices
      # working in free space only
      C <- G
      # what diagonal am I on?
      dg <- (row(C)-col(C)) %% n
      # split by diagonal
      C <- split(C,dg)
      # average diagonals -> circulant diagonals
      C <- sapply(C,mean)
      # diagonalize the circulant matrix
      C <- FFT(C)
      # inverse matrix
      C <- 1/Re(C)
      PC.UPDATE <- function()
      {
        # working in free space only
        C <- G[FREE,FREE]
        # what diagonal am I on
        dg <- (row(C)-col(C)) %% m
        # split by diagonal
        C <- split(C,dg)
        # average diagonals - circulant diagonals
        C <- sapply(C,mean)
        # diagonalize the circulant matrix
        C <- FFT(C)
        # inverse matrix
        C <<- 1/Re(C)
      }
    }
    else if(PC=="Toeplitz")
    {
      stop("!fast & Topelitz not yet implemented.")
    }
  }
  else ### using FFT tricks ###########################################################
  {
    n <- length(FLOOR) # number of sampled times
    N <- length(G) # numer of grid times
    q <- 1-p

    # stuff needed to apply preconditioners
    if(PC=="IID")
    {
      # # pre-conditioner eigen-values
      G1 <- G[1] ; GN <- G[N]
      # conditioning transformation eigen-values
      PC0 <- (G1 + (n-1)*GN)
      PC1 <- (G1 - GN)
      PC.UPDATE <- function()
      {
        # conditioning transformation eigen-values
        PC0 <<- (G1 + (m-1)*GN)
        PC1 <<- (G1 - GN)
      }
    }
    else if(PC=="banded")
    {
      stop("fast banded not yet implemented")
    }
    else if(PC=="circulant")
    {
      # unmodified copy of 'G' corresponding to lag
      g <- G
      # sparser lags approximating free space
      LAG <- seq(0,last(lag),length.out=n)
      # interpolate over approximate free space
      C <- stats::approx(lag,g,LAG,rule=2)$y
      R <- (0:(n-1))/n
      C <- (1-R)*C + R*c(0,rev(C[-1]))
      C <- FFT(C)
      C <- 1/Re(C)

      PC.UPDATE <- function()
      {
        LAG <- seq(0,last(lag),length.out=m)
        C <- stats::approx(lag,g,LAG,rule=2)$y
        R <- (0:(m-1))/m
        C <- (1-R)*C + R*c(0,rev(C[-1]))
        C <- FFT(C)
        C <<- 1/Re(C)
      }
    }

    # circulant embedding for FFT matrix multiplication
    T.C.embed <- function(M)
    {
      M <- c(M,0,rev(M[-1]))
      # Fourier transform
      M <- FFT(M)
      M <- Re(M)
      return(M)
    }
    G <- T.C.embed(G)

    if(PC=="Toeplitz")
    {
      # separate out the initial discontinuity
      IG0 <- IG[1] - IG[2]
      IG[1] <- IG[2]
      # now the preconditioner is IG0*Id + IG
      IG <- T.C.embed(IG)
    }

    # O(n log n) matrix.vector multiplication
    G.VEC <- function(V,inverse=FALSE)
    {
      # grid embedding
      Vgrid <- array(0,N)
      Vgrid[FLOOR+1] <- q*V
      Vgrid[FLOOR] <- Vgrid[FLOOR] + p*V

      # circulant embedding
      Vgrid <- c(Vgrid,rep(0,N))
      # Fourier transform / diagonalize
      Vgrid <- FFT(Vgrid)

      # matrix multiplication in diagonalized basis
      if(!inverse)
      { Vgrid <- G * Vgrid }
      else
      { Vgrid <- IG * Vgrid }

      # inverse Fourier transform
      Vgrid <- IFFT(Vgrid)
      # Toepltiz part
      Vgrid <- Vgrid[1:N]
      Vgrid <- Re(Vgrid)

      # map back from grid
      V <- q*Vgrid[FLOOR+1]
      V <- V + p*Vgrid[FLOOR]

      return(V)
    }

    if(PC=="Toeplitz")
    {
      # inverse approximation normalization
      IGM <- rep(1,n) # perhaps the most important vector to fix the norm with
      IGM <- (IG0 * IGM) + G.VEC(IGM,inverse=TRUE)
      IGM <- G.VEC(IGM,inverse=FALSE)
      IGM <- abs(IGM)
      IGM <- mean(IGM) # average and quadratic form now
      IGM <- 1/IGM
      PC.UPDATE <- function()
      {
        IGM <- FREE
        IGM <- (IG0 * IGM) + G.VEC(IGM,inverse=TRUE)
        IGM <- G.VEC(IGM,inverse=FALSE)
        IGM <- IGM[FREE]
        IGM <- abs(IGM)
        IGM <- mean(IGM) # average and quadratic form now
        IGM <- 1/IGM
      }
    }

    if(PC=="direct")
    {
      # calculate n*n G matrix
      G <- sapply(1:n,function(i){ V <- rep(0,n) ; V[i] <- 1 ; G.VEC(V) })
      # somehow this ends up missing a sparse number of times
      G <- pmax(G,t(G))

      G.VEC <- function(V,inverse=FALSE) { c(G %*% V) }
    }
  }

  ########################################################
  ######### DEFINE PRECONDITIONERS #######################
  # preconditioner application code
  if(PC=="IID")
  {
    # 1/sqrt(G_IID) pre-conditioning matrix multiplication
    # this is supposed to be a rough approximation of 1/G
    PC.VEC <- function(V)
    { (1/PC1)*V + (1/PC0-1/PC1)*mean(V[FREE])*FREE }
  }
  else if(PC=="banded")
  {
    # tri-band solver without constant
    BAND.SOLVE <- function(VEC)
    {
      # forward elimination
      B <- rep(BAND[1]/DIAG,m)
      V <- rep(VEC[1]/DIAG,m)
      for(i in 2:(m-1))
      {
        # can precalculate some of this
        A <- DIAG - BAND[i-1]*B[i-1]
        B[i] <- BAND[i]/A
        V[i] <- (VEC[i] - BAND[i-1]*V[i-1])/A
      }
      A <- DIAG - BAND[m-1]*B[m-1]
      V[m] <- (VEC[m] - BAND[m-1]*V[m-1])/A

      # backwards substitution
      for(i in (m-1):1)
      {
        V[i] <- V[i] - B[i]*V[i+1]
      }

      return(V)
    }

    # Woodbury matrix identity to account for constant atop tri-band
    PC.VEC <- function(V)
    {
      ONE <- rep(1,m)
      TBF <- BAND.SOLVE(ONE)
      TBV <- BAND.SOLVE(V[FREE])

      V[FREE] <- TBV - c(ONE %*% TBV)/(1/CONS + c(ONE %*% TBF))*TBF

      return(V)
    }
  }
  else if(PC=="circulant")
  {
    PC.VEC <- function(V)
    {
      W <- V[FREE]
      W <- FFT(W)
      W <- C * W
      W <- IFFT(W)
      W <- Re(W)
      V[FREE] <- W
      return(V)
    }
  }
  else if(PC=="Toeplitz")
  {
    PC.VEC <- function(V)
    {
      V[FREE] <- IGM*( IG0*V[FREE] + G.VEC(V,inverse=TRUE)[FREE] )
      return(V)
    }
  }
  else if(PC=="Markov")
  {
    # This is all exactly the same, fast and slow
    # following the notation of Rybicki (1994)
    PC.R <- exp(-diff(MARKOV$t))
    PC.E <- PC.R/(1-PC.R^2)
    PC.D <- PC.R*PC.E
    PC.D <- 1 + c(0,PC.D) + c(PC.D,0)
    PC.UPDATE <- function()
    {
      PC.R <<- exp(-diff(MARKOV$t[FREE]))
      PC.E <<- PC.R/(1-PC.R^2)
      PC.D <<- PC.R*PC.E
      PC.D <<- 1 + c(0,PC.D) + c(PC.D,0)
    }

    # tri-banded inverse
    MARKOV.SOLVE <- function(V)
    {
      W <- V[FREE]
      W <- PC.D*W - c(PC.E*W[-1],0) - c(0,PC.E*W[-m])
      V[FREE] <- W/MARKOV$VAR
      return(V)
    }

    # use Woodbury matrix identity to add constant term
    PC.VEC <- function(V)
    {
      MSF <- MARKOV.SOLVE(FREE)
      MSV <- MARKOV.SOLVE(V)

      V <- MSV - c(FREE %*% MSV)/(1/MARKOV$ERROR + c(FREE %*% MSF))*MSF

      return(V)
    }
  }

  # find the best linear combination of vectors for solution to V=solve(G,1)
  L.SEARCH <- function(V,G.V)
  {
    k <- length(V)
    k.V <- length(G.V)

    # populate remainder of G.V if necessary
    if(k.V < k)
    {
      for(i in (k.V+1):k)
      {
        G.V[[i]] <- G.VEC(V[[i]])
        G.V[[i]][ACTIVE] <- 0
      }
    }

    # matrix of inner productions WRT G
    M <- array(0,c(k,k))
    for(i in 1:k)
    {
      # fill upper triangle
      for(j in i:k)
      { M[i,j] <- c(V[[i]] %*% G.V[[j]]) }
      # copy lower triangle
      for(j in 1:i)
      { M[j,i] <- M[i,j] }
    }

    # vector of inner products WRT 1
    B <- sapply(V,sum)

    # safe solver in case of colinearity
    A <- PDsolve(M,pseudo=TRUE) %*% B

    V <- lapply(1:k,function(i) { A[i] * V[[i]] })
    G.V <- lapply(1:k,function(i) { A[i] * G.V[[i]] })

    V <- Reduce('+',V)
    G.V <- Reduce('+',G.V)

    return(list(V=V,G.V=G.V))
  }

  #####################
  # CHECK KKT CONDITIONS TO FIX ACTIVE/FREE DIMENSIONS
  ACTIVE <- rep(FALSE,n) # current active constraints
  FREE <- !ACTIVE # current feasible dimensions
  m <- n # number of free dimensions
  CHANGE <- TRUE # did our constraints change?, then keep working
  MISE <- rep(Inf,3) # last 3 MISE objectives - the 3>1 is arbitrary
  LAMBDA <- rep(1/n,3) # last 3 normalization constants
  KKT <- function()
  {
    # we're going to toss the negative probabilities out of the feasible region
    ZERO <- (P<=0)
    P[ZERO] <<- 0

    # update past esitmates (un-normalized)
    P.OLD <<- rbind(P,P.OLD[1:2,])

    # normalization constant
    LAMBDA[2:3] <<- LAMBDA[1:2] # update past lambdas
    LAMBDA[1] <<- 1/sum(P) # the current Lagrange multiplier

    G.P <<- G.VEC(P)

    # update MISE objective
    MISE[2:3] <<- MISE[1:2]
    MISE[1] <<- LAMBDA[1]^2 * c(P %*% G.P)

    # optimization on the probability simplex determines the active constraints from KKT equations
    GRAD <- 1 - G.P # -GRAD/LAMBDA, where GRAD == slack at solution (KKT)

    # currently active inequality constraints
    ACTIVE -> OLD
    ACTIVE <<- ZERO & (GRAD <= 0) # in these dimensions, the gradient is trying to push through the boundary, creating slack (KKT)
    CHANGE <<- any(ACTIVE-OLD) # did our active constraints change?
    FREE <<- !ACTIVE # free space to work in
    m <<- sum(FREE) # size of current free space

    return(NULL)
  }

  ############################
  # INITIALIZE SOLUTION
  # use previous solution
  if(get("EMPTY",pos=PQP.env))
  {
    # not yet a normalized probability
    P <- rep(1,n)

    # approximate solve(G,1) consistent with preconditioner
    if(PC!="direct")
    { P <- PC.VEC(P) }
    else # might as well do something
    { P[FREE] <- qr.solve(G[FREE,FREE],rep(1,m),tol=.Machine$double.eps) }
  }
  else # use preconditioner solution
  {
    P <- get("P",pos=PQP.env)
  }
  G.P <- P # initializing variable with junk
  P.OLD <- rbind(P,P,P) # more junk x3
  # check KKT conditions, incase we have a non-trivial initial guess
  KKT() # this returns a number now... does that matter in R?
  # as a side effect this also evaluates G.P, which we need

  CHANGE <- TRUE
  STEPS <- 0
  CHANGES <- 0
  while(CHANGE)
  {
    CHANGES <- CHANGES + 1

    # DIRECT SOLVER #############
    if(PC=="direct")
    {
      STEPS <- STEPS + 1
      P[FREE] <- qr.solve(G[FREE,FREE],rep(1,m),tol=.Machine$double.eps)
    }
    else # CONJUGATE GRADIENT SOLVER ########################
    {
      PC.UPDATE()
      G.P[ACTIVE] <- 0

      # we can do a better job of fixing G[FREE,FREE]*P[FREE] == 1[FREE] than just tossing out the zeros from P
      # P[FREE] <- A[1]*P[FREE] + A[2]*1[FREE] + A[3]*PC(1[FREE]) + A[4]*dP[FREE]
      P <- list(P)
      P[[2]] <- as.numeric(FREE)
      if(PC!="IID" & (!get("EMPTY",pos=PQP.env) || CHANGES>1)) # otherwise this is the same as [[2]] || [[1]]
      {
        P[[3]] <- PC.VEC(P[[2]])
        P[[3]][ACTIVE] <- 0
      }
      G.P <- list(G.P)

      P <- L.SEARCH(P,G.P)
      G.P <- P$G.V
      P <- P$V
      # The improvement with this seems rather marginal, but worth keeping...

      # STEEPEST DESCENT INITIALIZATION INTO CONJUGATE GRADIENT ON !CHANGE
      # optimization in the free space to actually solve for Q=solve(G_FREE,1_FREE) proportional to solution by P=LAMBDA*Q
      RES <- FREE - G.P # this is -GRAD in the conditioned free space, using conjugate-gradient notation of residual
      Z <- PC.VEC(RES) # preconditioned conjugate gradient variable
      CONJ <- Z # this is often called 'p', but that's confusing because we are working with probabilities
      # optimize step length in free space assumption

      # CONJUGATE GRADIENT (STEPEST DESCENT ON FIRST STEP)
      ERROR <- Inf
      CG.STEPS <- 0
      while(ERROR > error && CG.STEPS < n)
      {
        # optimize step length in free space assumption
        G.CONJ <- G.VEC(CONJ)
        G.CONJ[ACTIVE] <- 0
        C.G.C <- c(CONJ %*% G.CONJ)
        if(abs(C.G.C)<=.Machine$double.eps) { if(!silent) { message("PCCG divergence in alpha") } ; break }
        Z.RES <- c(Z %*% RES)
        A <- Z.RES / C.G.C
        dP <- A*CONJ # need this for stopping condition

        # make sure Lagrangian/objective is decreasing
        dL <- c((P+dP/2) %*% G.CONJ)*A - sum(dP)
        # the solution going forward is actually worse
        if(dL>=0) { if(!silent) { message("PCCG stopped converging") } ; break }

        P <- P + dP
        dRES <- -A*G.CONJ # need this for Polak-Ribiere formula
        RES <- RES + dRES
        Z <- PC.VEC(RES)
        Z.dRES <- c(Z %*% dRES) # Polak-Ribiere formula
        Z.dRES <- max(0,Z.dRES) # reset to steepest descent if necessary
        if(abs(Z.RES)<=.Machine$double.eps) { if(!silent) { message("PCCG divergence in beta") } ; break }
        B <- Z.dRES / Z.RES
        CONJ <- Z + B*CONJ

        ####################################
        # make this the decrease in Lagrangian or something ?
        # weights are all on the same scale.................
        ERROR <- max(abs(dP)) * n

        CG.STEPS <- CG.STEPS + 1
      }
      STEPS <- STEPS + CG.STEPS
    }

    # CHECK KKT CONDITIONS AND RESET ACTIVE/FREE SPACE
    KKT()

    # make sure the MISE is decreasing and we are not stuck in a rare numerical-indiscrimination / solution-boundary loop
    if(all(MISE[1] >= MISE[2:3])) { break }
  }

  # select best solution of past few (in case of rare boundary loop)
  MIN <- which.min(MISE)
  MISE <- MISE[MIN]
  P <- P.OLD[MIN,]

  # STORE SOLUTION FOR LATER OPTIM EVALUATIONS
  assign("P",P,pos=PQP.env)
  assign("EMPTY",FALSE,pos=PQP.env) # SET TO FALSE TO ACTIVATE MEMORY !!!!!!!!!!!!!!!!!!!!!!!!!!

  P <- LAMBDA[MIN] * P

  RETURN <- list(P=P,MISE=MISE,STEPS=STEPS,CHANGES=CHANGES)
  return(RETURN)
}

# store solutions to pass between optimize evaluations
PQP.env <- new.env()

# empty out this environment when we are done
empty.env <- function(ENV)
{
  STUFF <- ls(ENV)
  rm(list=STUFF,pos=ENV)
  assign("EMPTY",TRUE,pos=ENV)
}
