# Generic parallel lapply based on ctmmweb/R/5_parallel.R
#########################################################
# If cores=1, vanilla lapply
# If UNIX & cores>1, mclapply
# If windows & fast, vanilla lapply
# If windows & !fast, parLapplyLB


# detect relevant core number given fast boolean
detectCores <- function(...,fast=TRUE)
{
  if(fast && .Platform$OS.type=="windows") { return(1) } # Windows cannot fork
  else { return(parallel::detectCores(logical=FALSE,...)) }
}


# resolve number of cores to use given user input
# NULL uses all cores
# non-positive values reserves that many cores
resolveCores <- function(cores=1,fast=TRUE)
{
  if(is.null(cores) || is.na(cores)) { cores <- detectCores(fast=fast) }
  else if(cores<1) { cores <-  max(1,detectCores(fast=fast) + cores) }
  # Windows can't fork
  if(fast && .Platform$OS.type=="windows") { cores <- 1 }

  return(cores)
}


# smart parallel lapply
plapply <- function(X,FUN,...,cores=1,fast=TRUE)
{
  WINDOWS <- (.Platform$OS.type=="windows")
  cores <- resolveCores(cores,fast=fast)
  cores <- min(cores,length(X)) # cap cores

  if(cores==1 || (fast && WINDOWS)) { return(lapply(X,FUN,...)) }
  else if(!WINDOWS) { return(parallel::mclapply(X,FUN,...,mc.cores=cores)) }

  ### Windows parallel code below ###
  win_init = expression({requireNamespace("ctmm",quietly=TRUE)})

  cl <- parallel::makeCluster(cores,outfile="")
  # have to export parameter too because it's not available in remote
  parallel::clusterExport(cl,c("win_init"),envir=environment())
  parallel::clusterEvalQ(cl,eval(win_init))
  RESULT <- parallel::parLapplyLB(cl,X,FUN)
  parallel::stopCluster(cl)

  return(RESULT)
}
