getMethod <- function(fn,signature,...)
{
  meth <- NULL
  # S3 and S4 method dispatching is incompatible for no reason
  meth <- methods::getMethod(fn,signature=signature,optional=TRUE,...) # try S4 first
  if(is.null(meth)) { meth <- utils::getS3method(fn,signature[1],optional=TRUE,...) } # then try S3
  # due to new CRAN policy, we can no longer have internal S3 methods
  # work around here with '_' dispatch method names
  if(is.null(meth)) { try( meth <- get(paste0(fn,"_",signature)) , silent=TRUE) }
  if(is.null(meth)) { stop('Cannot find method ',fn,' for class ',signature) }
  return(meth)
}


# base match.arg cannot match NA ???
match.arg <- function(arg,choices,...)
{
  if(is.na(arg) && any(is.na(choices))) { return(NA) }
  else { return(base::match.arg(arg,choices,...)) }
}

# sort arguments by class
# sort.arg <- function(arg,sig)
# {
#
# }


# does this thing exist and, if so, is it true
is.good <- function(x) { !is.null(x) & !is.na(x) & x }
is.bad <- function(x) { is.null(x) | is.na(x) | !x }

# not in #
"%nin%" <- function(x, table) { match(x, table, nomatch = 0) == 0 }


# positive only sequence for for() loops so that for(i in 1:0) does nothing
"%:%" <- function(x,y)
{
  if(x<=y) { x <- x:y } else { x <- NULL }
  return(x)
}


composite <- function(n) { 2^ceiling(log(n,2)) }


# sinc functions
sinc <- Vectorize( function(x,SIN=sin(x))
{
  if(x==0)
  { return(1) }
  else
  { return(SIN/x) }
} )

sinch <- Vectorize( function(x,SINH=sinh(x))
{
  if(x==0)
  { return(1) }
  else
  { return(SINH/x) }
} )


# multivariate polygamma function
mpsigamma <- function(x,deriv=0,dim=1)
{
  PSI <- 1 - 1:dim
  PSI <- x + PSI/2
  if(deriv>=0) { PSI <- sapply(PSI,function(p) psigamma(p,deriv=deriv)) }
  else if(deriv==-1) { PSI <- sapply(PSI,function(p) lgamma(p)) }
  else { stop("Derivative ",deriv+1," of log(Gamma(x)) not supported.") }
  PSI <- sum(PSI)
  return(PSI)
}

### S4 ###

# forwarding function for list of a particular datatype
zoom.list <- function(x,...)
{
  CLASS <- class(x[[1]])[1]
  #utils::getS3method("zoom",CLASS)(x,...)
  getMethod("zoom",signature=CLASS)(x,...)
}
methods::setMethod("zoom",signature(x="list"), function(x,...) zoom.list(x,...))


### S3 ###

# forwarding function for list of a particular datatype
log.list <- function(x,...)
{
  CLASS <- class(x[[1]])[1]
  getMethod("log",CLASS)(x,...)
}
# this doesn't work outside of ctmm


# forwarding function for list of a particular datatype
mean.list <- function(x,...)
{
  CLASS <- class(x[[1]])[1]
  getMethod("mean",CLASS)(x,...)
}
#methods::setMethod("mean",signature(x="list"), function(x,...) mean.list(x,...))

# forwarding function for list of a particular datatype
median.list <- function(x,na.rm=FALSE,...)
{
  CLASS <- class(x[[1]])[1]
  getMethod("median",CLASS)(x,na.rm=na.rm,...)
}

# forwarding function for list of a particular datatype
plot.list <- function(x,...)
{
  CLASS <- class(x[[1]])[1]
  getMethod("plot",CLASS)(x,...)
}
#methods::setMethod("plot",signature(x="list"), function(x,...) plot.list(x,...))

# forwarding function for list of a particular datatype
summary.list <- function(object,...)
{
  # recurse if necessary
  CLASS <- "list"
  DATA <- object
  while(CLASS=="list")
  {
    DATA <- DATA[[1]]
    CLASS <- class(DATA)
  }

  getMethod("summary",CLASS)(object,...)
}

# forwarding function for list of a particular datatype
writeVector.list <- function(x,filename,...)
{
  CLASS <- class(x[[1]])[1]
  if(missing(filename))
  { getMethod("writeVector",methods::signature(x=CLASS,filename="missing"))(x,filename=filename,...) }
  else
  { getMethod("writeVector",methods::signature(x=CLASS,filename="character"))(x,filename=filename,...) }
}
methods::setMethod("writeVector",methods::signature(x="list",filename="character"), function(x,filename,...) writeVector.list(x,filename,...) )
methods::setMethod("writeVector",methods::signature(x="list",filename="missing"), function(x,filename,...) writeVector.list(x,filename,...) )


# replace NA elements
na.replace <- function(x,rep)
{
  REP <- is.na(x)
  x[REP] <- rep[REP]
  return(x)
}

########################
# 0/0 -> NaN -> to
# fixes a priori known limits
nant <- function(x,to)
{
  NAN <- is.na(x) # TRUE for NaN and NA
  if(any(NAN))
  {
    to <- array(to,length(x))
    x[NAN] <- to[NAN]
  }
  return(x)
}

# fix for infite PD matrix
# useful after nant(x,Inf)
inft <- function(x,to=0)
{
  INF <- diag(x)==Inf
  if(any(INF))
  {
    # force positive definite
    x[INF,] <- x[,INF] <- 0
    # restore infinite variance
    diag(x)[INF] <- Inf
  }
  return(x)
}


# parity tests
is.even <- Vectorize(function(x) {x %% 2 == 0})

is.odd <- Vectorize(function(x) {x %% 2 != 0})


# last element of array
last <- function(vec) { vec[length(vec)] }
first <- function(vec) { vec[1] }
# assign to last element... doesn't work
# "last<-" <- function(vec,ass)
# {
#   vec[length(vec)] <- ass
#   return(vec)
# }


# prepend to a vector
prepend <- function(x,values,before=1)
{ append(x,values,after=before-1) }


# CLAMP A NUMBER
clamp <- function(num,min=0,max=1)
{ ifelse(num<min,min,ifelse(num>max,max,num)) }


# PAD VECTOR
pad <- function(vec,size=length(vec),diff=size-length(vec),padding=0,side=+1)
{
  # this is now the pad length instead of total length
  padding <- array(padding,diff)

  if(side>0)
  { vec <- c(vec,padding) }
  else if(side<0)
  { vec <- c(padding,vec) }

  return(vec)
}

# row pad for data frames / matrices
rpad <- function(mat,size=nrow(mat),diff=size-nrow(mat),padding=0,side=+1)
{
  mat <- cbind(mat)
  COL <- ncol(mat)
  padding <- array(padding,c(diff,COL))
  colnames(padding) <- colnames(mat)

  if(side>0)
  { mat <- rbind(mat,padding) }
  else if(side<0)
  { mat <- rbind(padding,mat) }

  return(mat)
}

# pad both sides of a matrix
mpad <- function(mat,size=max(dim(mat)),diff=size-max(dim(mat)),padding=0,side=+1,padname=NULL,...)
{
  DIM <- dim(mat)

  mat <- rpad(mat,size,diff,padding=padding,side=side)
  mat <- t(mat)
  mat <- rpad(mat,size,diff,padding=padding,side=side)
  mat <- t(mat)

  if(!is.null(padname))
  {
    if(side<0)
    { SUB <- 1:diff }

    if(side>0)
    { SUB <- DIM[1] + 1:diff }

    rownames(mat)[SUB] <- padname

    if(side>0)
    { SUB <- DIM[2] + 1:diff }

    colnames(mat)[SUB] <- padname
  }

  return(mat)
}


#remove rows and columns by name
rm.name <- function(object,name)
{
  if(length(dim(object))==2)
  { object <- object[! rownames(object) %in% name,! colnames(object) %in% name,drop=FALSE] }
  else
  { object <- object[! names(object) %in% name] }
  return(object)
}


# put in a list if not in a list already
listify <- function(x)
{
  if(is.null(x)) { return(x) }

  if(class(x)[1] != "list")
  {
    x <- list(x)
    names(x) <- attr(x[[1]],'info')$identity
  }
  return(x)
}


# rename elements of an object
rename <- function(object,name1,name2)
{
  IND <- which(names(object)==name1)
  names(object)[IND] <- name2
  return(object)
}

rename.matrix <- function(object,name1,name2)
{
  NAMES <- dimnames(object)[[1]]
  IND <- which(NAMES==name1)
  NAMES[IND] <- name2
  dimnames(object) <- list(NAMES,NAMES)
  return(object)
}


# glue strings together if they different
glue <- function(...)
{
  x <- c(...) # removes NULLs
  x <- unique(x) # removes matchs
  x <- paste(x) # paste different strings together
  return(x)
}


# from
capitalize <- function(s)
{
  substr(s,1,1) <- toupper(substr(s,1,1))
  s
}


simplify.formula <- function(x)
{
  #x <- as.character(x)[2]
  x <- stats::terms(x,simplify=TRUE)
  x <- as.character(x)[2]
  x <- paste("~",x)
  x <- eval(parse(text=x))
  return(x)
}


copy <- function(from,to)
{
  NAMES <- names(from)
  for(n in NAMES){ to[[n]] <- from[[n]] }
  return(to)
}

mid <- function(x)
{
  n <- length(x)
  (x[-1]+x[-n])/2
}


name.list <- function(x)
{
  if(class(x)[1]=="list" && !length(names(x)))
  {
    NAMES <- sapply(x,function(y){attr(y,"info")$identity})
    if(class(NAMES)[1]=="character") { names(x) <- NAMES }
  }
  return(x)
}
