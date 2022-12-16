# base match.arg cannot match NA ???
match.arg <- function(arg,choices,...)
{
  if(is.na(arg) && any(is.na(choices))) { return(NA) }
  else { return(base::match.arg(arg,choices,...)) }
}


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
  methods::getMethod("zoom",signature=CLASS)(x,...)
}
methods::setMethod("zoom",signature(x="list"), function(x,...) zoom.list(x,...))


### S3 ###

# forwarding function for list of a particular datatype
log.list <- function(x,...)
{
  CLASS <- class(x[[1]])[1]
  utils::getS3method("log",CLASS)(x,...)
}
# this doesn't work outside of ctmm


# forwarding function for list of a particular datatype
mean.list <- function(x,...)
{
  CLASS <- class(x[[1]])[1]
  utils::getS3method("mean",CLASS)(x,...)
}
#methods::setMethod("mean",signature(x="list"), function(x,...) mean.list(x,...))

# forwarding function for list of a particular datatype
median.list <- function(x,na.rm=FALSE,...)
{
  CLASS <- class(x[[1]])[1]
  utils::getS3method("median",CLASS)(x,na.rm=na.rm,...)
}

# forwarding function for list of a particular datatype
plot.list <- function(x,...)
{
  CLASS <- class(x[[1]])[1]
  utils::getS3method("plot",CLASS)(x,...)
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

  utils::getS3method("summary",CLASS)(object,...)
}

# forwarding function for list of a particular datatype
writeShapefile.list <- function(object,folder,file=NULL,...)
{
  CLASS <- class(object[[1]])[1]
  utils::getS3method("writeShapefile",CLASS)(object,folder,file=file,...)
}


# replace NA elements
na.replace <- function(x,rep)
{
  REP <- is.na(x)
  x[REP] <- rep[REP]
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
pad <- function(vec,size,padding=0,side="right")
{
  # this is now the pad length instead of total length
  size <- size - length(vec)
  padding <- array(padding,size)

  if(side=="right"||side=="r")
  { vec <- c(vec,padding) }
  else if(side=="left"||side=="l")
  { vec <- c(padding,vec) }

  return(vec)
}

# row pad for data frames / matrices
rpad <- function(mat,size,padding=0,side="right")
{
  mat <- cbind(mat)
  size <- size - nrow(mat)
  COL <- ncol(mat)
  padding <- array(padding,c(size,COL))
  colnames(padding) <- colnames(mat)

  if(side=="right"||side=="r")
  { mat <- rbind(mat,padding) }
  else if(side=="left" || side=="l")
  { mat <- rbind(padding,mat) }

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
