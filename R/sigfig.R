sigfig <- function(est,VAR=NULL,SD=NULL,level=0.95,digits=2,...)
{
  if(!is.null(VAR)) { SD <- sqrt(VAR) }

  if(is.null(SD)) # CIs
  {
    est <- rbind(est) # expecting (low,est,high) or table thereof
    NAMES <- rownames(est)

    DIFF <- cbind(est[,2]-est[,1],est[,3]-est[,2]) # sub-interval differences
    DIFF <- cbind(DIFF[,1], pmin(DIFF[,1],DIFF[,2]), DIFF[,2])
    POW <- floor(log10(DIFF) + .Machine$double.eps) - digits + 1 # final decimal place

    FIRST <- floor(log10(abs(est)))
    SIG <- FIRST-POW + 1
    ZERO <- SIG<1 # estimate is below sigfig threshold

    FORMAT <- function(i)
    {
      if(est[i]==Inf) { return("\u221E") }
      if(est[i]==-Inf) { return("-\u221E") }

      if(ZERO[i]) # will turn to zeroes
      { S <- DIFF[i] }
      else
      {
        S <- est[i]
        digits <- SIG[i]
      }

      S <- signif(S,digits=digits)
      S <- formatC(S,format="fg",digits=digits,flag="#")

      if(ZERO[i]) # turn S into zeroes
      {
        if(grepl("[eE]",S))
        { S <- strsplit(S,"[eE]") }

        S[1] <- gsub('[1-9]','0',S[1])

        S <- paste(S,collapse="e")
      }

      # not sure what whitespace is for
      S <- gsub(" ","",S)

      # strip trailing decimal
      n <- nchar(S)
      if(substr(S,n,n)==".") { S <- substr(S,1,n-1) }

      return(S)
    }

    est[] <- sapply(1:length(est),FORMAT)
    CI <- paste0("(",est[,1],"\u2014",est[,3],")")

    est <- paste(est[,2],CI)
  }
  else # mean and VAR
  {
    NAMES <- names(est)

    alpha <- 1-level
    z <- stats::qnorm(1-alpha/2)
    est <- cbind(est-z*SD,est,est+z*SD)
    colnames(est) <- NAMES.CI
    est <- sigfig(est)[,2]
  }

  names(est) <- NAMES
  return(est)
}
