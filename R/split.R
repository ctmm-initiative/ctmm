cleave <- function(object,fraction=0.5,name="CLEFT",...)
{
  data <- object
  n <- nrow(data)
  SUB <- round(n*fraction)
  STORE <- FALSE

  manipulate::manipulate(
    {
      DATA <- list()
      DATA[[1]] <- if(SUB>=1) { data[1:SUB,] } else { data[NULL,] }
      DATA[[2]] <- if(SUB<n) { data[(SUB+1):n,] } else { data[NULL,] }
      names(DATA) <- c("before","after")

      if(STORE) # store to global variable `name`
      {
        envir <- .GlobalEnv
        assign(name,DATA,envir=envir)
        # message(paste("First segment ends at index",SUB))
      }
      else # don't store, update plot
      {
        plot(DATA,col=c('red','blue'),...)
        TITLE1 <- paste0("1:",SUB)
        TITLE2 <- paste0(SUB+1,":",n)

        if(SUB==0) { title(TITLE2,col.main="blue") }
        else if(SUB==n) { title(TITLE1,col.main="red") }
        else
        {
          # this was amazingly annoying to figure out
          eval( bquote(title(expression(.(TITLE1) * "   " * phantom(.(TITLE2))),col.main="red")) )
          eval( bquote(title(expression(phantom(.(TITLE1)) * "   " * .(TITLE2)),col.main="blue")) )
        }
      }
    },
    SUB=manipulate::slider(0,n,initial=SUB,step=1,label="Index",ticks=TRUE),
    STORE=manipulate::button(paste("Save to",name))
    )
}
