cleave <- function(object,fraction=0.5,name="CLEFT",...)
{
  C1 <- 'blue'
  C2 <- 'red'

  data <- object
  n <- nrow(data)
  SUB <- object$t[1] + fraction*(last(object$t)-first(object$t))
  SUB <- last( which(SUB > object$t) )
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
        plot(DATA,col=c(C1,C2),...)
        TITLE1 <- c( paste0("1:",SUB) , as.character(last(DATA[[1]]$timestamp)) )
        TITLE2 <- c( paste0(SUB+1,":",n) , as.character(DATA[[2]]$timestamp[1]) )

        if(SUB==0)
        {
          graphics::title(paste0(TITLE2[1],),col.main=C2,line=1)
          graphics::title(paste0(TITLE2[2],),col.main=C2)
        }
        else if(SUB==n)
        {
          graphics::title(TITLE1[1],col.main=C1,line=1)
          graphics::title(TITLE1[2],col.main=C1)
        }
        else
        {
          t <- c( last(DATA[[1]]$timestamp) , first(DATA[[2]]$timestamp) )
          t <- t[1] + diff(t)/2
          t <- as.character(t)

          # this was amazingly annoying to figure out
          eval( bquote(graphics::title(expression(.(TITLE1[1]) * "   " * phantom(.(TITLE2[1]))),col.main=C1,line=1)) )
          eval( bquote(graphics::title(expression(.(TITLE1[2]) * "   " * phantom(.(TITLE2[2]))),col.main=C1)) )
          eval( bquote(graphics::title(expression(phantom(.(TITLE1[1])) * "   " * .(TITLE2[1])),col.main=C2,line=1)) )
          eval( bquote(graphics::title(expression(phantom(.(TITLE1[2])) * "   " * .(TITLE2[2])),col.main=C2)) )
        }
      }
    },
    SUB=manipulate::slider(0,n,initial=SUB,step=1,label="Index",ticks=TRUE),
    STORE=manipulate::button(paste("Save to",name))
    )
}
