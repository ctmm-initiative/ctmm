str.custom <- function(object,...)
{
  info <- object@info
  STR <- paste0("Formal class '",class(object)[1],"' [package \"ctmm\"]")

  # print $ slots
  NAMES <- names(object)
  object <- object@.Data
  names(object) <- NAMES
  STR <- c(STR,utils::capture.output(utils::str(object,...))[-1])

  # print @ slots
  END <- utils::capture.output(utils::str(info,...))
  END[1] <- "@ info"
  STR <- c(STR,END)

  # format and print
  cat(STR,sep="\n")
}

str.ctmm <- function(object,...) { str.custom(object,...) }

str.UERE <- function(object,...) { str.custom(object,...) }

# don't export this, the margins don't shift properly
str.covm <- function(object,...)
{
  STR <- "Formal class 'covm' [package \"ctmm\"]"
  STR <- c(STR,utils::capture.output(utils::str(object@.Data,...)))

  END <- utils::capture.output(utils::str(object@par,...))
  END[1] <- paste0("@ par       : ",END[1])
  END[2] <- paste0("              ",END[2])
  STR <- c(STR,END)

  END <- utils::capture.output(utils::str(object@isotropic,...))
  END <-    paste0("@ isotropic : ",END)
  STR <- c(STR,END)

  cat(STR,sep="\n")
}


