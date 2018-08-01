anonymize <- function(data)
{
  data <- listify(data)

  for(i in 1:length(data))
  {
    projection(data[[i]]) <- median(data[[i]],k=2)
    data[[i]]$t <- data[[i]]$t - data[[i]]$t[1]

    COLS <- c("timestamp","longitude","latitude")
    for(cl in COLS) { data[[i]][[cl]] <- NULL }

    SLOTS <- c("timezone","projection")
    for(sl in SLOTS) { attr(data[[i]],'info')[[sl]] <- NULL }
  }

  if(length(data)==1) { data <- data[[1]] }
  return(data)
}
