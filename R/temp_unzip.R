# To decompress zip, gz, bzip2, xz into temp file, run function then remove temp file.
# more detail in https://gist.github.com/xhdong-umd/6429e7f96735142fa467f3b1daa91a2c
temp_unzip <- function(filename, fun, ...){
  BFR.SIZE <- 1e7
  if (!file.exists(filename)) {
    stop("No such file: ", filename);
  }
  if (!is.function(fun)) {
    stop(sprintf("Argument 'fun' is not a function: %s", mode(fun)));
  }
  temp_dir <- tempdir()
  # test if it's zip
  files_in_zip <- try(unzip(filename, list = TRUE)$Name, silent = TRUE)
  if (class(files_in_zip) == "character") {
    if(length(files_in_zip)>1) { warning(paste0(
      "  Zip file contains multiple files.\n  Mac OS built in zip compressor will add hidden folder for even single file zip.\nUsing the first file: ", 
      files_in_zip[1])) }
    unzip(filename, exdir = temp_dir, overwrite = TRUE)
    dest_file <- file.path(temp_dir, files_in_zip[1])
  } else {
    dest_file <- tempfile()
    # Setup input and output connections
    inn <- gzfile(filename, open = "rb")
    out <- file(description = dest_file, open = "wb")
    # Process
    nbytes <- 0
    repeat {
      bfr <- readBin(inn, what=raw(0L), size=1L, n=BFR.SIZE)
      n <- length(bfr)
      if (n == 0L) break;
      nbytes <- nbytes + n
      writeBin(bfr, con=out, size=1L)
      bfr <- NULL  # Not needed anymore
    }
    close(inn)
    close(out)
  }
  # call fun with temp file
  res <- fun(dest_file, ...)
  file.remove(dest_file)
  return(res)
}
