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
  files_in_zip <- try(utils::unzip(filename, list = TRUE)$Name, silent = TRUE)
  if (class(files_in_zip)[1] == "character") {
    # hidden files can be ignored: starting with ., ending with $, __MACOSX folder
    visible_files <- files_in_zip[!grepl("((^__MACOSX\\/.*)|(^\\..*)|(^.*\\$$))",
                                          files_in_zip)]
    # will not continue for multiple non-hidden files since behavior is not well defined.
    if(length(visible_files)>1) {
      stop(paste0("Zip file contains multiple visible files:\n",
                  paste0("    ", visible_files, collapse = "\n")))
    }
    if(length(visible_files) == 0) { stop("\n  No visible file found in Zip file")}
    # proceed with single non-hidden file
    utils::unzip(filename, files = visible_files[1], exdir = temp_dir, overwrite = TRUE)
    dest_file <- file.path(temp_dir, visible_files[1])
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
