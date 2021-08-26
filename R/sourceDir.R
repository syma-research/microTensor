# source all files in a directory
sourceDir <- function (path, pattern = "[.][Rr]$", trace = TRUE, recursive = FALSE, 
          ...) 
{
  for (nm in list.files(path, pattern = pattern, recursive = recursive, 
                        ...)) {
    source(file.path(path, nm), ...)
    if (trace) 
      cat(nm, "\n")
  }
}