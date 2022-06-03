# Cache the value of an expression on the file system.
#
# The expression will be evaluated if there is no matching cache file found.
# The cache files will be located in a directory "cache" located in the current
# working directory.
#
# @param name Human readable part of the cache file name.
# @param objects A vector of objects that this expression depends on. The hash
#   of those objects will be used for identifying the cache file.
cached <- function(name, objects, expr) {
  if (!dir.exists("cache")) {
    dir.create("cache")
  }

  id <- rlang::hash(objects)
  cache_file <- sprintf("cache/%s_%s.rda", name, id)

  if (file.exists(cache_file)) {
    # If the cache file exists, we restore the data from it.
    load(cache_file)
  } else {
    # If the cache file doesn't exist, we have to do the computation.
    data <- expr

    # The results are cached for the next run.
    save(data, file = cache_file, compress = "xz")
  }

  data
}

#' Format and round a numeric value.
#'
#' @param number The number to use.
#' @param digits Number of decimal places.
#'
#' @return A character value.
#' @noRd
num <- function(number, digits) {
  format(round(number, digits = digits), nsmall = digits)
}

# This is needed to make data.table's symbols available within the package.
#' @import data.table
NULL
