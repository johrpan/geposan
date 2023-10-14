#' Describe a new method for analyzing gene position data.
#'
#' @param id Unique identifier for the method.
#' @param name Human readable name.
#' @param description Slightly longer description.
#' @param func Function to apply the method. The function should accept two
#'   parameters: an object of class `geposan_preset` as input and a function to
#'   report progress information to as a numeric value. The return value should
#'   be an object of class `geposan_result`.
#'
#' @return An object of class `geposan_method`.
#'
#' @export
method <- function(id, name, description, func) {
  stopifnot(is.character(id) & length(id) == 1)
  stopifnot(is.character(name) & length(name) == 1)
  stopifnot(is.character(description) & length(description) == 1)
  stopifnot(is.function(func))

  structure(
    list(
      id = glue::glue("geposan_method_{id}"),
      name = name,
      description = description,
      func = func
    ),
    class = "geposan_method"
  )
}

#' Get a list of all available methods.
#'
#' @export
all_methods <- function() {
  list(
    distance(),
    adjacency(),
    clustering(),
    correlation(),
    random_forest()
  )
}

#' Print a method object.
#'
#' @param x The method to print.
#' @param ... Other parameters.
#'
#' @seealso [method()]
#'
#' @export
print.geposan_method <- function(x, ...) {
  cat(sprintf(
    paste0(
      "geposan method:",
      "\n  Method ID: %s",
      "\n  Name: %s",
      "\n  Description: %s",
      "\n"
    ),
    x$id,
    x$name,
    x$description
  ))

  invisible(x)
}
