#' Result of applying a method on gene position data.
#'
#' @param method_id ID of the method that produced this result.
#' @param scores A `data.frame` mapping gene IDs (`gene`) to computed scores
#'   between 0.0 and 1.0 (`score`).
#' @param details Optional details that may contain intermediate results as
#'   well as other information on the method application.
#'
#' @return An object of class `geposan_result`.
#'
#' @export
result <- function(method_id, scores, details = list()) {
    stopifnot(is.data.frame(scores) &
        c("gene", "score") %chin% colnames(scores))
    stopifnot(is.list(details))

    structure(
        list(
            method_id = method_id,
            scores = scores,
            details = details
        ),
        class = "geposan_result"
    )
}

#' Print a result object.
#'
#' @param x The result to print.
#' @param ... Other parameters.
#'
#' @seealso [result()]
#'
#' @export
print.geposan_result <- function(x, ...) {
    cat(sprintf(
        paste0(
            "geposan result:",
            "\n  Method: %s",
            "\n  Number of genes: %i",
            "\n  Available details: %s",
            "\n"
        ),
        x$method_id,
        nrow(x$scores),
        paste(names(x$details), collapse = ", ")
    ))

    invisible(x)
}
