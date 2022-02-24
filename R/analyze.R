#' Analyze genes based on position data.
#'
#' @param preset The preset to use which should be created using [preset()].
#' @param progress A function to be called for progress information. The
#'   function should accept a number between 0.0 and 1.0 for the current
#'   progress. If no function is provided, a simple text progress bar will be
#'   shown.
#' @param include_results Whether to include the detailed results. If this is
#'   set to `FALSE`, only the `scores` are available.
#'
#' @returns An object containing the results of the analysis with the following
#'   items:
#'   \describe{
#'     \item{`preset`}{The preset that was used.}
#'     \item{`scores`}{Table containing all scores for each gene.}
#'     \item{`results`}{Results from the different methods including details.}
#'   }
#'
#' @export
analyze <- function(preset, progress = NULL, include_results = TRUE) {
    if (!inherits(preset, "geposan_preset")) {
        stop("Preset is invalid. Use geposan::preset() to create one.")
    }

    if (is.null(progress)) {
        progress_bar <- progress::progress_bar$new()
        progress_bar$update(0.0)

        progress <- function(progress_value) {
            if (!progress_bar$finished) {
                progress_bar$update(progress_value)
                if (progress_value >= 1.0) {
                    progress_bar$terminate()
                }
            }
        }
    }

    progress_buffer <- 0.0
    method_count <- length(preset$methods)

    method_progress <- function(progress_value) {
        progress(progress_buffer + progress_value / method_count)
    }

    scores <- data.table(gene = preset$gene_id)
    results <- list()

    for (method in preset$methods) {
        method_results <- method$func(preset, method_progress)

        scores <- merge(scores, method_results$scores, by = "gene")
        setnames(scores, "score", method$id)

        results <- c(results, list(method_results))

        progress_buffer <- progress_buffer + 1 / method_count
        progress(progress_buffer)
    }

    structure(
        list(
            preset = preset,
            scores = scores,
            results = if (include_results) results else NULL
        ),
        class = "geposan_analysis"
    )
}

#' Print an analysis object.
#'
#' @param x The analysis to print.
#' @param ... Other parameters.
#'
#' @seealso [analyze()]
#'
#' @export
print.geposan_analysis <- function(x, ...) {
    cat("geposan analysis:\n\n")
    print(x$preset)
    cat("\n")

    for (result in x$results) {
        print(result)
        cat("\n")
    }

    invisible(x)
}
