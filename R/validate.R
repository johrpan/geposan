#' Perform cross-validation for the ranking.
#'
#' This function reoptimizes the ranking leaving out one of the original
#' reference genes at a time.
#'
#' @param ranking The ranking to validate.
#' @param reference_gene_ids The reference gene IDs whose ranking should be
#'   validated.
#' @param method_ids IDs of the methods that were used.
#' @param progress An optional progress function that should accept a single
#'   value between 0.0 and 1.0 for progress information.
#'
#' @returns A validation object with the following items:
#'   \describe{
#'     \item{`validation`}{A `data.table` containing percentiles of the
#'       comparison genes from the original ranking as well as their validation.
#'     }
#'     \item{`mean_score`}{The mean score of the genes.}
#'     \item{`mean_percentile_original`}{The mean percentile of the genes in
#'       the original ranking.
#'     }
#'     \item{`mean_percentile_validation`}{The mean percentile of the genes
#'       when optimizing without themselves.
#'     }
#'     \item{`mean_error`}{The mean absolute error.}
#'   }
#'
#' @export
validate <- function(ranking, reference_gene_ids, method_ids, progress = NULL) {
    if (!inherits(ranking, "geposan_ranking")) {
        stop("Ranking is invalid. Use geposan::ranking().")
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

    progress_state <- 0.0
    progress_step <- 1.0 / length(reference_gene_ids)

    results <- ranking[gene %chin% reference_gene_ids, .(gene, percentile)]

    for (gene_id in reference_gene_ids) {
        included_gene_ids <- reference_gene_ids[
            reference_gene_ids != gene_id
        ]

        weights <- optimal_weights(
            ranking,
            method_ids,
            included_gene_ids
        )

        ranking_validation <- ranking(ranking, weights)

        results[
            gene == gene_id,
            percentile_validation := ranking_validation[
                gene == gene_id,
                percentile
            ]
        ]

        if (!is.null(progress)) {
            progress_state <- progress_state + progress_step
            progress(progress_state)
        }
    }

    results[, error := percentile - percentile_validation]
    setorder(results, error)

    structure(
        list(
            validation = results,
            mean_percentile_original = results[, mean(percentile)],
            mean_percentile_validation = results[, mean(percentile_validation)],
            mean_error = results[, mean(error)]
        ),
        class = "geposan_validation"
    )
}

#' S3 method to print a validation object.
#'
#' @param x The validation to print.
#' @param ... Other parameters.
#'
#' @seealso [validate()]
#'
#' @export
print.geposan_validation <- function(x, ...) {
    cat(sprintf(
        paste0(
            "geposan validation:",
            "\n  Mean percentile original: %.1f%%",
            "\n  Mean percentile validation: %.1f%%",
            "\n  Mean error: %.1f percent points",
            "\n"
        ),
        x$mean_percentile_original * 100,
        x$mean_percentile_validation * 100,
        x$mean_error * 100
    ))

    invisible(x)
}
