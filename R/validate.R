#' Perform cross-validation for the analysis.
#'
#' This function reoptimizes the analysis leaving out one of the original
#' reference genes at a time.
#'
#' @param analysis The analysis to validate.
#' @param progress An optional progress function that should accept a single
#'   value between 0.0 and 1.0 for progress information.
#'
#' @returns An object containing the mean absolute error and the mean percent
#'   rank for the original analysis as well as the validation.
#'
#' @export
validate <- function(analysis, progress = NULL) {
    if (class(analysis) != "geposan_analysis") {
        stop("Analysis is invalid. Use geposan::analyze().")
    }

    cached("validation", analysis$preset, {
        reference_gene_ids <- analysis$preset$reference_gene_ids
        n_references <- length(reference_gene_ids)
        methods <- analysis$preset$methods
        ranking_reference <- analysis$ranking
        n_ranks <- nrow(ranking_reference)

        mean_error_reference <- mean(
            1.0 - ranking_reference[gene %chin% reference_gene_ids, score]
        )

        mean_rank_reference <- mean(
            1.0 - ranking_reference[gene %chin% reference_gene_ids, rank] /
                n_ranks
        )

        mean_error_validation <- 0.0
        mean_rank_validation <- 0.0

        progress_state <- 0.0
        progress_step <- 1.0 / n_references

        for (validation_gene_id in reference_gene_ids) {
            included_gene_ids <- reference_gene_ids[
                reference_gene_ids != validation_gene_id
            ]

            weights <- optimal_weights(
                analysis,
                methods,
                included_gene_ids
            )

            ranking_validation <- ranking(analysis, weights)

            mean_error_validation <- mean_error_validation +
                (1.0 - ranking_validation[gene == validation_gene_id, score]) /
                    n_references

            mean_rank_validation <- mean_rank_validation +
                (1.0 - ranking_validation[gene == validation_gene_id, rank] /
                    n_ranks) / n_references

            if (!is.null(progress)) {
                progress_state <- progress_state + progress_step
                progress(progress_state)
            }
        }

        structure(
            list(
                mean_error_reference = mean_error_reference,
                mean_error_validation = mean_error_validation,
                mean_rank_reference = mean_rank_reference,
                mean_rank_validation = mean_rank_validation
            ),
            class = "geposan_validation"
        )
    })
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
    cat("geposan validation:\n")
    cat(sprintf(
        paste0(
            "\n  Absolute scores:",
            "\n  Mean error reference:  %.3f",
            "\n  Mean error validation: %.3f",
            "\n",
            "\n  Ranks:",
            "\n  Mean rank reference:   %.1f%%",
            "\n  Mean rank validation:  %.1f%%",
            "\n"
        ),
        x$mean_error_reference,
        x$mean_error_validation,
        x$mean_rank_reference * 100,
        x$mean_rank_validation * 100
    ))

    invisible(x)
}
