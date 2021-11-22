#' Analyze by applying the specified preset.
#'
#' @param preset The preset to use which should be created using [preset()].
#' @param progress A function to be called for progress information. The
#'   function should accept a number between 0.0 and 1.0 for the current
#'   progress.
#'
#' @returns An object containing the results of the analysis with the following
#'   items:
#'   \describe{
#'     \item{`preset`}{The preset that was used.}
#'     \item{`weights`}{The optimal weights for ranking the reference genes.}
#'     \item{`ranking`}{The optimal ranking created using the weights.}
#'   }
#'
#' @export
analyze <- function(preset, progress = NULL) {
    if (class(preset) != "geposan_preset") {
        stop("Preset is invalid. Use geposan::preset() to create one.")
    }

    # Available methods by ID.
    #
    # A method describes a way to perform a computation on gene distance data
    # that results in a single score per gene. The function should accept the
    # preset to apply (see [preset()]) and an optional progress function (that
    # may be called with a number between 0.0 and 1.0) as its parameters.
    #
    # The function should return a [data.table] with the following columns:
    #
    #  - `gene` Gene ID of the processed gene.
    #  - `score` Score for the gene between 0.0 and 1.0.
    methods <- list(
        "clusteriness" = clusteriness,
        "correlation" = correlation,
        "neural" = neural,
        "proximity" = proximity
    )

    analysis <- cached("analysis", preset, {
        total_progress <- 0.0
        method_count <- length(preset$methods)
        results <- data.table(gene = preset$gene_ids)

        for (method_id in preset$methods) {
            method_progress <- if (!is.null(progress)) {
                function(p) {
                    progress(total_progress + p / method_count)
                }
            }

            method_results <- methods[[method_id]](
                preset,
                progress = method_progress
            )

            setnames(method_results, "score", method_id)

            results <- merge(
                results,
                method_results,
                by = "gene"
            )

            total_progress <- total_progress + 1 / method_count
        }

        results <- structure(
            results,
            class = c("geposan_results", class(results))
        )

        weights <- optimal_weights(
            results,
            preset$methods,
            preset$reference_gene_ids,
            target = preset$optimization_target
        )

        ranking <- ranking(results, weights)

        structure(
            list(
                preset = preset,
                weights = weights,
                ranking = ranking
            ),
            class = "geposan_analysis"
        )
    })

    if (!is.null(progress)) {
        progress(1.0)
    }

    analysis
}
