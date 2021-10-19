#' Create a new preset.
#'
#' A preset is used to specify which methods and inputs should be used for an
#' analysis. Note that the genes to process should normally include the
#' reference genes to be able to assess the results later.
#'
#' Available methods are:
#'
#'  - `clusteriness` How much the gene distances cluster across species.
#'  - `correlation` The mean correlation with the reference genes.
#'  - `proximity` Mean proximity to telomeres.
#'  - `neural` Assessment by neural network.
#'
#' @param methods IDs of methods to apply.
#' @param species IDs of species to include.
#' @param genes IDs of genes to screen.
#' @param reference_genes IDs of reference genes to compare to.
#'
#' @return The preset to use with [analyze()].
#'
#' @export
preset <- function(methods, species, genes, reference_genes) {
    list(
        method_ids = methods,
        species_ids = species,
        gene_ids = genes,
        reference_gene_ids = reference_genes
    )
}

#' Analyze by applying the specified preset.
#'
#' @param preset The preset to use which can be created using [preset()].
#' @param progress A function to be called for progress information. The
#'   function should accept a number between 0.0 and 1.0 for the current
#'   progress.
#'
#' @return A [data.table] with one row for each gene identified by it's ID
#'   (`gene` column). The additional columns contain the resulting scores per
#'   method and are named after the method IDs.
#'
#' @export
analyze <- function(preset, progress = NULL) {
    # Available methods by ID.
    #
    # A method describes a way to perform a computation on gene distance data
    # that results in a single score per gene. The function should accept the
    # distances data, the preset to apply (see [preset()]) and an optional
    # progress function (that may be called with a number between 0.0 and 1.0)
    # as its parameters.
    #
    # The function should return a [data.table] with the following columns:
    #
    #  - `gene` Gene ID of the processed gene.
    #  - `score` Score for the gene between 0.0 and 1.0.
    methods <- list(
        "clusteriness" = clusteriness,
        "correlation" = correlation,
        "proximity" = proximity,
        "neural" = neural
    )

    total_progress <- 0.0
    method_count <- length(preset$method_ids)
    results <- data.table(gene = genes$id)

    for (method_id in preset$method_ids) {
        method_progress <- if (!is.null(progress)) function(p) {
            progress(total_progress + p / method_count)
        }

        method_results <- methods[[method_id]](
            distances,
            preset,
            method_progress
        )

        setnames(method_results, "score", method_id)

        results <- merge(
            results,
            method_results,
            by = "gene"
        )

        total_progress <- total_progress + 1 / method_count
    }

    if (!is.null(progress)) {
        progress(1.0)
    }

    results
}
