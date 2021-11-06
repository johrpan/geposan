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
#' @param methods Methods to apply.
#' @param species_ids IDs of species to include.
#' @param gene_ids IDs of genes to screen.
#' @param reference_gene_ids IDs of reference genes to compare to.
#'
#' @return The preset to use with [analyze()].
#'
#' @export
preset <- function(methods = c(
                       "clusteriness",
                       "correlation",
                       "neural",
                       "proximity"
                   ),
                   species_ids = NULL,
                   gene_ids = NULL,
                   reference_gene_ids = NULL) {
    # The included data gets sorted to be able to produce predictable hashes
    # for the object later.
    structure(
        list(
            methods = sort(methods),
            species_ids = sort(species_ids),
            gene_ids = sort(gene_ids),
            reference_gene_ids = sort(reference_gene_ids)
        ),
        class = "geposan_preset"
    )
}

#' S3 method to print a preset object.
#'
#' @param x The preset to print.
#' @param ... Other parameters.
#'
#' @seealso [preset()]
#'
#' @export
print.geposan_preset <- function(x, ...) {
    cat("geposan preset:")
    cat("\n  Included methods: ")
    cat(x$method_ids, sep = ", ")

    cat(sprintf(
        "\n  Input data: %i species, %i genes",
        length(x$species_ids),
        length(x$gene_ids)
    ))

    cat(sprintf(
        "\n  Comparison data: %i reference genes\n",
        length(x$reference_gene_ids)
    ))

    invisible(x)
}
