#' Create a new preset.
#'
#' A preset is used to specify which methods and inputs should be used for an
#' analysis. Note that the genes to process should normally include the
#' reference genes to be able to assess the results later. The genes will be
#' filtered based on how many species have data for them. Genes which only have
#' orthologs for less than 25% of the input species will be excluded from the
#' preset and the analyis.
#'
#' Available methods are:
#'
#'  - `clusteriness` How much the gene distances to the nearest telomere
#'    cluster across species.
#'  - `clusteriness_positions` The same as `clusteriness` but using absolute
#'    gene positions instead of distances.
#'  - `correlation` The mean correlation of gene distances to the nearest
#'    telomere across species.
#'  - `correlation_positions` Correlation using position data.
#'  - `neural` Assessment by neural network trained using distances.
#'  - `neural_positions` Assessment by neural network trained using absolute
#'    position data.
#'  - `proximity` Mean proximity to telomeres.
#'
#' Available optimization targets are:
#'
#'  - `mean` Mean rank of the reference genes.
#'  - `max` First rank of the reference genes.
#'  - `min` Last rank of the reference genes.
#'
#' @param methods Methods to apply.
#' @param species_ids IDs of species to include.
#' @param gene_ids IDs of genes to screen.
#' @param reference_gene_ids IDs of reference genes to compare to.
#' @param optimization_target Parameter of the reference genes that the ranking
#'   should be optimized for.
#'
#' @return The preset to use with [analyze()].
#'
#' @export
preset <- function(methods = c(
                       "clusteriness",
                       "clusteriness_positions",
                       "correlation",
                       "correlation_positions",
                       "neural",
                       "proximity"
                   ),
                   species_ids = NULL,
                   gene_ids = NULL,
                   reference_gene_ids = NULL,
                   optimization_target = "mean_rank") {
    # Count included species per gene.
    genes_n_species <- geposan::distances[
        species %chin% species_ids,
        .(n_species = .N),
        by = "gene"
    ]

    # Filter out genes with less than 25% existing orthologs.
    gene_ids_filtered <- genes_n_species[
        n_species >= 0.25 * length(species_ids),
        gene
    ]

    # The included data gets sorted to be able to produce predictable hashes
    # for the object later.
    structure(
        list(
            methods = sort(methods),
            species_ids = sort(species_ids),
            gene_ids = sort(gene_ids_filtered),
            reference_gene_ids = sort(reference_gene_ids),
            optimization_target = optimization_target
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
    cat(x$methods, sep = ", ")

    cat(sprintf(
        "\n  Input data: %i species, %i genes",
        length(x$species_ids),
        length(x$gene_ids)
    ))

    cat(sprintf(
        "\n  Comparison data: %i reference genes",
        length(x$reference_gene_ids)
    ))

    cat(sprintf(
        "\n  Optimization target: %s\n",
        x$optimization_target
    ))

    invisible(x)
}
