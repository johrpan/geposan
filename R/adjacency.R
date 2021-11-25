# Score genes based on their proximity to the reference genes.
#
# This method finds the distance value with the maximum density for each gene
# (i.e. the mode of its estimated distribution). Genes are scored by comparing
# those distance values with the values of the reference genes.
adjacency <- function(preset, progress = NULL) {
    species_ids <- preset$species_ids
    gene_ids <- preset$gene_ids
    reference_gene_ids <- preset$reference_gene_ids

    cached("adjacency", c(species_ids, gene_ids, reference_gene_ids), {
        # Get the virtual distance value with the highest density.
        compute_densest_distance <- function(distances) {
            if (length(distances <= 2)) {
                mean(distances)
            } else {
                d <- stats::density(distances)
                d$x[which.max(d$y)]
            }
        }

        # Filter distances by species and gene and find the distance with the
        # highest density of values for each gene.
        data <- geposan::distances[
            species %chin% species_ids & gene %chin% gene_ids,
            .(densest_distance = compute_densest_distance(distance)),
            by = gene
        ]

        # Compute the absolute value of the difference between the provided
        # densest distance value in comparison to the mean of the densest
        # distances of the comparison genes.
        compute_difference <- function(densest_distance, comparison_ids) {
            # Get the mean of the densest distances of the reference genes.
            mean_densest_distance <- data[
                gene %chin% comparison_ids,
                mean(densest_distance)
            ]

            abs(densest_distance - mean_densest_distance)
        }

        # Compute the differences to the reference genes.
        data[
            !gene %chin% reference_gene_ids,
            difference := compute_difference(
                densest_distance,
                reference_gene_ids
            )
        ]

        if (!is.null(progress)) {
            progress(0.5)
        }

        # Exclude the reference gene itself when computing its difference.
        data[
            gene %chin% reference_gene_ids,
            difference := compute_difference(
                densest_distance,
                reference_gene_ids[reference_gene_ids != gene]
            ),
            by = gene
        ]

        # Compute the final score by normalizing the difference.
        data[, score := 1 - difference / max(difference)]

        if (!is.null(progress)) {
            progress(1.0)
        }

        structure(
            list(
                results = data[, .(gene, score)],
                details = data
            ),
            class = "geposan_method_results"
        )
    })
}
