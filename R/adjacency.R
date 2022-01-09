#' Score genes based on their proximity to the reference genes.
#'
#' @param estimate A function that will be used to summarize the distance
#'   values for each gene. By default, [median()] is used.
#' @param combination A function that will be used to combine the different
#'   distances to the reference genes. By default [min()] is used. That means
#'   the distance to the nearest reference gene will be scored.
#'
#' @return An object of class `geposan_method`.
#'
#' @export
adjacency <- function(estimate = stats::median, combination = min) {
    method(
        id = "adjacency",
        name = "Adjacency",
        description = "Adjacency to reference genes",
        function(preset, progress) {
            species_ids <- preset$species_ids
            gene_ids <- preset$gene_ids
            reference_gene_ids <- preset$reference_gene_ids

            cached(
                "adjacency",
                c(
                    species_ids,
                    gene_ids,
                    reference_gene_ids,
                    estimate,
                    combination
                ),
                { # nolint
                    # Filter distances by species and gene and summarize each
                    # gene's distance values using the estimation function.
                    data <- geposan::distances[
                        species %chin% species_ids & gene %chin% gene_ids,
                        .(distance = as.numeric(estimate(distance))),
                        by = gene
                    ]

                    # Compute the absolute value of the difference between the
                    # estimated distances of each gene to the reference genes.
                    compute_difference <- function(distance_value,
                                                   comparison_ids) {
                        differences <- data[
                            gene %chin% comparison_ids,
                            .(difference = abs(distance_value - distance))
                        ]

                        combination(differences$difference)
                    }

                    # Compute the differences to the reference genes.
                    data[
                        !gene %chin% reference_gene_ids,
                        difference := compute_difference(
                            distance,
                            reference_gene_ids
                        ),
                        by = gene
                    ]

                    progress(0.5)

                    # Exclude the reference gene itself when computing its
                    # difference.
                    data[
                        gene %chin% reference_gene_ids,
                        difference := compute_difference(
                            distance,
                            reference_gene_ids[reference_gene_ids != gene]
                        ),
                        by = gene
                    ]

                    # Compute the final score by normalizing the difference.
                    data[, score := 1 - difference / max(difference)]

                    progress(1.0)

                    result(
                        method = "adjacency",
                        scores = data[, .(gene, score)],
                        details = list(data = data)
                    )
                }
            )
        }
    )
}
