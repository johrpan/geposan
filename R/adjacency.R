#' Find the densest value in the data.
#'
#' This function assumes that data represents a continuous variable and finds
#' a single value with the highest estimated density. This can be used to
#' estimate the mode of the data. If there is only one value that value is
#' returned. If multiple density maxima with the same density exist, their mean
#' is returned.
#'
#' @param data The input data.
#'
#' @return The densest value of data.
#'
#' @export
densest <- function(data) {
    as.numeric(if (length(data) <= 0) {
        NULL
    } else if (length(data) == 1) {
        data
    } else {
        density <- stats::density(data)
        mean(density$x[density$y == max(density$y)])
    })
}

#' Score genes based on their proximity to the reference genes.
#'
#' @param estimate A function that will be used to summarize the distance
#'   values for each gene. See [densest()] for the default implementation.
#'
#' @return An object of class `geposan_method`.
#'
#' @export
adjacency <- function(estimate = densest) {
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
                c(species_ids, gene_ids, reference_gene_ids, estimate),
                { # nolint
                    # Filter distances by species and gene and summarize each
                    # gene's distance values using the estimation function.
                    data <- geposan::distances[
                        species %chin% species_ids & gene %chin% gene_ids,
                        .(distance = estimate(distance)),
                        by = gene
                    ]

                    # Compute the absolute value of the difference between the
                    # estimated distances of each gene to the reference genes.
                    compute_difference <- function(distance,
                                                   comparison_ids) {
                        reference_distance <- data[
                            gene %chin% comparison_ids,
                            mean(distance)
                        ]

                        abs(distance - reference_distance)
                    }

                    # Compute the differences to the reference genes.
                    data[
                        !gene %chin% reference_gene_ids,
                        difference := compute_difference(
                            distance,
                            reference_gene_ids
                        )
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
