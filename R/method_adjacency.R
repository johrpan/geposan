#' Score genes based on their proximity to the reference genes.
#'
#' In this case, the distance data that is available for one gene is first
#' combined. The resulting value is compared to the reference genes and
#' determines the gene's score in relation to other genes.
#'
#' @param id Unique ID for the method and its results.
#' @param name Human readable name for the method.
#' @param description Method description.
#' @param distance_estimate A function that will be used to summarize the
#'   distance values for each gene. See [densest()] for the default
#'   implementation.
#'
#' @return An object of class `geposan_method`.
#'
#' @export
adjacency <- function(id = "adjacency",
                      name = "Adjacency",
                      description = "Adjacency to reference genes",
                      distance_estimate = densest) {
  method(
    id = id,
    name = name,
    description = description,
    help = paste0(
      "Adjacency to the reference genes across species. This method penalizes ",
      "genes that do not occur in the region typical for the reference genes, ",
      "without artificially defining a fixed boundary."
    ),
    function(preset, progress) {
      species_ids <- preset$species_ids
      gene_ids <- preset$gene_ids
      reference_gene_ids <- preset$reference_gene_ids

      cached(
        id,
        c(
          species_ids,
          gene_ids,
          reference_gene_ids,
          distance_estimate
        ),
        { # nolint
          # Filter distances by species and gene and summarize each
          # gene's distance values using the estimation function.
          data <- geposan::distances[
            species %chin% species_ids & gene %chin% gene_ids,
            .(distance = as.numeric(distance_estimate(distance))),
            by = gene
          ]

          # Compute the absolute value of the difference between the
          # estimated distances of each gene to the reference genes.
          compute_difference <- function(distance_values,
                                         comparison_ids) {
            comparison_distance <- data[
              gene %chin% comparison_ids,
              distance_estimate(distance)
            ]

            abs(distance_values - comparison_distance)
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
            )
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
