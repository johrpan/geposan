#' Score the distance of genes to the telomeres across species.
#'
#' A score will be given to each gene such that 0.0 corresponds to the maximal
#' distance across all genes and 1.0 corresponds to a distance of 0.
#'
#' @param id Unique ID for the method and its results.
#' @param name Human readable name for the method.
#' @param description Method description.
#' @param summarize A function for combining the different proximities into one
#'   metric. By default, [stats::median()] is used. Other suggested options
#'   include [min()] and [mean()].
#'
#' @return An object of class `geposan_method`.
#'
#' @export
distance <- function(id = "distance",
                     name = "Distance",
                     description = "Distance to telomeres",
                     summarize = stats::median) {
  method(
    id = id,
    name = name,
    description = description,
    function(preset, progress) {
      species_ids <- preset$species_ids
      gene_ids <- preset$gene_ids

      cached(id, c(species_ids, gene_ids), {
        # Prefilter distances by species and gene.
        data <- geposan::distances[
          species %chin% preset$species_ids &
            gene %chin% preset$gene_ids
        ]

        # Compute the score as described above.
        data <- data[,
          .(combined_distance = as.double(summarize(distance))),
          by = "gene"
        ]

        # Normalize scores.
        data[, score := combined_distance / max(combined_distance)]

        progress(1.0)

        result(
          method = "distance",
          scores = data[, .(gene, score)],
          details = list(data = data)
        )
      })
    }
  )
}
