#' Score the mean distance of genes to the telomeres across species.
#'
#' A score will be given to each gene such that 0.0 corresponds to the maximal
#' mean distance across all genes and 1.0 corresponds to a distance of 0.
#'
#' @return An object of class `geposan_method`.
#'
#' @export
proximity <- function() {
    method(
        id = "proximity",
        name = "Proximity",
        description = "Proximity to telomeres",
        function(preset, progress) {
            species_ids <- preset$species_ids
            gene_ids <- preset$gene_ids

            cached("proximity", c(species_ids, gene_ids), {
                # Prefilter distances by species and gene.
                data <- geposan::distances[
                    species %chin% preset$species_ids &
                        gene %chin% preset$gene_ids
                ]

                # Compute the score as described above.
                data <- data[, .(mean_distance = mean(distance)), by = "gene"]
                max_distance <- data[, max(mean_distance)]
                data[, score := 1 - mean_distance / max_distance]

                progress(1.0)

                result(
                    method = "proximity",
                    scores = data[, .(gene, score)]
                )
            })
        }
    )
}
