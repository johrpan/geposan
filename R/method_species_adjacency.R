#' Score genes based on their adjacency to the reference genes within species.
#'
#' For each gene and species, the method will first combine the gene's distances
#' to the reference genes within that species. Afterwards, the results are
#' summarized across species and determine the gene's score.
#'
#' @param id Unique ID for the method and its results.
#' @param name Human readable name for the method.
#' @param description Method description.
#' @param distance_estimate Function for combining the distance differences
#'   within one species.
#' @param summarize Function for summarizing the distance values across species.
#'
#' @return An object of class `geposan_method`.
#'
#' @seealso [adjacency()]
#'
#' @export
species_adjacency <- function(id = "species_adjacency",
                              name = "Species adj.",
                              description = "Species adjacency",
                              distance_estimate = stats::median,
                              summarize = stats::median) {
  method(
    id = id,
    name = name,
    description = description,
    function(preset, progress) {
      species_ids <- preset$species_ids
      gene_ids <- preset$gene_ids
      reference_gene_ids <- preset$reference_gene_ids

      cached(
        "species_adjacency",
        c(
          species_ids,
          gene_ids,
          reference_gene_ids,
          distance_estimate,
          summarize
        ),
        { # nolint
          # Prefilter distances.
          data <- geposan::distances[
            species %chin% species_ids & gene %chin% gene_ids
          ]

          progress_state <- 0.0
          progress_step <- 0.9 / length(species_ids)

          # Iterate through all species and find the distance
          # estimates within that species.
          for (species_id in species_ids) {
            # For all genes, compute the distance to one reference
            # gene at a time in one go.
            for (reference_gene_id in reference_gene_ids) {
              comparison_distance <- data[
                species == species_id &
                  gene == reference_gene_id,
                distance
              ]

              column <- quote(reference_gene_id)

              if (length(comparison_distance) != 1) {
                # If we don't have a comparison distance, we
                # can't compute a difference. This happens, if
                # the species doesn't have the reference gene.
                data[
                  species == species_id &
                    gene %chin% gene_ids,
                  eval(column) := NA_integer_
                ]
              } else {
                data[
                  species == species_id &
                    gene %chin% gene_ids,
                  eval(column) :=
                    abs(distance - comparison_distance)
                ]
              }
            }

            # Combine the distances to the different reference genes
            # into one value using the provided function.
            data[
              species == species_id &
                gene %chin% gene_ids,
              combined_distance := as.numeric(
                distance_estimate(stats::na.omit(
                  # Convert the data.table subset into a
                  # vector to get the correct na.omit
                  # behavior.
                  as.matrix(.SD)[1, ]
                ))
              ),
              .SDcols = reference_gene_ids,
              by = gene
            ]

            progress_state <- progress_state + progress_step
            progress(progress_state)
          }

          progress(0.9)

          # Remove the distances between the reference genes.
          for (reference_gene_id in reference_gene_ids) {
            column <- quote(reference_gene_id)
            data[gene == reference_gene_id, eval(column) := NA]
          }

          # Recompute the combined distance for the reference genes.
          data[
            gene %chin% reference_gene_ids,
            combined_distance := as.numeric(
              distance_estimate(stats::na.omit(
                as.matrix(.SD)[1, ]
              ))
            ),
            .SDcols = reference_gene_ids,
            by = list(species, gene)
          ]

          # Combine the distances into one value.
          results <- data[,
            .(
              summarized_distances = as.numeric(
                summarize(stats::na.omit(combined_distance))
              )
            ),
            by = gene
          ]

          # Compute the final score by normalizing the difference.
          results[
            ,
            score := 1 - summarized_distances /
              max(summarized_distances)
          ]

          progress(1.0)

          result(
            method = "species_adjacency",
            scores = results[, .(gene, score)],
            details = list(
              data = data,
              results = results
            )
          )
        }
      )
    }
  )
}
