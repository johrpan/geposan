#' Score genes based on their correlation with the reference genes.
#'
#' @param summarize A function for combining the different correlation
#'   coefficients into one metric. By default, [stats::median()] is used. Other
#'   suggested options include [max()] and [mean()].
#'
#' @return An object of class `geposan_method`.
#'
#' @export
correlation <- function(summarize = stats::median) {
  method(
    id = "correlation",
    name = "Correlation",
    description = "Correlation with reference genes",
    function(preset, progress) {
      species_ids <- preset$species_ids
      gene_ids <- preset$gene_ids
      reference_gene_ids <- preset$reference_gene_ids

      cached(
        "correlation",
        c(species_ids, gene_ids, reference_gene_ids, summarize),
        { # nolint
          # Prefilter distances by species.
          distances <- geposan::distances[species %chin% species_ids]

          # Tranform data to get species as rows and genes as columns.
          # We construct columns per species, because it requires
          # fewer iterations, and transpose the table afterwards.

          data <- data.table(gene = gene_ids)

          # Make a column containing distance data for each species.
          for (species_id in species_ids) {
            species_data <- distances[
              species == species_id,
              .(gene, distance)
            ]

            data <- merge(data, species_data, all.x = TRUE)
            setnames(data, "distance", species_id)
          }

          # Transpose to the desired format.
          data <- transpose(data, make.names = "gene")

          progress(0.33)

          # Take the reference data.
          reference_data <- data[, ..reference_gene_ids]

          # Perform the correlation between all possible pairs.
          results <- stats::cor(
            data[, ..gene_ids],
            reference_data,
            use = "pairwise.complete.obs",
            method = "spearman"
          )

          results <- data.table(results, keep.rownames = TRUE)
          setnames(results, "rn", "gene")

          # Remove correlations between the reference genes
          # themselves.
          for (reference_gene_id in reference_gene_ids) {
            column <- quote(reference_gene_id)
            results[gene == reference_gene_id, eval(column) := NA]
          }

          progress(0.66)

          # Combine the correlation coefficients.
          results[,
            max_correlation := as.double(summarize(stats::na.omit(
              # Convert the data.table subset into a
              # vector to get the correct na.omit
              # behavior.
              as.matrix(.SD)[1, ]
            ))),
            .SDcols = reference_gene_ids,
            by = gene
          ]

          # Normalize scores.
          results[
            ,
            score := (max_correlation - min(max_correlation)) /
              (max(max_correlation) - min(max_correlation))
          ]

          # Normalize scores.

          results[, .(gene, score)]

          result(
            method = "correlation",
            scores = results[, .(gene, score)],
            details = list(all_correlations = results)
          )
        }
      )
    }
  )
}
