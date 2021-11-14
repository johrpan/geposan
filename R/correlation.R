# Compute the mean correlation coefficient comparing gene distances with a set
# of reference genes.
correlation <- function(preset, use_positions = FALSE, progress = NULL) {
    species_ids <- preset$species_ids
    gene_ids <- preset$gene_ids
    reference_gene_ids <- preset$reference_gene_ids

    cached(
        "correlation",
        c(species_ids, gene_ids, reference_gene_ids, use_positions),
        { # nolint
            # Prefilter distances by species.
            distances <- geposan::distances[species %chin% species_ids]

            # Tranform data to get species as rows and genes as columns. We
            # construct columns per species, because it requires fewer
            # iterations, and transpose the table afterwards.

            data <- data.table(gene = gene_ids)

            # Make a column containing distance data for each species.
            for (species_id in species_ids) {
                species_data <- if (use_positions) {
                    setnames(distances[
                        species == species_id,
                        .(gene, position)
                    ], "position", "distance")
                } else {
                    distances[
                        species == species_id,
                        .(gene, distance)
                    ]
                }

                data <- merge(data, species_data, all.x = TRUE)
                setnames(data, "distance", species_id)
            }

            # Transpose to the desired format.
            data <- transpose(data, make.names = "gene")

            if (!is.null(progress)) progress(0.33)

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

            # Remove correlations between the reference genes themselves.
            for (reference_gene_id in reference_gene_ids) {
                column <- quote(reference_gene_id)
                results[gene == reference_gene_id, eval(column) := NA]
            }

            if (!is.null(progress)) progress(0.66)

            # Compute the final score as the mean of known correlation scores.
            # Negative correlations will correctly lessen the score, which will
            # be clamped to zero as its lower bound. Genes with no possible
            # correlations at all will be assumed to have a score of 0.0.

            compute_score <- function(scores) {
                score <- mean(scores, na.rm = TRUE)

                if (is.na(score) | score < 0.0) {
                    score <- 0.0
                }

                score
            }

            results[,
                score := compute_score(as.matrix(.SD)),
                .SDcols = reference_gene_ids,
                by = gene
            ]

            results[, .(gene, score)]
        }
    )
}
