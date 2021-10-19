# Compute the mean correlation coefficient comparing gene distances with a set
# of reference genes.
correlation <- function(distances, preset, progress = NULL) {
    results <- data.table(gene = preset$gene_ids)
    reference_gene_ids <- preset$reference_gene_ids
    reference_count <- length(reference_gene_ids)

    # Prefilter distances by species.
    distances <- distances[species %chin% preset$species_ids]

    # Add an index for quickly accessing data per gene.
    setkey(distances, gene)

    # Prepare the reference genes' data.
    reference_distances <- distances[gene %chin% reference_gene_ids]

    genes_done <- 0
    genes_total <- length(preset$gene_ids)

    # Perform the correlation for one gene.
    compute <- function(gene_id) {
        gene_distances <- distances[gene_id]
        gene_species_count <- nrow(gene_distances)

        # Return a score of 0.0 if there is just one or no value at all.
        if (gene_species_count <= 1) {
            return(0.0)
        }

        # Buffer for the sum of correlation coefficients.
        correlation_sum <- 0

        # Correlate with all reference genes but not with the gene itself.
        gene_reference_gene_ids <- reference_gene_ids[
            reference_gene_ids != gene_id
        ]

        for (reference_gene_id in gene_reference_gene_ids) {
            data <- merge(
                gene_distances,
                reference_distances[reference_gene_id],
                by = "species"
            )

            # Skip this reference gene, if there are not enough value pairs.
            # This will lessen the final score, because it effectively
            # represents a correlation coefficient of 0.0.
            if (nrow(data) <= 1) {
                next
            }

            # Order data by the reference gene's distance to get a monotonic
            # relation.
            setorder(data, distance.y)

            correlation <- abs(stats::cor(
                data[, distance.x], data[, distance.y],
                method = "spearman"
            ))

            # If the correlation is NA, this will effectively mean 0.0.
            if (!is.na(correlation)) {
                correlation_sum <- correlation_sum + correlation
            }
        }

        # Compute the score as the mean correlation coefficient.
        score <- correlation_sum / length(gene_reference_gene_ids)

        if (!is.null(progress)) {
            genes_done <<- genes_done + 1
            progress(genes_done / genes_total)
        }

        score
    }

    results[, score := compute(gene), by = 1:nrow(results)]
}
