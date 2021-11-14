# Score the mean distance of genes to the telomeres across species.
#
# A score will be given to each gene such that 0.0 corresponds to the maximal
# mean distance across all genes and 1.0 corresponds to a distance of 0.
proximity <- function(preset, progress = NULL) {
    species_ids <- preset$species_ids
    gene_ids <- preset$gene_ids

    cached("proximity", c(species_ids, gene_ids), {
        # Prefilter distances by species and gene.
        data <- geposan::distances[
            species %chin% preset$species_ids & gene %chin% preset$gene_ids
        ]

        # Compute the score as described above.
        data <- data[, .(mean_distance = mean(distance)), by = "gene"]
        max_distance <- data[, max(mean_distance)]
        data[, score := 1 - mean_distance / max_distance]

        if (!is.null(progress)) {
            # We do everything in one go, so it's not possible to report
            # detailed progress information. As the method is relatively quick,
            # this should not be a problem.
            progress(1.0)
        }

        data[, .(gene, score)]
    })
}
