# Score the mean distance of genes to the telomeres across species.
#
# A score will be given to each gene such that 0.0 corresponds to the maximal
# mean distance across all genes and 1.0 corresponds to a distance of 0.
proximity <- function(preset, progress = NULL) {
    # Prefilter distances by species and gene.
    distances <- geposan::distances[
        species %chin% preset$species_ids & gene %chin% preset$gene_ids
    ]

    # Compute the score as described above.

    distances <- distances[, .(mean_distance = mean(distance)), by = "gene"]
    max_distance <- distances[, max(mean_distance)]
    distances[, score := 1 - mean_distance / max_distance]

    if (!is.null(progress)) {
        # We do everything in one go, so it's not possible to report detailed
        # progress information. As the method is relatively quick, this should
        # not be a problem.
        progress(1.0)
    }

    distances[, .(gene, score)]
}
