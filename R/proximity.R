# Score the mean distance of genes to the telomeres across species.
#
# A score will be given to each gene such that 0.0 corresponds to the maximal
# mean distance across all genes and 1.0 corresponds to a distance of 0.
proximity <- function(distances, preset) {
    # Prefilter distances by species and gene.
    distances <- distances[
        species %chin% preset$species_ids & gene %chin% preset$gene_ids
    ]

    # Compute the score as described above.

    distances <- distances[, .(mean_distance = mean(distance)), by = "gene"]
    max_distance <- distances[, max(mean_distance)]
    distances[, score := 1 - mean_distance / max_distance]

    distances[, .(gene, score)]
}
