# Perform a cluster analysis.
#
# This function will cluster the data using `hclust` and `cutree` (with the
# specified height). Every cluster with at least two members qualifies for
# further analysis. Clusters are then ranked based on their size in relation
# to the number of values. The return value is a final score between zero and
# one. Lower ranking clusters contribute less to this score.
clusteriness_priv <- function(data, height = 1000000) {
    n <- length(data)

    # Return a score of 0.0 if there is just one or no value at all.
    if (n < 2) {
        return(0.0)
    }

    # Cluster the data and compute the cluster sizes.

    tree <- stats::hclust(stats::dist(data))
    clusters <- stats::cutree(tree, h = height)
    cluster_sizes <- sort(tabulate(clusters), decreasing = TRUE)

    # Compute the "clusteriness" score.

    score <- 0.0

    for (i in seq_along(cluster_sizes)) {
        cluster_size <- cluster_sizes[i]

        if (cluster_size >= 2) {
            cluster_score <- cluster_size / n
            score <- score + cluster_score / i
        }
    }

    score
}

# Process genes clustering their distance to telomeres.
clusteriness <- function(preset, progress = NULL) {
    species_ids <- preset$species_ids
    gene_ids <- preset$gene_ids

    cached("clusteriness", c(species_ids, gene_ids), {
        results <- data.table(gene = gene_ids)

        # Prefilter the input data by species.
        distances <- geposan::distances[species %chin% species_ids]

        # Add an index for quickly accessing data per gene.
        setkey(distances, gene)

        genes_done <- 0
        genes_total <- length(gene_ids)

        # Perform the cluster analysis for one gene.
        compute <- function(gene_id) {
            data <- distances[gene_id, distance]
            score <- clusteriness_priv(data)

            if (!is.null(progress)) {
                genes_done <<- genes_done + 1
                progress(genes_done / genes_total)
            }

            score
        }

        results[, score := compute(gene), by = 1:nrow(results)]
    })
}
