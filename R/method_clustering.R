#' Perform a cluster analysis.
#'
#' This function will cluster the data using [stats::hclust()] and
#' [stats::cutree()]. Every cluster with at least two members qualifies for
#' further analysis. Clusters are then ranked based on their size in relation
#' to the total number of values. The return value is a final score between
#' 0.0 and 1.0. Lower ranking clusters contribute less to this score.
#'
#' @param data The values that should be scored.
#' @param span The maximum span of values considered to be in one cluster.
#' @param weight The weight that will be given to the next largest cluster in
#'   relation to the previous one. For example, if `weight` is 0.5 (the
#'   default), the first cluster will weigh 1.0, the second 0.5, the third 0.25
#'   etc.
#'
#' @return A score between 0.0 and 1.0 summarizing how much the data clusters.
#'
#' @export
clusteriness <- function(data, span = 1000000, weight = 0.5) {
  n <- length(data)

  # Return a score of 0.0 if there is just one or no value at all.
  if (n < 2) {
    return(0.0)
  }

  # Cluster the data and compute the cluster sizes.

  tree <- stats::hclust(stats::dist(data))
  clusters <- stats::cutree(tree, h = span)
  cluster_sizes <- sort(tabulate(clusters), decreasing = TRUE)

  # Compute the "clusteriness" score.

  score <- 0.0

  for (i in seq_along(cluster_sizes)) {
    cluster_size <- cluster_sizes[i]

    if (cluster_size >= 2) {
      cluster_score <- cluster_size / n
      score <- score + weight^(i - 1) * cluster_score # nolint
    }
  }

  score
}

#' Process genes clustering their distance to telomeres.
#'
#' The result will be cached and can be reused for different presets, because
#' it is independent of the reference genes in use. Most parameters are exposed
#' for the [clusteriness()] function. See its documentation for more
#' information.
#'
#' @param id Unique ID for the method and its results.
#' @param name Human readable name for the method.
#' @param description Method description.
#' @param span See [clusteriness()].
#' @param weight See [clusteriness()].
#'
#' @return An object of class `geposan_method`.
#'
#' @seealso [clusteriness()]
#'
#' @export
clustering <- function(id = "clustering",
                       name = "Clustering",
                       description = "Clustering of genes",
                       span = 1000000,
                       weight = 0.5) {
  method(
    id = id,
    name = name,
    description = description,
    help = paste0(
      "Proportion of orthologs that have similar telomeric distances across ",
      "species. This favors genes whose position is evolutionarily conserved."
    ),
    function(preset, progress) {
      species_ids <- preset$species_ids
      gene_ids <- preset$gene_ids

      cached(
        id,
        c(species_ids, gene_ids, span, weight),
        { # nolint
          scores <- data.table(gene = gene_ids)

          # Prefilter the input data by species.
          distances <- geposan::distances[species %chin% species_ids]

          genes_done <- 0
          genes_total <- length(gene_ids)

          # Perform the cluster analysis for one gene.
          compute <- function(gene_id) {
            data <- distances[gene == gene_id, distance]

            score <- clusteriness(data, span = span, weight = weight)

            genes_done <<- genes_done + 1
            progress(genes_done / genes_total)

            score
          }

          scores[, score := compute(gene), by = gene]

          result(
            method = "clustering",
            scores = scores
          )
        }
      )
    }
  )
}
