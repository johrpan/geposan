#' Plot a ranking as a scatter plot of scores.
#'
#' This function requires the package `plotly`.
#'
#' @param ranking The ranking to visualize.
#' @param gene_sets A list of gene sets (containing vectors of gene IDs) that
#'   will be highlighted in the plot.
#' @param labels Labels for the gene sets. This is required if gene sets are
#'   given and has to have the same length.
#'
#' @seealso ranking()
#'
#' @export
plot_scores <- function(ranking, gene_sets = NULL, labels = NULL) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
        stop("Please install \"plotly\" to use this function.")
    }

    plot <- plotly::plot_ly() |>
        plotly::add_trace(
            data = ranking,
            x = ~rank,
            y = ~score,
            color = "All genes",
            type = "scatter",
            mode = "markers",
            hoverinfo = "skip"
        ) |>
        plotly::layout(
            xaxis = list(title = "Rank"),
            yaxis = list(title = "Score")
        )

    if (length(gene_sets) > 0) {
        # Take out the genes to be highlighted.
        gene_set_data <- ranking[gene %chin% unlist(gene_sets)]

        # Add labels for each gene set.
        for (i in seq_along(gene_sets)) {
            gene_set_data[gene %chin% gene_sets[[i]], label := labels[i]]
        }

        # Include gene information which will be used for labeling
        gene_set_data <- merge(gene_set_data, genes, by.x = "gene", by.y = "id")

        plot <- plot |> plotly::add_trace(
            data = gene_set_data,
            x = ~rank,
            y = ~score,
            color = ~label,
            text = ~name,
            type = "scatter",
            mode = "markers",
            marker = list(size = 20)
        )
    }

    plot
}

#' Visualize a ranking by comparing gene sets in a boxplot.
#'
#' This function requires the package `plotly`.
#'
#' @param ranking The ranking to visualize.
#' @param gene_sets A list of gene sets (containing vectors of gene IDs) that
#'   will be shown as separate boxes.
#' @param labels Labels for the gene sets. This is required if gene sets are
#'   given and has to have the same length.
#'
#' @seealso ranking()
#'
#' @export
plot_boxplot <- function(ranking, gene_sets = NULL, labels = NULL) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
        stop("Please install \"plotly\" to use this function.")
    }

    data <- copy(ranking)

    # Add labels for each gene set.
    for (i in seq_along(gene_sets)) {
        data[gene %chin% gene_sets[[i]], label := labels[i]]
    }

    # Label the other genes.
    data[!gene %chin% unlist(gene_sets), label := "Other genes"]

    plotly::plot_ly(
        data = data,
        y = ~score,
        color = ~label,
        type = "box"
    )
}
