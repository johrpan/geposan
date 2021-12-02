#' Plot gene positions.
#'
#' This function requires the package `plotly`.
#'
#' @param species_ids IDs of species to show in the plot.
#' @param gene_sets A list of gene sets (containing vectors of gene IDs) that
#'   will be highlighted in the plot.
#' @param labels Labels for the gene sets. This is required if gene sets are
#'   given and has to have the same length.
#' @param use_positions Whether to display positions instead of distances.
#'
#' @export
plot_positions <- function(species_ids,
                           gene_sets,
                           labels,
                           use_positions = FALSE) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
        stop("Please install \"plotly\" to use this function.")
    }

    data <- merge(
        geposan::distances[gene %chin% unlist(gene_sets) &
            species %chin% species_ids],
        geposan::genes[, .(id, name)],
        by.x = "gene", by.y = "id"
    )

    if (use_positions) {
        data[, value := position]
    } else {
        data[, value := distance]
    }

    # Add labels for each gene set.
    for (i in seq_along(gene_sets)) {
        data[gene %chin% gene_sets[[i]], label := labels[i]]
    }

    species <- geposan::species[id %chin% species_ids]

    yaxis_title <- if (use_positions) {
        "Position [Bp]"
    } else {
        "Distance to telomeres [Bp]"
    }

    plotly::plot_ly(
        data = data,
        x = ~species,
        y = ~value,
        color = ~label,
        text = ~name,
        type = "scatter",
        mode = "markers"
    ) |> plotly::layout(
        xaxis = list(
            title = "Species",
            tickvals = species$id,
            ticktext = species$name
        ),
        yaxis = list(title = yaxis_title)
    )
}


#' Plot a side-by-side comparison of multiple rankings.
#'
#' Each ranking's scores will be shown as a vertical violin plot without any
#' additional markings. The gene sets will be shown as markers on top of the
#' density visualization.
#'
#' This function requires the package `plotly`.
#'
#' @param rankings A named list of rankings to display. The names will be shown
#'   as labels in the plot.
#' @param gene_sets A named list of vectors of gene IDs to highlight. The names
#'   will be used to distinguish the sets and in the legend.
#'
#' @export
plot_rankings <- function(rankings, gene_sets) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
        stop("Please install \"plotly\" to use this function.")
    }

    plot <- plotly::plot_ly(colors = "Set2") |>
        plotly::layout(
            xaxis = list(
                title = "Ranking",
                tickvals = names(rankings)
            ),
            yaxis = list(title = "Score")
        )

    is_first <- TRUE

    for (ranking_name in names(rankings)) {
        ranking <- rankings[[ranking_name]]

        data <- merge(
            ranking,
            geposan::genes,
            by.x = "gene",
            by.y = "id"
        )

        plot <- plot |> plotly::add_trace(
            data = ranking,
            x = ranking_name,
            y = ~score,
            color = "All genes",
            type = "violin",
            spanmode = "hard",
            points = FALSE,
            showlegend = is_first,
            hoverinfo = "skip"
        )

        for (gene_set_name in names(gene_sets)) {
            gene_set <- gene_sets[[gene_set_name]]

            plot <- plot |> plotly::add_markers(
                data = data[gene %chin% gene_set],
                x = ranking_name,
                y = ~score,
                text = ~name,
                color = gene_set_name,
                showlegend = is_first,
                marker = list(size = 20, opacity = 0.66)
            )
        }

        is_first <- FALSE
    }

    plot
}


#' Plot a ranking as a scatter plot of scores.
#'
#' This function requires the package `plotly`.
#'
#' @param ranking The ranking to visualize.
#' @param gene_sets A list of gene sets (containing vectors of gene IDs) that
#'   will be highlighted in the plot.
#' @param labels Labels for the gene sets. This is required if gene sets are
#'   given and has to have the same length.
#' @param max_rank The maximum rank of the highlighted genes. All genes that
#'   are ranked lower will appear greyed out.
#'
#' @seealso ranking()
#'
#' @export
plot_scores <- function(ranking,
                        gene_sets = NULL,
                        labels = NULL,
                        max_rank = NULL) {
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
        gene_set_data <- merge(
            gene_set_data,
            geposan::genes,
            by.x = "gene",
            by.y = "id"
        )

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


    if (!is.null(max_rank)) {
        first_not_included_rank <- max_rank + 1
        last_rank <- ranking[, .N]

        if (first_not_included_rank <= last_rank) {
            plot <- plot |> plotly::layout(
                shapes = list(
                    type = "rect",
                    fillcolor = "black",
                    opacity = 0.1,
                    x0 = first_not_included_rank,
                    x1 = last_rank,
                    y0 = 0.0,
                    y1 = 1.0
                )
            )
        }
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

#' Show the distribution of scores across chromosomes.
#'
#' This function requires the package `plotly`.
#'
#' @param ranking The ranking to visualize.
#'
#' @seealso ranking()
#'
#' @export
plot_chromosomes <- function(ranking) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
        stop("Please install \"plotly\" to use this function.")
    }

    data <- merge(ranking, geposan::genes, by.x = "gene", by.y = "id")
    data <- data[, .(score = mean(score)), by = "chromosome"]

    # Get an orderable integer from a chromosome name.
    chromosome_index <- function(chromosome) {
        index <- suppressWarnings(as.integer(chromosome))

        ifelse(
            !is.na(index),
            index,
            ifelse(
                chromosome == "X",
                998,
                999
            )
        )
    }

    data[, index := chromosome_index(chromosome)]
    setorder(data, "index")

    plotly::plot_ly(
        data = data,
        x = ~chromosome,
        y = ~score,
        type = "bar"
    ) |>
        plotly::layout(
            xaxis = list(
                title = "Chromosome",
                categoryorder = "array",
                categoryarray = ~chromosome
            ),
            yaxis = list(title = "Mean score")
        )
}
