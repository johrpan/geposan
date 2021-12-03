#' Plot gene positions.
#'
#' This function requires the package `plotly`.
#'
#' @param species_ids IDs of species to show in the plot.
#' @param gene_sets A list of gene sets (containing vectors of gene IDs) that
#'   will be highlighted in the plot. The names will be used as labels.
#'
#' @export
plot_positions <- function(species_ids, gene_sets) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
        stop("Please install \"plotly\" to use this function.")
    }

    # Prefilter data by species.
    data <- geposan::distances[species %chin% species_ids]

    species_max_distance <- data[,
        .(max_distance = max(distance)),
        by = species
    ]

    # Prefilter species.
    species <- geposan::species[id %chin% species_ids]

    plot <- plotly::plot_ly(colors = "Set2") |>
        plotly::layout(
            xaxis = list(
                title = "Species",
                tickvals = species$id,
                ticktext = species$name
            ),
            yaxis = list(title = "Distance to telomeres [Bp]"),
            bargap = 0.9
        ) |> plotly::add_bars(
            data = species_max_distance,
            x = ~species,
            y = ~max_distance,
            color = "All genes"
        )

    if (length(gene_sets) > 0) {
        # Include gene information which will be used for labeling
        gene_set_data <- merge(
            data[gene %chin% unlist(gene_sets)],
            geposan::genes,
            by.x = "gene",
            by.y = "id"
        )

        for (gene_set_name in names(gene_sets)) {
            gene_set <- gene_sets[[gene_set_name]]

            plot <- plot |> plotly::add_markers(
                data = gene_set_data[gene %chin% gene_set],
                x = ~species,
                y = ~distance,
                text = ~name,
                color = gene_set_name,
                marker = list(size = 10, opacity = 0.66)
            )
        }
    }

    plot
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
            xaxis = list(tickvals = names(rankings)),
            yaxis = list(title = "Score")
        )

    is_first <- TRUE

    for (ranking_name in names(rankings)) {
        ranking <- rankings[[ranking_name]]

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

        if (length(gene_sets) > 0) {
            gene_set_data <- merge(
                ranking[gene %chin% unlist(gene_sets)],
                geposan::genes,
                by.x = "gene",
                by.y = "id"
            )

            for (gene_set_name in names(gene_sets)) {
                gene_set <- gene_sets[[gene_set_name]]

                plot <- plot |> plotly::add_markers(
                    data = gene_set_data[gene %chin% gene_set],
                    x = ranking_name,
                    y = ~score,
                    text = ~name,
                    color = gene_set_name,
                    showlegend = is_first,
                    marker = list(size = 20, opacity = 0.66)
                )
            }
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
#' @param gene_sets A named list of gene sets (containing vectors of gene IDs)
#'   that will be highlighted in the plot. The names will be used in the legend.
#' @param max_rank The maximum rank of included genes. All genes that are ranked
#'   lower will appear greyed out.
#'
#' @seealso ranking()
#'
#' @export
plot_scores <- function(ranking, gene_sets = NULL, max_rank = NULL) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
        stop("Please install \"plotly\" to use this function.")
    }

    plot <- plotly::plot_ly(colors = "Set2") |>
        plotly::add_markers(
            data = ranking,
            x = ~rank,
            y = ~score,
            color = "All genes",
            hoverinfo = "skip"
        ) |>
        plotly::layout(
            xaxis = list(title = "Rank"),
            yaxis = list(title = "Score")
        )

    if (length(gene_sets) > 0) {
        # Include gene information which will be used for labeling
        gene_set_data <- merge(
            ranking[gene %chin% unlist(gene_sets)],
            geposan::genes,
            by.x = "gene",
            by.y = "id"
        )

        for (gene_set_name in names(gene_sets)) {
            gene_set <- gene_sets[[gene_set_name]]

            plot <- plot |> plotly::add_markers(
                data = gene_set_data[gene %chin% gene_set],
                x = ~rank,
                y = ~score,
                text = ~name,
                color = gene_set_name,
                marker = list(size = 20, opacity = 0.66)
            )
        }
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
#' @param gene_sets A named list of gene sets (containing vectors of gene IDs)
#'   that will be shown as separate boxes. The names will be used as labels.
#'
#' @seealso ranking()
#'
#' @export
plot_boxplot <- function(ranking, gene_sets = NULL) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
        stop("Please install \"plotly\" to use this function.")
    }

    plot <- plotly::plot_ly(colors = "Set2") |>
        plotly::add_boxplot(
            data = ranking,
            x = "All genes",
            y = ~score,
            color = "All genes",
            showlegend = FALSE
        ) |>
        plotly::layout(
            xaxis = list(tickvals = c("All genes", names(gene_sets))),
            yaxis = list(title = "Score")
        )

    if (length(gene_sets) > 0) {
        for (gene_set_name in names(gene_sets)) {
            gene_set <- gene_sets[[gene_set_name]]

            plot <- plot |> plotly::add_boxplot(
                data = ranking[gene %chin% gene_set],
                x = gene_set_name,
                y = ~score,
                color = gene_set_name,
                showlegend = FALSE
            )
        }
    }

    plot
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
