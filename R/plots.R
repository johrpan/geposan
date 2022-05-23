#' Base color for the plots.
#' @noRd
base_color <- function() "#1964bf"

#' Transparent version of the base color.
#' @noRd
base_color_transparent <- function() "#1964bf80"

#' Color palette for gene sets.
#' @noRd
gene_set_color <- function(index) {
    c("#FF7F00", "#4DAF4A", "#984EA3")[index]
}

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

    plot <- plotly::plot_ly() |>
        plotly::layout(
            xaxis = list(
                title = "Species",
                tickvals = species$id,
                ticktext = species$name
            ),
            yaxis = list(title = "Distance to telomeres [Bp]"),
            bargap = 0.9
        ) |>
        plotly::add_bars(
            data = species_max_distance,
            x = ~species,
            y = ~max_distance,
            name = "All genes",
            marker = list(color = base_color())
        )

    if (length(gene_sets) > 0) {
        # Include gene information which will be used for labeling
        gene_set_data <- merge(
            data[gene %chin% unlist(gene_sets)],
            geposan::genes,
            by.x = "gene",
            by.y = "id"
        )

        index <- 1

        for (gene_set_name in names(gene_sets)) {
            gene_set <- gene_sets[[gene_set_name]]

            plot <- plot |> plotly::add_markers(
                data = gene_set_data[gene %chin% gene_set],
                x = ~species,
                y = ~distance,
                name = gene_set_name,
                text = ~ glue::glue(
                    "<b>{name}</b><br>",
                    "{round(distance / 1000000, digits = 2)} MBp"
                ),
                hoverinfo = "text",
                marker = list(
                    size = 10,
                    color = gene_set_color(index)
                )
            )

            index <- index + 1
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

    plot <- plotly::plot_ly() |>
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
            name = "All genes",
            type = "violin",
            spanmode = "hard",
            points = FALSE,
            showlegend = is_first,
            hoverinfo = "skip",
            line = list(color = base_color()),
            fillcolor = base_color_transparent()
        )

        if (length(gene_sets) > 0) {
            gene_set_data <- merge(
                ranking[gene %chin% unlist(gene_sets)],
                geposan::genes,
                by.x = "gene",
                by.y = "id"
            )

            index <- 1

            for (gene_set_name in names(gene_sets)) {
                gene_set <- gene_sets[[gene_set_name]]

                plot <- plot |> plotly::add_markers(
                    data = gene_set_data[gene %chin% gene_set],
                    x = ranking_name,
                    y = ~score,
                    name = gene_set_name,
                    text = ~ glue::glue(
                        "<b>{name}</b><br>",
                        "Score: {round(score, digits = 2)}<br>",
                        "Rank: {rank}<br>",
                        "Percentile: {round(percentile * 100, digits = 2)}%"
                    ),
                    hoverinfo = "text",
                    showlegend = is_first,
                    marker = list(
                        size = 10,
                        color = gene_set_color(index)
                    )
                )

                index <- index + 1
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

    # To speed up rendering, don't show every single gene.
    n_ranks <- nrow(ranking)
    sample_ranking <- ranking[seq(1, n_ranks, 5)]

    plot <- plotly::plot_ly() |>
        plotly::add_lines(
            data = sample_ranking,
            x = ~percentile,
            y = ~score,
            name = "All genes",
            hoverinfo = "skip",
            line = list(width = 4, color = base_color())
        ) |>
        plotly::layout(
            xaxis = list(
                title = "Percentile",
                tickformat = ".0%"
            ),
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

        index <- 1

        for (gene_set_name in names(gene_sets)) {
            gene_set <- gene_sets[[gene_set_name]]

            plot <- plot |> plotly::add_markers(
                data = gene_set_data[gene %chin% gene_set],
                x = ~percentile,
                y = ~score,
                name = gene_set_name,
                text = ~ glue::glue(
                    "<b>{name}</b><br>",
                    "Score: {round(score, digits = 2)}<br>",
                    "Rank: {rank}<br>",
                    "Percentile: {round(percentile * 100, digits = 2)}%"
                ),
                hoverinfo = "text",
                marker = list(
                    size = 10,
                    color = gene_set_color(index)
                )
            )

            index <- index + 1
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
                    x0 = 1 - first_not_included_rank / n_ranks,
                    x1 = 1 - last_rank / n_ranks,
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

    plot <- plotly::plot_ly() |>
        plotly::add_boxplot(
            data = ranking,
            x = "All genes",
            y = ~score,
            name = "All genes",
            showlegend = FALSE,
            line = list(color = base_color())
        ) |>
        plotly::layout(
            xaxis = list(tickvals = c("All genes", names(gene_sets))),
            yaxis = list(title = "Score")
        )

    if (length(gene_sets) > 0) {
        index <- 1

        for (gene_set_name in names(gene_sets)) {
            gene_set <- gene_sets[[gene_set_name]]

            plot <- plot |> plotly::add_boxplot(
                data = ranking[gene %chin% gene_set],
                x = gene_set_name,
                y = ~score,
                name = gene_set_name,
                showlegend = FALSE,
                line = list(color = gene_set_color(index))
            )

            index <- index + 1
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
        type = "bar",
        marker = list(color = base_color())
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

#' Plot scores in relation to chromosomal position of genes.
#'
#' @param ranking The ranking to visualize.
#' @param chromosome_name The chromosome to visualize. If this is `NULL` all,
#'   chromosomes will be included and the x-axis will show distances instead of
#'   positions.
#' @param gene_sets Named list of vectors of genes to highlight. The list names
#'   will be used as labels.
#'
#' @return A `plotly` figure.
#' @seealso ranking()
#'
#' @export
plot_scores_by_position <- function(ranking,
                                    chromosome_name = NULL,
                                    gene_sets = NULL) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
        stop("Please install \"plotly\" to use this function.")
    }

    distance_data <- if (!is.null(chromosome_name)) {
        chromosome_name_ <- chromosome_name
        geposan::distances[
            species == "hsapiens" &
                chromosome_name == chromosome_name_
        ]
    } else {
        geposan::distances[species == "hsapiens"]
    }

    data <- merge(ranking, distance_data, by = "gene")

    data <- merge(
        data,
        geposan::genes,
        by.x = "gene",
        by.y = "id"
    )

    data[, `:=`(gene_set = "All genes", color = base_color())]

    index <- 1
    for (gene_set_name in names(gene_sets)) {
        gene_set_genes <- gene_sets[[gene_set_name]]
        data[
            gene %chin% gene_set_genes,
            `:=`(
                gene_set = gene_set_name,
                color = gene_set_color(index)
            )
        ]

        index <- index + 1
    }

    # Use distances instead of positions in case all chromosomes are included.
    if (is.null(chromosome_name)) {
        data[, x := distance]
    } else {
        data[, x := start_position]
    }

    plotly::plot_ly() |>
        plotly::add_markers(
            data = data,
            x = ~x,
            y = ~score,
            name = ~gene_set,
            text = ~ glue::glue(
                "<b>{name}</b><br>",
                if (is.null(chromosome_name)) "Distance: " else "Position: ",
                "{round(x / 1000000, digits = 2)} MBp<br>",
                "Score: {round(score, digits = 2)}<br>",
                "Rank: {rank}<br>",
                "Percentile: {round(percentile * 100, digits = 2)}%"
            ),
            hoverinfo = "text",
        ) |>
        plotly::layout(
            xaxis = list(title = if (is.null(chromosome_name)) {
                "Distance (Bp)"
            } else {
                "Position (Bp)"
            }),
            yaxis = list(title = "Score")
        )
}
