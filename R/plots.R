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
#' @param reference_gene_ids Optionally, a set of reference genes that will be
#'   used to reorder the species.
#'
#' @export
plot_positions <- function(species_ids, gene_sets, reference_gene_ids = NULL) {
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

  # Sort species if reference genes have been provided.
  if (!is.null(reference_gene_ids)) {
    species_median_distance <- data[
      gene %chin% reference_gene_ids,
      .(median_distance = as.numeric(stats::median(distance))),
      by = species
    ]

    species <- merge(
      species,
      species_median_distance,
      by.x = "id",
      by.y = "species",
      all.x = TRUE
    )

    setorder(species, median_distance)
  }

  plot <- plotly::plot_ly() |>
    plotly::layout(
      xaxis = list(title = "Distance to telomeres [Bp]"),
      yaxis = list(
        title = "Species",
        type = "category",
        categoryorder = "array",
        categoryarray = species$id,
        tickmode = "array",
        tickvals = species$id,
        ticktext = species$name
      ),
      bargap = 0.9
    ) |>
    plotly::add_bars(
      data = species_max_distance,
      x = ~max_distance,
      y = ~species,
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
        x = ~distance,
        y = ~species,
        name = gene_set_name,
        text = ~ glue::glue(
          "<b>{name}</b><br>",
          "{round(distance / 1000000, digits = 2)} MBp"
        ),
        hoverinfo = "text",
        marker = list(
          size = 5,
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

#' Plot a scatter plot to compare two rankings.
#'
#' This function requires the package `plotly`.
#'
#' @param ranking_x The ranking to be shown on the X axis.
#' @param ranking_y The ranking to be shown on the Y axis.
#' @param name_x Title of the X axis.
#' @param name_y Title of the Y axis.
#' @param gene_sets A named list of vectors of gene IDs to highlight. The names
#'   will be used to distinguish the sets and in the legend.
#' @param use_ranks Show ranks instead of scores.
#' @param use_sample Limit genes outside of the gene sets to a random sample.
#'
#' @export
plot_rankings_correlation <- function(ranking_x,
                                      ranking_y,
                                      name_x,
                                      name_y,
                                      gene_sets = NULL,
                                      use_ranks = TRUE,
                                      use_sample = TRUE) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Please install \"plotly\" to use this function.")
  }

  data <- merge(ranking_x, ranking_y, by = "gene")
  data <- merge(data, geposan::genes, by.x = "gene", by.y = "id")

  data[, `:=`(
    x = if (use_ranks) rank.x else score.x,
    y = if (use_ranks) rank.y else score.y
  )]

  model <- stats::lm(y ~ x, data)
  model_data <- data.table(x = seq(min(data$x), max(data$x), length = 100))
  model_data[, c("y", "lower", "upper") := data.table(
    stats::predict(model, model_data, interval = "confidence")
  )]

  fig <- plotly::plot_ly()

  if (use_sample) {
    fig <- fig |> plotly::add_markers(
      data = data[!gene %chin% unlist(gene_sets)][sample(.N, 1000)],
      x = ~x,
      y = ~y,
      name = "All genes",
      marker = list(
        color = base_color(),
        size = 5
      ),
      hoverinfo = "skip"
    )
  } else {
    fig <- fig |> plotly::add_markers(
      data = data,
      x = ~x,
      y = ~y,
      name = "All genes",
      marker = list(
        color = base_color(),
        size = 5
      ),
      text = ~ glue::glue(
        "<b>{name}</b>",
        "<br>",
        "{name_x}: {round(x, digits = 2)} ",
        "({round(percentile.x * 100, digits = 2)}%)<br>",
        "{name_y}: {round(y, digits = 2)} ",
        "({round(percentile.y * 100, digits = 2)}%)"
      ),
      hoverinfo = "text"
    )
  }

  fig <- fig |>
    plotly::add_lines(
      data = model_data,
      x = ~x,
      y = ~y,
      line = list(color = base_color()),
      showlegend = FALSE,
      hoverinfo = "skip"
    ) |>
    plotly::add_ribbons(
      data = model_data,
      x = ~x,
      ymin = ~lower,
      ymax = ~upper,
      fillcolor = base_color_transparent(),
      line = list(width = 0),
      showlegend = FALSE,
      hoverinfo = "skip"
    )

  gene_set_index <- 1

  for (gene_set_name in names(gene_sets)) {
    gene_set <- gene_sets[[gene_set_name]]

    fig <- fig |>
      plotly::add_markers(
        data = data[gene %chin% gene_set],
        x = ~x,
        y = ~y,
        name = gene_set_name,
        marker = list(
          color = gene_set_color(gene_set_index),
          size = 8
        ),
        text = ~ glue::glue(
          "<b>{name}</b>",
          "<br>",
          "{name_x}: {round(x, digits = 2)} ",
          "({round(percentile.x * 100, digits = 2)}%)<br>",
          "{name_y}: {round(y, digits = 2)} ",
          "({round(percentile.y * 100, digits = 2)}%)"
        ),
        hoverinfo = "text"
      )

    gene_set_index <- gene_set_index + 1
  }

  fig <- fig |> plotly::layout(
    xaxis = list(title = name_x),
    yaxis = list(title = name_y)
  )

  if (use_ranks) {
    fig <- fig |> plotly::layout(
      xaxis = list(autorange = "reversed"),
      yaxis = list(autorange = "reversed")
    )
  }

  fig
}


#' Plot a ranking as a scatter plot of scores.
#'
#' This function requires the package `plotly`.
#'
#' @param ranking The ranking to visualize.
#' @param gene_sets A named list of gene sets (containing vectors of gene IDs)
#'   that will be highlighted in the plot. The names will be used in the legend.
#'
#' @seealso ranking()
#'
#' @export
plot_scores <- function(ranking, gene_sets = NULL) {
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
      yaxis = list(title = "Score"),
      shapes = list(
        vline(0.95),
        vline(0.75),
        vline(0.50),
        vline(0.25),
        vline(0.05)
      ),
      annotations = list(
        vlineannotation(0.95),
        vlineannotation(0.75),
        vlineannotation(0.50),
        vlineannotation(0.25),
        vlineannotation(0.05)
      )
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
      x = ~score,
      y = "All genes",
      name = "All genes",
      showlegend = FALSE,
      boxpoints = FALSE,
      line = list(color = base_color())
    ) |>
    plotly::layout(
      xaxis = list(title = "Score"),
      yaxis = list(tickvals = c("All genes", names(gene_sets)))
    )

  if (length(gene_sets) > 0) {
    index <- 1

    for (gene_set_name in names(gene_sets)) {
      gene_set <- gene_sets[[gene_set_name]]

      plot <- plot |> plotly::add_boxplot(
        data = ranking[gene %chin% gene_set],
        x = ~score,
        y = gene_set_name,
        name = gene_set_name,
        showlegend = FALSE,
        boxpoints = FALSE,
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

  data <- merge(
    geposan::distances[
      species == "9606",
      .(
        gene,
        chromosome,
        position = (start_position + end_position) / 2
      )
    ],
    ranking,
    by = "gene"
  )[, .(chromosome, position, score)]

  data <- merge(data, geposan::chromosomes, by.x = "chromosome", by.y = "id")

  global_max_position <- data[, max(position)]

  centromeres <- fread("
    name, start, end
    1, 122026460, 125184587
    2, 92188146, 94090557
    3, 90772459, 93655574
    4, 49708101, 51743951
    5, 46485901, 50059807
    6, 58553889, 59829934
    7, 58169654, 60828234
    8, 44033745, 45877265
    9, 43236168, 45518558
    10, 39686683, 41593521
    11, 51078349, 54425074
    12, 34769408, 37185252
    13, 16000001, 18051248
    14, 16000001, 18173523
    15, 17000001, 19725254
    16, 36311159, 38280682
    17, 22813680, 26885980
    18, 15460900, 20861206
    19, 24498981, 27190874
    20, 26436233, 30038348
    21, 10864561, 12915808
    22, 12954789, 15054318
    X, 58605580, 62412542
    Y, 10316945, 10544039")

  rows_chromosomes <- list(
    list("1", "2", "3", "4", "5"),
    list("6", "7", "8", "9", "10", "11", "12"),
    list("13", "14", "15", "16", "17", "18"),
    list("19", "20", "21", "22", "X", "Y")
  )

  cell_width <- floor(100 / max(lengths(rows_chromosomes)))
  row_height <- floor(100 / length(rows_chromosomes))
  rows_figs <- list()

  for (chromosomes in rows_chromosomes) {
    row_figs <- list()
    n_chromosomes <- length(chromosomes)
    first <- TRUE

    for (chromosome_name in chromosomes) {
      chromosome_data <- data[name == chromosome_name]

      # Center chromosomes horizontally
      offset <- chromosome_data[, (global_max_position - max(position)) / 2]
      chromosome_data[, position := position + offset]

      model <- stats::loess(
        score ~ position,
        chromosome_data
      )

      positions <- seq(0, global_max_position, by = global_max_position / 100)
      scores <- stats::predict(model, positions)
      centromere_position <- offset + centromeres[
        name == chromosome_name,
        (start + end) / 2
      ]

      centromere_score <- scores[!is.na(scores)][
        which.min(abs(positions[!is.na(scores)] - centromere_position))
      ]

      row_figs <- c(row_figs, list(htmltools::div(
        style = glue::glue("width: {cell_width}%; height: {row_height}%"),
        plotly::plot_ly(
          width = 150,
          height = 150
        ) |>
          plotly::add_lines(
            x = positions,
            y = scores,
            line = list(width = 4, color = base_color())
          ) |>
          plotly::add_markers(
            x = centromere_position,
            y = centromere_score,
            marker = list(
              size = 12,
              color = "white",
              line = list(width = 4, color = base_color())
            )
          ) |>
          plotly::layout(
            annotations = list(
              x = 0.5,
              y = 0.8,
              text = chromosome_name,
              xref = "paper",
              yref = "paper",
              xanchor = "center",
              showarrow = FALSE,
              font = list(size = 16)
            ),
            xaxis = list(
              range = c(0, global_max_position),
              showgrid = FALSE,
              zeroline = FALSE,
              visible = FALSE
            ),
            yaxis = list(
              range = c(0, 1),
              showgrid = FALSE,
              zeroline = FALSE,
              visible = FALSE
            ),
            showlegend = FALSE,
            margin = list(
              b = 0,
              l = 0,
              r = 0,
              t = 0
            )
          ) |>
          plotly::config(
            staticPlot = TRUE
          )
      )))

      first <- FALSE
    }

    rows_figs <- c(rows_figs, list(htmltools::div(
      style = "display: flex; justify-content: space-between",
      row_figs
    )))
  }

  htmltools::browsable(htmltools::div(
    style = "display: flex; flex-direction: column",
    rows_figs
  ))
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
    chromosome_id <- geposan::chromosomes[
      species == "9606" & name == chromosome_name,
      id
    ]

    geposan::distances[species == "9606" & chromosome == chromosome_id]
  } else {
    geposan::distances[species == "9606"]
  }

  data <- merge(ranking, distance_data, by = "gene")

  data <- merge(
    data,
    geposan::genes,
    by.x = "gene",
    by.y = "id"
  )

  # Use distances instead of positions in case all chromosomes are included.
  if (is.null(chromosome_name)) {
    data[, x := distance]
  } else {
    data[, x := (start_position + end_position) / 2]
  }

  fig <- plotly::plot_ly() |>
    plotly::add_markers(
      data = data[!gene %chin% unlist(gene_sets)],
      x = ~x,
      y = ~score,
      name = "All genes",
      text = ~ glue::glue(
        "<b>{name}</b><br>",
        if (is.null(chromosome_name)) "Distance: " else "Position: ",
        "{round(x / 1000000, digits = 2)} MBp<br>",
        "Score: {round(score, digits = 2)}<br>",
        "Rank: {rank}<br>",
        "Percentile: {round(percentile * 100, digits = 2)}%"
      ),
      marker = list(
        color = base_color(),
        size = 5
      ),
      hoverinfo = "text"
    ) |>
    plotly::layout(
      xaxis = list(title = if (is.null(chromosome_name)) {
        "Distance (Bp)"
      } else {
        "Position (Bp)"
      }),
      yaxis = list(title = "Score")
    )

  index <- 1

  for (gene_set_name in names(gene_sets)) {
    gene_set_genes <- gene_sets[[gene_set_name]]

    fig <- fig |>
      plotly::add_markers(
        data = data[gene %chin% gene_set_genes],
        x = ~x,
        y = ~score,
        name = gene_set_name,
        text = ~ glue::glue(
          "<b>{name}</b><br>",
          if (is.null(chromosome_name)) "Distance: " else "Position: ",
          "{round(x / 1000000, digits = 2)} MBp<br>",
          "Score: {round(score, digits = 2)}<br>",
          "Rank: {rank}<br>",
          "Percentile: {round(percentile * 100, digits = 2)}%"
        ),
        marker = list(
          color = gene_set_color(index),
          size = 8
        ),
        hoverinfo = "text"
      )

    index <- index + 1
  }

  fig
}

#' Helper function for creating a vertical line for plotly.
#' @noRd
vline <- function(x) {
  list(
    type = "line",
    y0 = 0,
    y1 = 1,
    yref = "paper",
    x0 = x,
    x1 = x,
    line = list(
      color = "#00000080",
      opacity = 0.5,
      dash = "dot"
    )
  )
}

#' Helper function for creating annotations for lines created using [vline()].
#' @noRd
vlineannotation <- function(x) {
  list(
    text = glue::glue("{round(x * 100)}%"),
    showarrow = FALSE,
    yref = "paper",
    x = x,
    y = 1,
    xanchor = "left",
    xshift = 4,
    align = "left",
    font = list(
      color = "#00000080"
    )
  )
}
