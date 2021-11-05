#' Rank the results by computing a score.
#'
#' This function takes the result of [analyze()] and creates a score by
#' computing a weighted mean across the different methods' results.
#'
#' @param analysis Analysis object resulting from [analyze()].
#' @param weights Named list pairing method names with weighting factors. Only
#'   methods that are contained within this list will be included.
#'
#' @returns A ranking object. The object extends the analysis with additional
#'   columns containing the `score` and the `rank` of each gene. It will be
#'   ordered by rank.
#'
#' @export
ranking <- function(analysis, weights) {
    if (!"geposan_analysis" %chin% class(analysis)) {
        stop("Invalid analyis. Use geposan::analyze().")
    }

    ranking <- copy(analysis)
    ranking[, score := 0.0]

    for (method in names(weights)) {
        weighted <- weights[[method]] * ranking[, ..method]
        ranking[, score := score + weighted]
    }

    # Normalize scores to be between 0.0 and 1.0.
    ranking[, score := score / sum(unlist(weights))]

    setorder(ranking, -score)
    ranking[, rank := .I]

    structure(
        ranking,
        class = c("geposan_ranking", "geposan_analysis", class(ranking))
    )
}

#' S3 method for plotting a ranking.
#'
#' @param gene_sets A list of gene sets (containing vectors of gene IDs) that
#'   will be highlighted in the plot.
#' @param labels Labels for the gene sets.
#'
#' @seealso ranking()
#'
#' @export
plot.geposan_ranking <- function(ranking, gene_sets = NULL, labels = NULL) {
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

        # Include gene information which will be used for laebling
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

#' Find the best weights to rank the results.
#'
#' This function finds the optimal parameters to [ranking()] that result in the
#' reference genes ranking particulary high.
#'
#' @param analysis Results from [analyze()] or [ranking()].
#' @param methods Methods to include in the score.
#' @param reference_gene_ids IDs of the reference genes.
#' @param target The optimization target. It may be one of "mean", "min" or
#'   "max" and results in the respective rank being optimized.
#'
#' @returns Named list pairing method names with their optimal weights. This
#'   can be used as an argument to [ranking()].
#'
#' @export
optimal_weights <- function(analysis, methods, reference_gene_ids,
                            target = "mean") {
    if (!"geposan_analysis" %chin% class(analysis)) {
        stop("Invalid analyis. Use geposan::analyze().")
    }

    # Create the named list from the factors vector.
    weights <- function(factors) {
        result <- NULL

        mapply(function(method, factor) {
            result[[method]] <<- factor
        }, methods, factors)

        result
    }

    # Compute the target rank of the reference genes when applying the weights.
    target_rank <- function(factors) {
        data <- ranking(analysis, weights(factors))

        data[gene %chin% reference_gene_ids, if (target == "min") {
            min(rank)
        } else if (target == "max") {
            max(rank)
        } else {
            mean(rank)
        }]
    }

    factors <- stats::optim(rep(1.0, length(methods)), target_rank)$par
    factors[factors < 0.0] <- 0.0
    total_weight <- sum(factors)

    weights(factors / total_weight)
}
