#' Rank the results by computing a score.
#'
#' This function takes the result of [analyze()] and creates a score by
#' computing a weighted mean across the different methods' results.
#'
#' @param analysis Analysis object resulting from [analyze()].
#' @param weights Named list pairing method names with weighting factors. Only
#'   methods that are contained within this list will be included.
#'
#' @returns A ranking object. The object extends the analysis result with
#'   additional columns containing the `score` and the `rank` of each gene. It
#'   will be ordered by rank.
#'
#' @export
ranking <- function(analysis, weights) {
    if ("geposan_analysis" %chin% class(analysis)) {
        ranking <- copy(analysis$ranking)
    } else if ("geposan_results" %chin% class(analysis)) {
        ranking <- copy(analysis)
    } else {
        stop("Invalid analyis. Use geposan::analyze().")
    }

    ranking[, score := 0.0]

    for (method in names(weights)) {
        weighted <- weights[[method]] * ranking[, ..method]
        ranking[, score := score + weighted]
    }

    # Normalize scores to be between 0.0 and 1.0.
    min_score <- ranking[, min(score)]
    max_score <- ranking[, max(score)]
    score_range <- max_score - min_score
    ranking[, score := (score - min_score) / score_range]

    setorder(ranking, -score)
    ranking[, rank := .I]

    structure(
        ranking,
        class = c("geposan_ranking", "geposan_results", class(ranking))
    )
}

#' Find the best weights to rank the results.
#'
#' This function finds the optimal parameters to [ranking()] that result in the
#' reference genes ranking particulary high.
#'
#' @param analysis Results from [analyze()] or [ranking()].
#' @param methods Methods to include in the score.
#' @param reference_gene_ids IDs of the reference genes.
#' @param target The optimization target. It may be one of "mean", "median",
#'   "min" or "max" and results in the respective rank being optimized.
#'
#' @returns Named list pairing method names with their optimal weights. This
#'   can be used as an argument to [ranking()].
#'
#' @export
optimal_weights <- function(analysis, methods, reference_gene_ids,
                            target = "mean") {
    if (!any(c("geposan_analysis", "geposan_results") %chin% class(analysis))) {
        stop("Invalid analyis. Use geposan::analyze().")
    }

    # Compute the target rank of the reference genes when applying the weights.
    target_rank <- function(factors) {
        data <- ranking(analysis, as.list(factors))

        result <- data[gene %chin% reference_gene_ids, if (target == "min") {
            min(rank)
        } else if (target == "max") {
            max(rank)
        } else if (target == "mean") {
            mean(rank)
        } else {
            median(rank)
        }]

        if (result > 0) {
            result
        } else {
            Inf
        }
    }

    initial_factors <- rep(1.0, length(methods))
    names(initial_factors) <- methods

    optimal_factors <- stats::optim(initial_factors, target_rank)$par

    as.list(optimal_factors / max(abs(optimal_factors)))
}
