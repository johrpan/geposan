#' Rank the results by computing a score.
#'
#' This function takes the result from [analyze()] and creates a score by
#' computing a weighted mean across the different methods' results.
#'
#' @param results Results from [analyze()].
#' @param weights Named list pairing method names with weighting factors.
#'
#' @result The input data with an additional column containing the score and
#'   another column containing the rank.
#'
#' @export
ranking <- function(results, weights) {
    results <- copy(results)
    results[, score := 0.0]

    for (method in names(weights)) {
        weighted <- weights[[method]] * results[, ..method]
        results[, score := score + weighted]
    }

    # Normalize scores to be between 0.0 and 1.0.
    results[, score := score / sum(unlist(weights))]

    setorder(results, -score)
    results[, rank := .I]
}

#' Find the best weights to rank the results.
#'
#' This function finds the optimal parameters to [ranking()] that result in the
#' reference genes ranking particulary high.
#'
#' @param results Results from [analyze()] or [ranking()].
#' @param methods Methods to include in the score.
#' @param reference_gene_ids IDs of the reference genes.
#' @param target The optimization target. It may be one of "mean", "min" or
#'   "max" and results in the respective rank being optimized.
#'
#' @returns Named list pairing method names with their optimal weights.
#'
#' @export
optimize_weights <- function(results, methods, reference_gene_ids,
                             target = "mean") {
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
        data <- ranking(results, weights(factors))

        data[gene %chin% reference_gene_ids, if (target == "min") {
            min(rank)
        } else if (target == "max") {
            max(rank)
        } else {
            mean(rank)
        }]
    }

    factors <- stats::optim(rep(1.0, length(methods)), target_rank)$par
    total_weight <- sum(factors)

    weights(factors / total_weight)
}
