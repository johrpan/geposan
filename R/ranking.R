#' Rank the results by computing a score.
#'
#' This function takes the result of [analyze()] and creates a score by
#' computing a weighted mean across the different methods' results.
#'
#' @param analysis Analysis object resulting from [analyze()].
#' @param weights Named list pairing method names with weighting factors. Only
#'   methods that are contained within this list will be included.
#' @param min_n_species Minimum number of required species per gene. Genes that
#'   have fewer species will not be included in the ranking.
#'
#' @returns A ranking object. The object extends the analysis result with
#'   additional columns containing the `score` and the `rank` of each gene. It
#'   will be ordered by rank.
#'
#' @export
ranking <- function(analysis, weights, min_n_species = 10) {
    if (!"geposan_analysis" %chin% class(analysis)) {
        stop("Invalid analyis. Use geposan::analyze().")
    }

    # Count included species from the preset per gene.
    genes_n_species <- geposan::distances[
        species %chin% analysis$preset$species_ids,
        .(n_species = .N),
        by = "gene"
    ]

    setkey(genes_n_species, gene)

    # Exclude genes with too few species.
    ranking <- analysis$results[
        genes_n_species[gene, n_species] >= min_n_species
    ]

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
        class = c("geposan_ranking", "geposan_analysis", class(ranking))
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
#' @param target The optimization target. It may be one of "mean", "min" or
#'   "max" and results in the respective rank being optimized.
#' @param min_n_species Minimum number of required species per gene. Genes that
#'   have fewer species will not be included in the rankings used to find the
#'   optimal weights.
#'
#' @returns Named list pairing method names with their optimal weights. This
#'   can be used as an argument to [ranking()].
#'
#' @export
optimal_weights <- function(analysis, methods, reference_gene_ids,
                            target = "mean", min_n_species = 10) {
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
        data <- ranking(
            analysis,
            weights(factors),
            min_n_species = min_n_species
        )

        result <- data[gene %chin% reference_gene_ids, if (target == "min") {
            min(rank)
        } else if (target == "max") {
            max(rank)
        } else {
            mean(rank)
        }]

        if (result > 0) {
            result
        } else {
            Inf
        }
    }

    factors <- stats::optim(
        rep(0.0, length(methods)),
        target_rank
    )$par

    weights(factors / max(abs(factors)))
}
