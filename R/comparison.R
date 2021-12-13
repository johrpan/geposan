#' Compare a set of genes with the ranking.
#'
#' @param ranking A ranking created using [ranking()].
#' @param comparison_gene_ids IDs of the genes of interest.
#'
#' @returns A comparison object with the following items:
#'   \describe{
#'     \item{`mean`}{The mean score of the genes.}
#'     \item{`min`}{The lowest score of the genes.}
#'     \item{`max`}{The highest score of the genes.}
#'     \item{`mean_rank`}{The mean rank of the genes.}
#'     \item{`first_rank`}{The first rank of the genes.}
#'     \item{`last_rank`}{The last rank of the genes.}
#'     \item{`p_value`}{p-value for the null hypothesis that the comparison
#'       genes do _not_ rank better than other genes. In other words: A low
#'       p-value means that the comparison genes rank particularly high.}
#'   }
#'
#' @export
compare <- function(ranking, comparison_gene_ids) {
    if (!inherits(ranking, "geposan_ranking")) {
        stop("Invalid ranking. Use geposan::ranking().")
    }

    comparison_ranking <- ranking[gene %chin% comparison_gene_ids]

    p_value <- stats::wilcox.test(
        x = comparison_ranking[, score],
        y = ranking[!gene %chin% comparison_gene_ids, score],
        alternative = "greater"
    )$p.value

    structure(
        list(
            mean = comparison_ranking[, mean(score)],
            min = comparison_ranking[, min(score)],
            max = comparison_ranking[, max(score)],
            mean_rank = comparison_ranking[, mean(rank)],
            first_rank = comparison_ranking[, min(rank)],
            last_rank = comparison_ranking[, max(rank)],
            p_value = p_value
        ),
        class = "geposan_comparison"
    )
}

#' S3 method to print a comparison object.
#'
#' @param x The comparison to print.
#' @param ... Other parameters.
#'
#' @seealso [compare()]
#'
#' @export
print.geposan_comparison <- function(x, ...) {
    cat("geposan comparison:")
    cat(sprintf(
        paste(
            "\n\n  Mean score: %.3f",
            "\n  Min score: %.3f",
            "\n  Max score: %.3f",
            "\n\n  Mean rank: %.1f",
            "\n  First rank: %i",
            "\n  Last rank: %i",
            "\n\n  p-value for better ranking: %.4f\n",
            sep = ""
        ),
        x$mean,
        x$min,
        x$max,
        x$mean_rank,
        x$first_rank,
        x$last_rank,
        x$p_value
    ))

    invisible(x)
}
