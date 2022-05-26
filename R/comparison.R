#' Compare a set of genes with the ranking.
#'
#' @param ranking A ranking created using [ranking()].
#' @param comparison_gene_ids IDs of the genes of interest.
#'
#' @returns A comparison object with the following items:
#'   \describe{
#'     \item{`quantiles`}{A `data.table` containing quantile values for the
#'       score, rank and percentile of the comparison genes.
#'     }
#'     \item{`mean_score`}{The mean score of the genes.}
#'     \item{`mean_rank`}{The mean rank of the genes.}
#'     \item{`mean_percentile`}{The mean percentile of the genes.}
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

  quantiles <- data.table(
    quantile = c("0%", "25%", "50%", "75%", "100%"),
    score = stats::quantile(comparison_ranking[, score]),
    rank = stats::quantile(
      comparison_ranking[, rank],
      probs = seq(1, 0, -0.25)
    ),
    percentile = stats::quantile(comparison_ranking[, percentile])
  )

  p_value <- stats::wilcox.test(
    x = comparison_ranking[, score],
    y = ranking[!gene %chin% comparison_gene_ids, score],
    alternative = "greater"
  )$p.value

  structure(
    list(
      quantiles = quantiles,
      mean_score = comparison_ranking[, mean(score)],
      mean_rank = comparison_ranking[, mean(rank)],
      mean_percentile = comparison_ranking[, mean(percentile)],
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
  cat("geposan comparison:\n\n")

  quantiles_formatted <- x$quantiles[, .(
    "Quantile" = quantile,
    "Score" = round(score, 3),
    "Rank" = rank,
    "Percentile" = paste0(
      format(round(percentile * 100, 1), nsmall = 1),
      "%"
    )
  )]

  print(quantiles_formatted, row.names = FALSE)

  cat(sprintf(
    paste0(
      "\n  Mean score: %.3f",
      "\n  Mean rank: %.1f",
      "\n  Mean percentile: %.1f%%",
      "\n  p-value for better scores: %.4f\n"
    ),
    x$mean_score,
    x$mean_rank,
    x$mean_percentile * 100,
    x$p_value
  ))

  invisible(x)
}
