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
#'     \item{`test_result`}{Results of applying a Wilcoxon rank sum test.}
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

  test <- stats::wilcox.test(
    x = comparison_ranking[, score],
    y = ranking[!gene %chin% comparison_gene_ids, score],
    conf.int = TRUE
  )

  structure(
    list(
      quantiles = quantiles,
      mean_score = comparison_ranking[, mean(score)],
      mean_rank = comparison_ranking[, mean(rank)],
      mean_percentile = comparison_ranking[, mean(percentile)],
      test_result = test
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

  cat(glue::glue(
    "\n",
    "\n Mean score: {num(x$mean_score, 3)}",
    "\n Mean rank: {num(x$mean_rank, 1)}",
    "\n Mean percentile: {num(x$mean_percentile * 100, 2)}",
    "\n",
    "\n Estimated difference in medians: ",
    "{num(x$test$conf.int[1], 2)} to {num(x$test$conf.int[2], 2)}",
    "\n Confidence level: 95%",
    "\n p-value: {num(x$test$p.value, 4)}"
  ))

  invisible(x)
}
