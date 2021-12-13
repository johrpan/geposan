#' Information on included species from the Ensembl database.
#'
#' @format A [data.table] with the following columns:
#' \describe{
#'   \item{id}{Unique species ID}
#'   \item{name}{Human readable species name}
#' }
"species"

#' Information on human genes within the Ensembl database.
#'
#' This includes only genes on the primary suggested assembly of the human
#' nuclear DNA.
#'
#' @format A [data.table] with the following columns:
#' \describe{
#'   \item{id}{Ensembl gene ID}
#'   \item{name}{The gene's HGNC name}
#'   \item{chrosome}{The human chromosome the gene is located on}
#' }
"genes"

#' Information on gene positions across species.
#'
#' This dataset contains each known value for a gene's distance to the telomeres
#' per species. The data is sourced from Ensembl.
#'
#' @format A [data.table] with the following columns:
#' \describe{
#'   \item{species}{Species ID}
#'   \item{gene}{Gene ID}
#'   \item{distance}{Distance to nearest telomere}
#' }
"distances"
