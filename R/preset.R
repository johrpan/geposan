#' Create a new preset.
#'
#' A preset is used to specify which methods and inputs should be used for an
#' analysis. Note that the genes to process should normally include the
#' reference genes to be able to assess the results later. The genes will be
#' filtered based on how many species have data for them. Afterwards, species
#' that still have many missing genes will also be excluded. See the different
#' method functions for the available methods: [clustering()], [correlation()],
#' [neural()], [adjacency()] and [species_adjacency()].
#'
#' @param reference_gene_ids IDs of reference genes to compare to.
#' @param methods List of methods to apply.
#' @param species_ids IDs of species to include.
#' @param gene_ids IDs of genes to screen.
#' @param species_requirement The proportion of species a gene has to have
#'  orthologs in in order for the gene to qualify.
#' @param gene_requirement The proportion of genes that a species has to have
#'  in order for the species to be included in the analysis.
#'
#' @return The preset to use with [analyze()].
#'
#' @export
preset <- function(reference_gene_ids,
                   methods = all_methods(),
                   species_ids = geposan::species$id,
                   gene_ids = geposan::genes$id,
                   species_requirement = 0.25,
                   gene_requirement = 0.5) {
  # Prefilter distances.
  distances <- geposan::distances[
    species %chin% species_ids & gene %chin% gene_ids
  ]

  # Count included species per gene.
  genes_n_species <- distances[, .(n_species = .N), by = "gene"]

  # Filter out genes with less too few existing orthologs.
  gene_ids_filtered <- genes_n_species[
    n_species >= species_requirement * length(species_ids),
    gene
  ]

  # Count included genes per species.
  species_n_genes <- geposan::distances[, .(n_genes = .N), by = "species"]

  # Filter out species that have too few of the genes.
  species_ids_filtered <- species_n_genes[
    n_genes >= gene_requirement * length(gene_ids_filtered),
    species
  ]

  reference_gene_ids_excluded <- reference_gene_ids[
    !reference_gene_ids %chin% gene_ids_filtered
  ]

  if (length(reference_gene_ids_excluded > 0)) {
    warning(paste0(
      "The following reference gene IDs are excluded from the preset ",
      "because they don't have enough data: ",
      paste(reference_gene_ids_excluded, collapse = ", ")
    ))
  }

  reference_gene_ids_included <- reference_gene_ids[
    reference_gene_ids %chin% gene_ids_filtered
  ]

  if (length(reference_gene_ids_included) < 1) {
    stop(paste0(
      "There has to be at least one reference gene for the preset to be ",
      "valid. Please note that some methods may require more reference ",
      "genes."
    ))
  }

  # The included data gets sorted to be able to produce predictable hashes
  # for the object later.
  structure(
    list(
      reference_gene_ids = sort(reference_gene_ids_included),
      methods = methods,
      species_ids = sort(species_ids_filtered),
      gene_ids = sort(gene_ids_filtered)
    ),
    class = "geposan_preset"
  )
}

#' S3 method to print a preset object.
#'
#' @param x The preset to print.
#' @param ... Other parameters.
#'
#' @seealso [preset()]
#'
#' @export
print.geposan_preset <- function(x, ...) {
  cat(sprintf(
    paste0(
      "geposan preset:",
      "\n  Reference genes: %i",
      "\n  Included methods: %s",
      "\n  Number of species: %i",
      "\n  Number of genes: %i",
      "\n"
    ),
    length(x$reference_gene_ids),
    paste(sapply(x$methods, function(m) m$id), collapse = ", "),
    length(x$species_ids),
    length(x$gene_ids)
  ))

  invisible(x)
}
