# This script does post processing on the data from Ensembl and imports it into
# the R package. Run this script after `ensembl_species.R` and
# `ensembl_species.R`.

library(data.table)

species <- fread("species.csv")
chromosomes <- fread("chromosomes.csv")
genes <- fread("genes.csv")

species_metadata <- chromosomes[,
  .(
    n_chromosomes = .N,
    median_chromosome_length = as.double(stats::median(length))
  ),
  by = species
]

species <- merge(
  species,
  species_metadata,
  by.x = "id",
  by.y = "species",
  sort = FALSE
)

# Remove duplicated genes within species.
genes <- genes[!duplicated(genes, by = c("species", "gene"))]

genes_chromosomes <- merge(
  genes,
  chromosomes,
  by.x = c("species", "chromosome"),
  by.y = c("species", "id"),
  sort = FALSE
)

genes_chromosomes[, distance := ifelse(
  start_position < length - end_position,
  start_position,
  length - end_position
)]

distances <- genes_chromosomes[, .(
  species,
  gene,
  chromosome,
  start_position,
  end_position,
  distance
)]

# This table will hold information on human genes.
genes <- genes_chromosomes[
  species == 9606,
  .(
    id = gene,
    chromosome = name
  )
]

genes[, name := gprofiler2::gconvert(
  id,
  target = "HGNC",
  mthreshold = 1,
  filter_na = FALSE
)$target]

# Previous versions of geposan used different species IDs. For backwards
# compatibility, convert integer IDs to character.

species[, id := as.character(id)]
distances[, species := as.character(species)]

usethis::use_data(species, overwrite = TRUE)
usethis::use_data(genes, overwrite = TRUE)
usethis::use_data(distances, overwrite = TRUE)
