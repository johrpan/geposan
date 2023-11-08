# This is an *interactive* script for retrieving information on species from the
# Ensembl database. There are taxons with more than one entry in the database.
# For each species that has already been seen, the script asks whether to keep
# it or replace it. We recommend to choose the most generic entry in most
# cases.

library(data.table)
library(DBI)
library(glue)

# These are the output tables of this script:

species <- data.table(
  id = integer(),
  name = character(),
  scientific_name = character(),
  table_name = character()
)

chromosomes <- data.table(
  species = integer(),
  id = integer(),
  name = character(),
  length = integer()
)

rlog::log_info("Connecting to Ensembl database server")
db <- dbConnect(
  RMariaDB::MariaDB(),
  host = "ensembldb.ensembl.org",
  port = 5306,
  user = "anonymous"
)

rlog::log_info("Retrieving list of databases")
tables <- dbGetQuery(db, "SHOW DATABASES LIKE '%_core_110_%'")[, 1]

# Populates the species and chromosomes tables using data from each species'
# table within the Ensembl database. Species without a karyotype will be skipped
# without adding any information to the tables.
for (table in tables) {
  rlog::log_info(glue("Reading species information from {table}"))
  dbExecute(db, glue_sql("USE {`table`}", .con = db))

  species_id <- db |>
    dbGetQuery("
      SELECT meta_value FROM meta
        WHERE meta_key = 'species.taxonomy_id'") |>
    as.integer()

  species_name <- db |>
    dbGetQuery("
      SELECT meta_value FROM meta
        WHERE meta_key = 'species.display_name'") |>
    as.character()

  species_scientific_name <- db |>
    dbGetQuery("
      SELECT meta_value FROM meta
        WHERE meta_key = 'species.scientific_name'") |>
    as.character()

  rlog::log_info(glue(
    "Found species {species_name} ({species_scientific_name})"
  ))

  if (species[id == species_id, .N] > 0) {
    old_name <- species[id == species_id, name]
    old_scientific_name <- species[id == species_id, scientific_name]
    input <- readline(glue("\\
      Taxon already present ({old_name}, {old_scientific_name}). \\
      Replace with {species_name} ({species_scientific_name})? [y/N] "))

    if (input == "y") {
      species <- species[id != species_id]
      chromosomes <- chromosomes[species != species_id]
    } else {
      next
    }
  }

  species_chromosomes <- db |>
    dbGetQuery(glue("
      SELECT seq_region_id, seq_region.name, length
      FROM seq_region
        JOIN seq_region_attrib USING (seq_region_id)
        JOIN attrib_type USING (attrib_type_id)
      WHERE code = 'karyotype_rank'
        AND NOT EXISTS
          (SELECT * FROM seq_region_attrib AS chromosome_attrib
            JOIN attrib_type USING (attrib_type_id)
            WHERE chromosome_attrib.seq_region_id = seq_region.seq_region_id
              AND code = 'sequence_location'
              AND chromosome_attrib.value != 'nuclear_chromosome');
    ")) |>
    as.data.table() |>
    setnames("seq_region_id", "id")

  species_chromosomes[, species := species_id]

  if (nrow(species_chromosomes) == 0) {
    rlog::log_info("Skipping (no karyotype)")
    next
  }

  species <- rbind(species, data.table(
    id = species_id,
    name = species_name,
    scientific_name = species_scientific_name,
    table_name = table
  ))

  chromosomes <- rbind(chromosomes, species_chromosomes)
}

dbDisconnect(db)

fwrite(species, "species.csv")
fwrite(chromosomes, "chromosomes.csv")