# This script retrieves genome data from the Ensembl database. Run
# `ensembl_species.R` first and keep its output files "species.csv" and
# "chromosomes.csv".

library(data.table)
library(DBI)
library(glue)

compara_table <- "ensembl_compara_110"

# This is the output table of this script:

genes <- data.table(
  species = integer(),
  gene = character(),
  chromosome = integer(),
  start_position = integer(),
  end_position = integer()
)

species <- fread("species.csv")
chromosomes <- fread("chromosomes.csv")

rlog::log_info("Connecting to Ensembl database server")

db <- dbConnect(
  RMariaDB::MariaDB(),
  host = "ensembldb.ensembl.org",
  port = 5306,
  user = "anonymous"
)

rlog::log_info("Retrieving human genes")

human_species_id <- 9606
human_present_row_count <- genes[species == human_species_id, .N]

if (human_present_row_count > 0) {
  rlog::log_info(glue("Skipping. Present rows: {human_present_row_count}"))
} else {
  human_table <- species[id == human_species_id, table_name]
  dbExecute(db, glue_sql("USE {`human_table`}", .con = db))

  human_chromosome_ids <- chromosomes[species == human_species_id, id]

  human_genes <- db |>
    dbGetQuery(glue_sql("
      SELECT stable_id, seq_region_id, seq_region_start, seq_region_end
        FROM gene WHERE seq_region_id IN ({human_chromosome_ids*})")) |>
    as.data.table() |>
    setnames(
      c("stable_id", "seq_region_id", "seq_region_start", "seq_region_end"),
      c("gene", "chromosome", "start_position", "end_position")
    )

  human_genes[, species := human_species_id]
  genes <- rbind(genes, human_genes)
}

dbExecute(db, glue_sql("USE {`compara_table`}", .con = db))

for (species_id in species[id != human_species_id, id]) {
  present_row_count <- genes[species == species_id, .N]
  species_name <- species[id == species_id, name]

  if (present_row_count > 0) {
    rlog::log_info(glue("Skipping species {species_id} ({species_name})"))
    rlog::log_info(glue("Present rows: {present_row_count}"))
    next
  }

  rlog::log_info(glue(
    "Retrieving genes for species {species_id} ({species_name})"
  ))

  table_name <- species[id == species_id, table_name]
  chromosome_ids <- chromosomes[species == species_id, id]

  species_genes <- db |>
    dbGetQuery(glue_sql("
      SELECT
        human.stable_id AS gene,
        species.seq_region_id AS chromosome,
        species.seq_region_start AS start_position,
        species.seq_region_end AS end_position
      FROM
        (
          SELECT
            homology_id,
            stable_id,
            seq_region_id,
            seq_region_start,
            seq_region_end
          FROM {`table_name`}.gene
            JOIN gene_member USING (stable_id)
            JOIN homology_member USING (gene_member_id)
            JOIN homology USING (homology_id)
          WHERE seq_region_id IN ({chromosome_ids*})
            AND homology.description IN (
              'ortholog_one2one',
              'ortholog_one2many',
              'ortholog_many2many'
            )
        ) AS species
        JOIN (
          SELECT
            homology_id,
            stable_id
          FROM homology_member
            JOIN gene_member USING (gene_member_id)
          WHERE taxon_id = {human_species_id}
        ) AS human ON species.homology_id = human.homology_id;
    ", .con = db)) |>
    as.data.table()

  if (nrow(species_genes) == 0) {
    rlog::log_info("No human homologs found")
  }

  species_genes[, species := species_id]
  genes <- rbind(genes, species_genes)
  fwrite(genes, "genes.csv")
}

dbDisconnect(db)