library(data.table)

rlog::log_info("Connecting to Ensembl API")

# Object to access the Ensembl API. We use the US east mirror to circumvent
# current issues with the main server being temporarily unreliable.
ensembl <- biomaRt::useEnsembl("ensembl", host = "useast.ensembl.org")

# Retrieve species information.

rlog::log_info("Retrieving species information")
ensembl_datasets <- data.table(biomaRt::listDatasets(ensembl))

# Filter out species ID and name from the result.
species <- ensembl_datasets[, .(
    id = stringr::str_match(dataset, "(.*)_gene_ensembl")[, 2],
    name = stringr::str_match(description, "(.*) genes \\(.*\\)")[, 2]
)]

#' Get all chromosome names for an Ensembl dataset.
#'
#' The following chromosome naming schemes will be recognized and have been
#' sourced from Ensembl by manually screening chromosome-level assemblies.
#'
#'  - a decimal number (most species' autosomes)
#'  - X, Y, W or Z (gonosomes)
#'  - LG followed by a decimal number (some fishes)
#'  - ssa/sgr followed by a number (Atlantic salmon/Turquoise killifish)
#'
#' The function tries to filter out those chromosome names from the available
#' assemblies in the dataset.
get_chromosome_names <- function(dataset) {
    chromosome_names <- biomaRt::listFilterOptions(dataset, "chromosome_name")
    chromosome_names[stringr::str_which(
        chromosome_names,
        "^(LG|sgr|ssa)?[0-9]+|[XYWZ]$"
    )]
}

# Retrieve information on human genes. This will only include genes on
# assembled chromosomes. Chromosomes are filtered using get_chromosome_names().

rlog::log_info("Retrieving information on human genes")
dataset <- biomaRt::useDataset("hsapiens_gene_ensembl", mart = ensembl)

human_data <- data.table(biomaRt::getBM(
    attributes = c(
        "ensembl_gene_id",
        "hgnc_symbol",
        "chromosome_name",
        "start_position",
        "end_position"
    ),
    filters = "chromosome_name",
    values = get_chromosome_names(dataset),
    mart = dataset
))

# Remove duplicated gene IDs (at the time of writing, there are a handful).
human_data <- unique(human_data, by = "ensembl_gene_id")

# Only keep relevant information on genes.
genes <- human_data[, .(
    id = ensembl_gene_id,
    name = hgnc_symbol,
    chromosome = chromosome_name
)]

# Retrieve gene distance data across species.

rlog::log_info("Retrieving distance data")

# Handle the human first, as we already retrieved the data and don't need to
# filter based on orthologies.

human_data[, chromosome_length := max(end_position), by = chromosome_name]

distances <- human_data[, .(
    species = "hsapiens",
    gene = ensembl_gene_id,
    position = start_position,
    distance = pmin(
        start_position,
        chromosome_length - end_position
    )
)]

# Iterate through all other species and retrieve their distance data.
for (species_id in species[!id == "hsapiens", id]) {
    rlog::log_info(sprintf("Loading species \"%s\"", species_id))

    dataset <- biomaRt::useDataset(
        sprintf("%s_gene_ensembl", species_id),
        mart = ensembl
    )

    # Besides the attributes that are always present, we need to check for
    # human orthologs. Some species don't have that information and will be
    # skipped.
    if (!"hsapiens_homolog_ensembl_gene" %chin%
        biomaRt::listAttributes(dataset, what = "name")) {
        rlog::log_info("No data on human orthologs")
        species <- species[id != species_id]

        next
    }

    chromosome_names <- get_chromosome_names(dataset)

    # Skip the species, if there are no assembled chromosomes.
    if (length(chromosome_names) <= 0) {
        rlog::log_info("No matching chromosome assemblies")
        species <- species[id != species_id]

        next
    }

    # Retrieve information on all genes of the current species, that have
    # human orthologs. This is called "homolog" in the Ensembl schema.
    species_distances <- data.table(biomaRt::getBM(
        attributes = c(
            "hsapiens_homolog_ensembl_gene",
            "chromosome_name",
            "start_position",
            "end_position"
        ),
        filters = c("with_hsapiens_homolog", "chromosome_name"),
        values = list(TRUE, chromosome_names),
        mart = dataset
    ))

    # Only include human genes that we have information on.
    species_distances <- species_distances[
        hsapiens_homolog_ensembl_gene %chin% genes$id
    ]

    # Only include one ortholog per human gene.
    species_distances <- unique(
        species_distances,
        by = "hsapiens_homolog_ensembl_gene"
    )

    # Precompute the genes' distance to the nearest telomere.

    species_distances[,
        chromosome_length := max(end_position),
        by = chromosome_name
    ]

    species_distances <- species_distances[, .(
        species = species_id,
        gene = hsapiens_homolog_ensembl_gene,
        position = start_position,
        distance = pmin(
            start_position,
            chromosome_length - end_position
        )
    )]

    distances <- rbindlist(list(distances, species_distances))
}

# Save data in the appropriate place.

usethis::use_data(species, overwrite = TRUE)
usethis::use_data(genes, overwrite = TRUE)
usethis::use_data(distances, overwrite = TRUE)
