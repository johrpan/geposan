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

# List of assemblies that the Ensembl Rest API advertises as chromosomes.
# Mitochondrial DNA has been manually removed. Unfortunately, species IDs from
# the Ensembl REST API don't map to dataset names in the BioMart interface.
# Because of that, we can't programatically filter chromosome names.
#
# See get_all_chromosomes()
valid_chromosome_names <- c(
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7",
    "8",
    "9",
    "10",
    "11",
    "12",
    "13",
    "14",
    "15",
    "16",
    "17",
    "18",
    "19",
    "20",
    "21",
    "22",
    "23",
    "24",
    "25",
    "26",
    "27",
    "28",
    "29",
    "Z",
    "1A",
    "4A",
    "30",
    "31",
    "32",
    "33",
    "34",
    "35",
    "36",
    "37",
    "38",
    "39",
    "40",
    "X",
    "25LG1",
    "25LG2",
    "LGE22",
    "Y",
    "41",
    "42",
    "43",
    "44",
    "45",
    "46",
    "47",
    "48",
    "49",
    "50",
    "LG34",
    "LG35",
    "2A",
    "2B",
    "LG1",
    "LG2",
    "LG3",
    "LG4",
    "LG5",
    "LG6",
    "LG7",
    "LG8",
    "LG9",
    "LG10",
    "LG11",
    "LG12",
    "LG13",
    "LG14",
    "LG15",
    "LG16",
    "LG17",
    "LG18",
    "LG19",
    "LG20",
    "LG21",
    "LG22",
    "LG23",
    "W",
    "LG24",
    "LG25",
    "LG26",
    "LG27",
    "LG28",
    "LG29",
    "LG30",
    "LG01",
    "LG02",
    "LG03",
    "LG04",
    "LG05",
    "LG06",
    "LG07",
    "LG08",
    "LG09",
    "A1",
    "A2",
    "A3",
    "B1",
    "B2",
    "B3",
    "B4",
    "C1",
    "C2",
    "D1",
    "D2",
    "D3",
    "D4",
    "E1",
    "E2",
    "E3",
    "F1",
    "F2",
    "LGE64",
    "LG7_11",
    "a",
    "b",
    "c",
    "d",
    "f",
    "g",
    "h",
    "LG28B",
    "LG30F",
    "LG36F",
    "LG37M",
    "LG42F",
    "LG44F",
    "LG45M",
    "LG48F",
    "LG49B",
    "ssa01",
    "ssa02",
    "ssa03",
    "ssa04",
    "ssa05",
    "ssa06",
    "ssa07",
    "ssa08",
    "ssa09",
    "ssa10",
    "ssa11",
    "ssa12",
    "ssa13",
    "ssa14",
    "ssa15",
    "ssa16",
    "ssa17",
    "ssa18",
    "ssa19",
    "ssa20",
    "ssa21",
    "ssa22",
    "ssa23",
    "ssa24",
    "ssa25",
    "ssa26",
    "ssa27",
    "ssa28",
    "ssa29",
    "2a",
    "2b",
    "7a",
    "7b",
    "I",
    "II",
    "III",
    "IV",
    "V",
    "VI",
    "VII",
    "VIII",
    "IX",
    "XI",
    "XII",
    "XIII",
    "XIV",
    "XV",
    "XVI",
    "LGE22C19W28_E50C23",
    "1a",
    "22a",
    "sgr01",
    "sgr02",
    "sgr03",
    "sgr04",
    "sgr05",
    "sgr06",
    "sgr07",
    "sgr08",
    "sgr09",
    "sgr10",
    "sgr11",
    "sgr12",
    "sgr13",
    "sgr14",
    "sgr15",
    "sgr16",
    "sgr17",
    "sgr18",
    "sgr19",
    "XVII",
    "XVIII",
    "XIX",
    "XX",
    "XXI",
    "XXII",
    "XXIII",
    "XXIV",
    "groupI",
    "groupII",
    "groupIII",
    "groupIV",
    "groupV",
    "groupVI",
    "groupVII",
    "groupVIII",
    "groupIX",
    "groupX",
    "groupXI",
    "groupXII",
    "groupXIII",
    "groupXIV",
    "groupXV",
    "groupXVI",
    "groupXVII",
    "groupXVIII",
    "groupXIX",
    "groupXX",
    "groupXXI",
    "2L",
    "2R",
    "3L",
    "3R",
    "MIC_1",
    "MIC_10",
    "MIC_11",
    "MIC_2",
    "MIC_3",
    "MIC_4",
    "MIC_5",
    "MIC_6",
    "MIC_7",
    "MIC_8",
    "MIC_9",
    "X1",
    "X2",
    "X3",
    "X4",
    "X5"
)

#' Get all chromosome names for an Ensembl dataset.
#'
#' The function tries to filter out valid chromosome names from the available
#' assemblies in the dataset.
get_chromosome_names <- function(dataset) {
    chromosome_names <- biomaRt::listFilterOptions(dataset, "chromosome_name")
    chromosome_names[chromosome_names %chin% valid_chromosome_names]
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
