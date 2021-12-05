library(data.table)
library(httr)

ensembl_api_url <- "https://rest.ensembl.org"

#' Perform a request to the Ensembl REST API.
ensembl_request <- function(api_path) {
    content(stop_for_status(GET(
        paste0(ensembl_api_url, api_path),
        content_type_json()
    )))
}

#' Get IDs of all available vertebrates.
get_species_ids <- function() {
    species <- ensembl_request("/info/species")$species
    sapply(species, function(species) species$name)
}

#' Get all chromosomes names for a species.
get_species_chromosomes <- function(species_id) {
    chromosomes <- unlist(ensembl_request(
        paste0("/info/assembly/", species_id)
    )$karyotype)
}

#' Get a vector of all available unqiue chromosome names.
#'
#' There are multiple names for mitochondrial DNA which have to be removed
#' manually, unfortunately.
get_all_chromosomes <- function() {
    chromosomes <- sapply(get_species_ids(), get_species_chromosomes)
    unique(unlist(chromosomes))
}
