# Find genes by training a neural network on reference position data.
#
# @param seed A seed to get reproducible results.
neural <- function(preset,
                   use_positions = FALSE,
                   progress = NULL,
                   seed = 448077) {
    species_ids <- preset$species_ids
    gene_ids <- preset$gene_ids
    reference_gene_ids <- preset$reference_gene_ids

    cached(
        "neural",
        c(species_ids, gene_ids, reference_gene_ids, use_positions),
        { # nolint
            set.seed(seed)
            gene_count <- length(gene_ids)

            # Prefilter distances by species.
            distances <- geposan::distances[species %chin% species_ids]

            # Input data for the network. This contains the gene ID as an
            # identifier as well as the per-species gene distances as input
            # variables.
            data <- data.table(gene = gene_ids)

            # Buffer to keep track of species included in the computation.
            # Species from `species_ids` may be excluded if they don't have
            # enough data.
            species_ids_included <- NULL

            # Make a column containing distance data for each species.
            for (species_id in species_ids) {
                species_data <- if (use_positions) {
                    setnames(distances[
                        species == species_id,
                        .(gene, position)
                    ], "position", "distance")
                } else {
                    distances[
                        species == species_id,
                        .(gene, distance)
                    ]
                }

                # Only include species with at least 25% known values.

                species_distances <- stats::na.omit(species_data)

                if (nrow(species_distances) >= 0.25 * gene_count) {
                    species_ids_included <- c(species_ids_included, species_id)
                    data <- merge(data, species_distances, all.x = TRUE)

                    # Replace missing data with mean values. The neural network
                    # can't handle NAs in a meaningful way. Choosing extreme
                    # values here would result in heavily biased results.
                    # Therefore, the mean value is chosen as a compromise.
                    # However, this will of course lessen the significance of
                    # the results.

                    mean_distance <- round(species_distances[, mean(distance)])
                    data[is.na(distance), distance := mean_distance]

                    # Name the new column after the species.
                    setnames(data, "distance", species_id)
                }
            }

            # Extract the reference genes.

            reference_data <- data[gene %chin% reference_gene_ids]
            reference_data[, neural := 1.0]

            # Take out random samples from the remaining genes. This is another
            # compromise with a negative impact on significance. Because there
            # is no information on genes with are explicitely *not* TPE-OLD
            # genes, we have to assume that a random sample of genes has a low
            # probability of including TPE-OLD genes.

            without_reference_data <- data[!gene %chin% reference_gene_ids]

            reference_samples <- without_reference_data[
                sample(
                    nrow(without_reference_data),
                    nrow(reference_data)
                )
            ]

            reference_samples[, neural := 0.0]

            # Merge training data. The training data includes all reference
            # genes as well as an equal number of random sample genes.
            training_data <- rbindlist(list(reference_data, reference_samples))

            # Construct and train the neural network.

            nn_formula <- stats::as.formula(sprintf(
                "neural~%s",
                paste(species_ids_included, collapse = "+")
            ))

            layer1 <- length(species_ids) * 0.66
            layer2 <- layer1 * 0.66
            layer3 <- layer2 * 0.66

            nn <- neuralnet::neuralnet(
                nn_formula,
                training_data,
                hidden = c(layer1, layer2, layer3),
                linear.output = FALSE
            )

            if (!is.null(progress)) {
                # We do everything in one go, so it's not possible to report
                # detailed progress information. As the method is relatively
                # quick, this should not be a problem.
                progress(0.5)
            }

            # Apply the neural network.
            data[, score := neuralnet::compute(nn, data)$net.result]

            if (!is.null(progress)) {
                # See above.
                progress(1.0)
            }

            data[, .(gene, score)]
        }
    )
}
