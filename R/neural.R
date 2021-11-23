# Find genes by training and applying a neural network.
neural <- function(preset, progress = NULL, seed = 49641) {
    species_ids <- preset$species_ids
    gene_ids <- preset$gene_ids
    reference_gene_ids <- preset$reference_gene_ids

    cached("neural", c(species_ids, gene_ids, reference_gene_ids), {
        tensorflow::set_random_seed(seed)

        gene_count <- length(gene_ids)

        progress_buffer <- 0
        progress_step <- 1 / (2 * length(reference_gene_ids) + 1)

        # Prefilter distances by species.
        distances <- geposan::distances[species %chin% species_ids]

        # Input data for the network. This contains the gene ID as an
        # identifier as well as the per-species gene distances as input
        # variables.
        data <- data.table(gene = gene_ids)

        # Buffer to keep track of the names of the input variables.
        input_vars <- NULL

        # Make a columns containing positions and distances for each
        # species.
        for (species_id in species_ids) {
            species_data <- distances[species == species_id, .(gene, distance)]

            # Only include species with at least 25% known values. As
            # positions and distances always coexist, we don't loose any
            # data here.

            species_data <- stats::na.omit(species_data)

            if (nrow(species_data) >= 0.25 * gene_count) {
                data <- merge(data, species_data, all.x = TRUE)

                # Replace missing data with mean values. The neural network
                # can't handle NAs in a meaningful way. Choosing extreme
                # values here would result in heavily biased results.
                # Therefore, the mean value is chosen as a compromise.
                # However, this will of course lessen the significance of
                # the results.

                mean_distance <- round(species_data[, mean(distance)])
                data[is.na(distance), `:=`(distance = mean_distance)]

                # Name the new column after the species.
                setnames(data, "distance", species_id)

                # Add the input variable to the buffer.
                input_vars <- c(input_vars, species_id)
            }
        }

        # Extract the reference genes.

        reference_data <- data[gene %chin% reference_gene_ids]
        reference_data[, score := 1.0]

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

        reference_samples[, score := 0.0]

        # Merge training data. The training data includes all reference
        # genes as well as an equal number of random sample genes.
        training_data <- rbindlist(list(reference_data, reference_samples))

        # Layers for the neural network.
        input_layer <- length(input_vars)
        layer1 <- input_layer
        layer2 <- 0.5 * input_layer
        layer3 <- 0.5 * layer2

        # Train the model using the specified subset of the training data and
        # apply it for predicting the genes.
        apply_network <- function(training_gene_ids, gene_ids) {
            # Create a new model for each training session, because the model
            # would keep its state across training sessions otherwise.
            model <- keras::keras_model_sequential() |>
                keras::layer_dense(
                    units = layer1,
                    activation = "relu",
                    input_shape = input_layer,
                ) |>
                keras::layer_dense(
                    units = layer2,
                    activation = "relu",
                    kernel_regularizer = keras::regularizer_l2()
                ) |>
                keras::layer_dense(
                    units = layer3,
                    activation = "relu",
                    kernel_regularizer = keras::regularizer_l2()
                ) |>
                keras::layer_dense(
                    units = 1,
                    activation = "sigmoid"
                ) |>
                keras::compile(loss = "binary_crossentropy")

            # Prepare training data by filtering it to the given genes and
            # converting it to a matrix.
            training_data <- training_data[gene %chin% training_gene_ids]
            training_matrix <- as.matrix(training_data[, ..input_vars])
            colnames(training_matrix) <- NULL
            training_matrix <- keras::normalize(training_matrix)

            fit <- keras::fit(
                model,
                x = training_matrix,
                y = training_data$score,
                epochs = 300,
                verbose = FALSE
            )

            # Convert the input data to a matrix.
            data_matrix <- as.matrix(data[gene %chin% gene_ids, ..input_vars])
            colnames(data_matrix) <- NULL
            data_matrix <- keras::normalize(data_matrix)

            data[
                gene %chin% gene_ids,
                score := stats::predict(model, data_matrix)
            ]

            if (!is.null(progress)) {
                progress_buffer <<- progress_buffer + progress_step
                progress(progress_buffer)
            }

            list(
                training_gene_ids = training_gene_ids,
                gene_ids = gene_ids,
                model = model,
                fit = fit
            )
        }

        # Apply the network to all non-training genes first.
        network <- apply_network(
            training_data$gene,
            gene_ids[!gene_ids %chin% training_data$gene]
        )

        cross_networks <- NULL

        # Apply the network to the training genes leaving out the gene itself.
        for (training_gene_id in training_data$gene) {
            cross_network <- apply_network(
                training_data[gene != training_gene_id, gene],
                training_gene_id
            )

            cross_networks <- c(cross_networks, cross_network)
        }

        structure(
            list(
                results = data[, .(gene, score)],
                network = network,
                cross_networks = cross_networks
            ),
            class = "geposan_method_results"
        )
    })
}
