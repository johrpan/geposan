# Find genes by training and applying a neural network.
#
# @param seed The seed will be used to make the results reproducible.
# @param n_models This number specifies how many sets of training data should
#   be created. For each set, there will be a model trained on the remaining
#   training data and validated using this set. For non-training genes, the
#   final score will be the mean of the result of applying the different
#   models.
neural <- function(preset, progress = NULL, seed = 49641, n_models = 5) {
    species_ids <- preset$species_ids
    gene_ids <- preset$gene_ids
    reference_gene_ids <- preset$reference_gene_ids

    cached(
        "neural",
        c(species_ids, gene_ids, reference_gene_ids, seed, n_models),
        { # nolint
            reference_count <- length(reference_gene_ids)
            if (!n_models %in% 2:reference_count) {
                stop(paste0(
                    "n_models has to be between 2 and the number of reference ",
                    "genes."
                ))
            }

            # Make results reproducible.
            tensorflow::set_random_seed(seed)

            # Step 1: Prepare input data.
            # ---------------------------

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
                species_data <- distances[
                    species == species_id,
                    .(gene, distance)
                ]

                # Only include species with at least 25% known values. As
                # positions and distances always coexist, we don't loose any
                # data here.

                species_data <- stats::na.omit(species_data)

                if (nrow(species_data) >= 0.25 * length(gene_ids)) {
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

            if (!is.null(progress)) {
                progress(0.1)
            }

            # Step 2: Prepare training data.
            # ------------------------------

            # Take out the reference data.

            reference_data <- data[gene %chin% reference_gene_ids]
            reference_data[, score := 1.0]

            # Take out random samples from the remaining genes. This is another
            # compromise with a negative impact on significance. Because there
            # is no information on genes with are explicitely *not* TPE-OLD
            # genes, we have to assume that a random sample of genes has a low
            # probability of including TPE-OLD genes.

            without_reference_data <- data[!gene %chin% reference_gene_ids]

            control_data <- without_reference_data[
                sample(
                    nrow(without_reference_data),
                    reference_count
                )
            ]

            control_data[, score := 0.0]

            # Split the training data into random sets to have validation data
            # for each model.

            # Scramble the source tables.
            reference_data <- reference_data[sample(reference_count)]
            control_data <- control_data[sample(reference_count)]

            networks <- list()

            indices <- seq_len(reference_count)
            indices_split <- split(indices, indices %% n_models)

            for (i in seq_len(n_models)) {
                training_data <- rbindlist(list(
                    reference_data[!indices_split[[i]]],
                    control_data[!indices_split[[i]]]
                ))

                validation_data <- rbindlist(list(
                    reference_data[indices_split[[i]]],
                    control_data[indices_split[[i]]]
                ))

                networks[[i]] <- list(
                    training_data = training_data,
                    validation_data = validation_data
                )
            }

            # Step 3: Create, train and apply neural network.
            # -----------------------------------------------

            # Layers for the neural network.
            input_layer <- length(input_vars)
            layer1 <- input_layer
            layer2 <- 0.5 * input_layer
            layer3 <- 0.5 * layer2

            # Convert data to matrix and normalize it.
            to_matrix <- function(data) {
                data_matrix <- as.matrix(data[, ..input_vars])
                colnames(data_matrix) <- NULL
                keras::normalize(data_matrix)
            }

            data_matrix <- to_matrix(data)
            output_vars <- NULL

            for (i in seq_along(networks)) {
                # Create a new model for each training session, because the
                # model would keep its state across training sessions otherwise.
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

                # Train the model.

                network <- networks[[i]]

                training_data <- network$training_data
                training_matrix <- to_matrix(training_data)
                validation_data <- network$validation_data
                validation_matrix <- to_matrix(validation_data)

                fit <- keras::fit(
                    model,
                    x = training_matrix,
                    y = training_data$score,
                    validation_data = list(
                        x_val = validation_matrix,
                        y_val = validation_data$score
                    ),
                    epochs = 300,
                    verbose = FALSE
                )

                # Apply the model.

                data[, new_score := stats::predict(model, data_matrix)]

                # Remove the values of the training data itself.
                data[gene %chin% training_data$gene, new_score := NA]

                output_var <- sprintf("score%i", i)
                setnames(data, "new_score", output_var)
                output_vars <- c(output_vars, output_var)


                # Store the details.

                networks[[i]]$model <- model
                networks[[i]]$fit <- fit

                if (!is.null(progress)) {
                    progress(0.1 + i * (0.9 / n_models))
                }
            }

            # Compute the final score as the mean score.
            data[,
                score := mean(as.numeric(.SD), na.rm = TRUE),
                .SDcols = output_vars,
                by = gene
            ]

            if (!is.null(progress)) {
                progress(1.0)
            }

            structure(
                list(
                    results = data[, .(gene, score)],
                    seed = seed,
                    n_models = n_models,
                    all_results = data[, !..input_vars],
                    networks = networks
                ),
                class = "geposan_method_results"
            )
        }
    )
}
