#' Find genes by training and applying a neural network.
#'
#' @param id Unique ID for the method and its results.
#' @param name Human readable name for the method.
#' @param description Method description.
#' @param seed The seed will be used to make the results reproducible.
#' @param n_models This number specifies how many sets of training data should
#'   be created. For each set, there will be a model trained on the remaining
#'   training data and validated using this set. For non-training genes, the
#'   final score will be the mean of the result of applying the different
#'   models. There should be at least two training sets. The analysis will only
#'   work, if there is at least one reference gene per training set.
#' @param control_ratio The proportion of random control genes that is included
#'   in the training data sets in addition to the reference genes. This should
#'   be a numeric value between 0.0 and 1.0.
#'
#' @return An object of class `geposan_method`.
#'
#' @export
neural <- function(id = "neural",
                   name = "Neural",
                   description = "Assessment by neural network",
                   seed = 180199,
                   n_models = 5,
                   control_ratio = 0.5) {
  method(
    id = id,
    name = name,
    description = description,
    function(preset, progress) {
      species_ids <- preset$species_ids
      gene_ids <- preset$gene_ids
      reference_gene_ids <- preset$reference_gene_ids

      cached(
        "neural",
        c(
          species_ids,
          gene_ids,
          reference_gene_ids,
          seed,
          n_models,
          control_ratio
        ),
        { # nolint
          reference_count <- length(reference_gene_ids)
          stopifnot(n_models %in% 2:reference_count)

          control_count <- ceiling(reference_count * control_ratio /
            (1 - control_ratio))

          # Make results reproducible.
          tensorflow::set_random_seed(seed)

          # Step 1: Prepare input data.
          # ---------------------------

          # Prefilter distances by species and gene.
          distances <- geposan::distances[species %chin% species_ids &
            gene %chin% gene_ids]

          # Reshape data to put species into columns.
          data <- dcast(
            distances,
            gene ~ species,
            value.var = "distance"
          )

          # Replace values that are still missing with mean values for the
          # species in question.
          data[, (species_ids) := lapply(species_ids, \(species) {
            species <- get(species)
            species[is.na(species)] <- mean(species, na.rm = TRUE)
            species
          })]

          progress(0.1)

          # Step 2: Prepare training data.
          # ------------------------------

          # Take out the reference data.
          reference_data <- data[gene %chin% reference_gene_ids]
          reference_data[, score := 1.0]

          # Draw control data from the remaining genes.
          control_data <- data[!gene %chin% reference_gene_ids][
            sample(.N, control_count)
          ]
          control_data[, score := 0.0]

          # Randomly distribute the indices of the reference and control genes
          # across one bucket per model.

          reference_sets <- split(
            sample(reference_count),
            seq_len(reference_count) %% n_models
          )

          control_sets <- split(
            sample(control_count),
            seq_len(control_count) %% n_models
          )

          # Prepare the data for each model. Each model will have one pair of
          # reference and control gene sets left out for validation. The
          # training data consists of all the remaining sets.
          networks <- lapply(seq_len(n_models), \(index) {
            training_data <- rbindlist(list(
              reference_data[!reference_sets[[index]]],
              control_data[!control_sets[[index]]]
            ))

            validation_data <- rbindlist(list(
              reference_data[reference_sets[[index]]],
              control_data[control_sets[[index]]]
            ))

            list(
              training_data = training_data,
              validation_data = validation_data
            )
          })

          # Step 3: Create, train and apply neural network.
          # -----------------------------------------------

          data_matrix <- prepare_data(data, species_ids)
          output_vars <- NULL

          for (i in seq_along(networks)) {
            network <- networks[[i]]

            # Create a new model for each training session, because
            # the model would keep its state across training
            # sessions otherwise.
            model <- create_model(length(species_ids))

            # Train the model.
            fit <- train_model(
              model,
              network$training_data,
              network$validation_data,
              species_ids
            )

            # Apply the model.
            data[, new_score := stats::predict(model, data_matrix)]

            # Remove the values of the training data itself.
            data[gene %chin% network$training_data$gene, new_score := NA]

            output_var <- sprintf("score%i", i)
            setnames(data, "new_score", output_var)
            output_vars <- c(output_vars, output_var)

            # Store the details.
            networks[[i]]$model <- keras::serialize_model(model)
            networks[[i]]$fit <- fit

            progress(0.1 + i * (0.9 / n_models))
          }

          # Compute the final score as the mean score.
          data[,
            score := mean(as.numeric(.SD), na.rm = TRUE),
            .SDcols = output_vars,
            by = gene
          ]

          progress(1.0)

          result(
            method = "neural",
            scores = data[, .(gene, score)],
            details = list(
              seed = seed,
              n_models = n_models,
              all_results = data[, !..species_ids],
              networks = networks
            )
          )
        }
      )
    }
  )
}

#' Create a `keras` model based on the number of input variables.
#'
#' @param n_input_vars Number of input variables (i.e. species).
#' @return A `keras` model.
#'
#' @noRd
create_model <- function(n_input_vars) {
  # Layers for the neural network.
  layer1 <- n_input_vars
  layer2 <- 0.5 * layer1
  layer3 <- 0.5 * layer2

  keras::keras_model_sequential() |>
    keras::layer_dense(
      units = layer1,
      activation = "relu",
      input_shape = n_input_vars,
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
    keras::compile(
      loss = keras::loss_mean_absolute_error(),
      optimizer = keras::optimizer_adam()
    )
}

#' Train a model on a specific training dataset.
#'
#' @param model The model created using [create_model()]. The model will be
#'   changed reflecting the state after training.
#' @param training_data Data to fit the model to.
#' @param validation_data Additional data to assess the model performance.
#' @param input_vars Character vector of input variables that should be
#'   included.
#'
#' @return The `keras` fit object describing the training process.
#' @noRd
train_model <- function(model, training_data, validation_data, input_vars) {
  training_matrix <- prepare_data(training_data, input_vars)
  validation_matrix <- prepare_data(validation_data, input_vars)

  keras::fit(
    model,
    x = training_matrix,
    y = training_data$score,
    validation_data = list(
      x_val = validation_matrix,
      y_val = validation_data$score
    ),
    epochs = 500,
    verbose = FALSE
  )
}

#' Convert data to a matrix and normalize it.
#'
#' @param data Input data.
#' @param input_vars Character vector of input variables that should be
#'   included.
#'
#' @return A data matrix that can be used within the models.
#' @noRd
prepare_data <- function(data, input_vars) {
  data_matrix <- as.matrix(data[, ..input_vars])
  colnames(data_matrix) <- NULL
  keras::normalize(data_matrix)
}
