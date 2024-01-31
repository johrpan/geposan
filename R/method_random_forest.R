#' Predict scores using a random forest.
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
#'   work, if there is at least one reference gene per training set. By default,
#'   one model per reference gene will be used.
#' @param control_ratio The proportion of random control genes that is included
#'   in the training data sets in addition to the reference genes. This should
#'   be a numeric value between 0.0 and 1.0.
#'
#' @return An object of class `geposan_method`.
#'
#' @export
random_forest <- function(id = "rforest",
                   name = "Random forest",
                   description = "Assessment by random forest",
                   seed = 180199,
                   n_models = NULL,
                   control_ratio = 0.75) {
  method(
    id = id,
    name = name,
    description = description,
    help = paste0(
      "Assessment of a random forest model trained on reference gene ",
      "distances. It may derive patterns in positional data that have not ",
      "been covered by other methods."
    ),
    function(preset, progress) {
      species_ids <- preset$species_ids
      gene_ids <- preset$gene_ids
      reference_gene_ids <- preset$reference_gene_ids

      reference_count <- length(reference_gene_ids)
      if (is.null(n_models)) {
        n_models = reference_count
      }

      cached(
        id,
        c(
          species_ids,
          gene_ids,
          reference_gene_ids,
          seed,
          n_models,
          control_ratio
        ),
        { # nolint
          stopifnot(n_models %in% 2:reference_count)

          control_count <- ceiling(reference_count * control_ratio /
            (1 - control_ratio))

          # Make results reproducible.
          set.seed(seed)

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
          models <- lapply(seq_len(n_models), \(index) {
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

          # Step 3: Create, train and apply the models.
          # -------------------------------------------

          output_vars <- NULL

          for (i in seq_along(models)) {
            model <- models[[i]]
            forest <- ranger::ranger(
              x = model$training_data[, ..species_ids],
              y = model$training_data$score
            )

            # TODO: Make use of validation data.

            # Apply the model.
            data[, new_score := stats::predict(forest, data)$predictions]

            # Remove the values of the training data itself.
            data[gene %chin% model$training_data$gene, new_score := NA]

            output_var <- sprintf("score%i", i)
            setnames(data, "new_score", output_var)
            output_vars <- c(output_vars, output_var)

            # Store the details.
            models[[i]]$forest <- forest

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
            method = id,
            scores = data[, .(gene, score)],
            details = list(
              seed = seed,
              n_models = n_models,
              all_results = data[, !..species_ids],
              models = models
            )
          )
        }
      )
    }
  )
}
