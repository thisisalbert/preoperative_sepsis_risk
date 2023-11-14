# Title: Model training preoperative
# Project: Risk assessment with gene expression markers in sepsis development
# Maintainer: Albert Garcia-Lopez
# Last update: June 8, 2023


# Define output directory and libraries -------------------------------------------------------

outDir <- dirname(rstudioapi::getActiveDocumentContext()$path)
stopifnot(dir.exists(outDir) & file.exists(rstudioapi::getActiveDocumentContext()$path))
pacman::p_load(tidyverse, rio, caret, Boruta, DescTools, InformationValue, parallel, doParallel)
n_cores <- detectCores()


# Load custom functions -----------------------------------------------------------------------

usethis::use_zip(
  url = "https://github.com/thisisalbert/AlbertML/archive/refs/heads/main.zip",
  destdir = outDir,
  cleanup = TRUE
)

outDir %>%
  dir("\\.R$", full.names = TRUE, recursive = TRUE) %>%
  str_subset("AlbertML-main/R") %>%
  sapply(source) %>%
  invisible()


# Analyze RFE results -------------------------------------------------------------------------

# For each comparison and for each model, select the amount of variables that accounted for at
# least 80% AUROC and pick the smallest set out of the possible options.

sel_degs_rfe <-
  map(
    .x = comparisons,
    .f = function(comp) {
      imap_dfr(
        .x = results_rfe[[comp]],
        .f = function(res, .y) {
          sel_subset <-
            res %>%
            chuck("model_rfe") %>%
            chuck("results") %>%
            as_tibble() %>%
            filter(AUROC >= 0.8) %>%
            slice_min(order_by = Variables, n = 1, with_ties = FALSE) %>%
            select(Variables, AUROC)

          res %>%
            chuck("model_rfe") %>%
            chuck("variables") %>%
            as_tibble() %>%
            filter(Variables == sel_subset$Variables) %>%
            group_by(var) %>%
            summarise(Overall = mean(Overall, na.rm = TRUE)) %>%
            slice_max(order_by = Overall, n = sel_subset$Variables, with_ties = FALSE) %>%
            mutate(
              model = .y,
              comparison = comp,
              AUROC = sel_subset$AUROC,
              n_var = nrow(.),
              .before = everything()
            ) %>%
            rename(importance = Overall) %>%
            arrange(desc(importance))
        }
      )
    }
  ) %>%
  set_names(comparisons)

# Select the model/variables with the least amount of variables smaller than 20.

selected_degs <-
  map(
    .x = comparisons,
    .f = function(comp) {
      id_model <-
        sel_degs_rfe[[comp]] %>%
        filter(n_var <= 20 & n_var > 1) %>%
        slice_max(order_by = AUROC, n = 1, with_ties = FALSE) %>%
        pull(model)

      sel_degs_rfe[[comp]] %>%
        filter(model == id_model) %>%
        pull(var)
    }
  ) %>%
  set_names(comparisons)


# Model training ------------------------------------------------------------------------------

doFinalTraining <- function(comparison) {

  res_final <- list()
  i <- 1
  total_iter <- 100
  pb <- txtProgressBar(min = i, max = total_iter, style = 3, char = "=")

  # Parallelization

  registerDoParallel(cores = n_cores)

  # Start iteration

  message("\n[", comparison, "]\n")

  while (i <= total_iter) {

    # Set progress bar

    setTxtProgressBar(pb, i)

    # Set fixed seed

    fix_seed <- sample(seq_seeds, 1)
    seq_seeds <- seq_seeds[seq_seeds != fix_seed]

    # Load selected DEG subset for this comparison

    degs <- selected_degs[[comparison]]

    # Filter and balance sexes in metadata file

    set.seed(fix_seed)
    metadata_balanced <-
      metadata_extended %>%
      select(id, sex, contrast = all_of(comparison)) %>%
      drop_na(contrast) %>%
      mutate(sex = factor(sex, c("Female", "Male"))) %>%
      group_by(sex) %>%
      slice_sample(n = min(table(.$sex)), replace = FALSE) %>%
      ungroup() %>%
      as.data.frame()

    # Define case and control

    group_levels <- unique(metadata_balanced$contrast)
    case <- str_subset(group_levels, "^inf_pos$|^sepsis_pos$")
    control <- group_levels[group_levels != case]

    # Combine metadata with counts matrix

    metadata_balanced_with_counts <-
      tpm_counts %>%
      filter(rownames(.) %in% degs) %>%
      select(all_of(metadata_balanced$id)) %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("id") %>%
      as_tibble() %>%
      inner_join(metadata_balanced, by = "id") %>%
      select(colnames(metadata_balanced), everything()) %>%
      mutate(contrast = factor(contrast, c(case, control))) %>%
      as.data.frame()

    # Sanity checks

    stopifnot(
      all(prop.table(table(metadata_balanced_with_counts$sex)) == 0.5) &
      all(str_subset(colnames(metadata_balanced_with_counts), "^ENSG") %in% degs) &
      all(as.character(unique(metadata_balanced_with_counts$contrast)) %in% c(case, control)) &
      is.data.frame(metadata_balanced_with_counts)
    )

    # Training

    set.seed(fix_seed)
    model_final <-
      try(
        expr = suppressWarnings(
          caret::train(
            x = metadata_balanced_with_counts %>% select(all_of(degs)),
            y = metadata_balanced_with_counts %>% pull(contrast),
            method = "ranger",
            preProcess = c("center", "scale", "nzv", "zv"),
            tuneLength = 25,
            metric = "AUROC",
            maximize = TRUE,
            trControl = trainControl(
              method = "cv",
              number = 10,
              search = "random",
              sampling = "down",
              returnData = TRUE,
              returnResamp = "final",
              savePredictions = "final",
              classProbs = TRUE,
              summaryFunction = improvedSummary,
              selectionFunction = "best",
              seeds = makeCVseeds(type = "tuneLength", tunes = 25, numCV = 10),
              allowParallel = TRUE
            )
          )
        ),
        silent = TRUE
      )

    if (inherits(model_final, "try-error")) next

    # Export

    res_final[[i]] <-
      list(
        i = i,
        fix_seed = fix_seed,
        metadata_balanced_with_counts = metadata_balanced_with_counts,
        case = case,
        control = control,
        degs = degs,
        model_final = model_final
      )

    i <- i + 1

  }

  close(pb)
  stopImplicitCluster()
  gc()
  return(res_final)

}

results_final <-
  map(
    .x = comparisons,
    .f = ~ doFinalTraining(comparison = .x)
  ) %>%
  set_names(comparisons)


# Analyze results -----------------------------------------------------------------------------

results_final_metrics <-
  map(
    .x = comparisons,
    .f = function(comp) {
      imap_dfr(
        .x = results_final[[comp]],
        .f = function(sub_results_final, .y) {
          sub_results_final %>%
            chuck("model_final") %>%
            chuck("resample") %>%
            as_tibble() %>%
            summarise(across(-Resample, ~ mean(.x, na.rm = TRUE))) %>%
            mutate(
              Model = .y,
              Comparison = comp,
              Levels = paste(sub_results_final$model_final$levels, collapse = ","),
              DEGs_n = length(predictors(sub_results_final$model_final)),
              DEGs = paste(predictors(sub_results_final$model_final), collapse = ",")
            ) %>%
            relocate(Model, Comparison, Levels, .before = everything()) %>%
            relocate(DEGs_n, DEGs, .after = last_col())
        }
      )
    }
  ) %>%
  set_names(comparisons)

results_final_metrics_list <-
  map(
    .x = results_final_metrics,
    .f = function(res) {
      res_all_models <-
        res %>%
        mutate(across(is.numeric, ~ na_if(.x, NaN))) %>%
        mutate(across(AUROC:F1, ~ calcConfInt(.x, perc = 0.95), .names = "{.col} 95CI")) %>%
        mutate(PrevalencePPVandNPV(
          Sens = Sensitivity, Spec = Specificity, Min_Prevalence = 0.1, Max_Prevalence = 0.5
        )) %>%
        select(
          Model, Comparison, Levels, starts_with("AUROC"), starts_with("Accuracy"),
          starts_with("Balanced_Accuracy"), starts_with("Kappa"),
          starts_with("Sensitivity"), starts_with("Specificity"),
          starts_with("F1"), starts_with("PPV"), starts_with("NPV"),
          DEGs_n, DEGs
        )

      res_mean <-
        res_all_models %>%
        mutate(across(c(
          AUROC, Accuracy, Balanced_Accuracy, Kappa, Sensitivity, Specificity, F1,
          PPV_Prevalence_0.1, PPV_Prevalence_0.2, PPV_Prevalence_0.3, PPV_Prevalence_0.4,
          PPV_Prevalence_0.5, NPV_Prevalence_0.1, NPV_Prevalence_0.2, NPV_Prevalence_0.3,
          NPV_Prevalence_0.4, NPV_Prevalence_0.5
        ), ~ mean(.x, na.rm = TRUE))) %>%
        distinct(AUROC, .keep_all = TRUE)

      list(
        res_all_models = res_all_models,
        res_mean = res_mean
      )
    }
  )

results_final_metrics_mean <-
  map_dfr(
    .x = results_final_metrics_list,
    .f = ~ .x %>% chuck("res_mean") %>% select(-Model)
  )

# Select reference models

reference_models <-
  map_dfr(
    .x = comparisons,
    .f = function(comp) {
      mean_values = mean_metrics %>% filter(Comparison == comp)

      all_models_metrics %>%
        filter(Comparison == comp) %>%
        mutate(Distance = sqrt(
          (AUROC - mean_values$AUROC)^2 +
            (Accuracy - mean_values$Accuracy)^2 +
            (Balanced_Accuracy - mean_values$Balanced_Accuracy)^2 +
            (Kappa - mean_values$Kappa)^2 +
            (Sensitivity - mean_values$Sensitivity)^2 +
            (Specificity - mean_values$Specificity)^2 +
            (F1 - mean_values$F1)^2
        )) %>%
        slice_min(Distance, n = 1)
    }
  )

