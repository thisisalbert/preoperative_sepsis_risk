# Title: Preoperative models validation on external dataset
# Project: Risk assessment with gene expression markers in sepsis development
# Maintainer: Albert Garcia-Lopez
# Last update: October 13, 2023


# Define output directory and libraries -------------------------------------------------------

outDir <- dirname(rstudioapi::getActiveDocumentContext()$path)
stopifnot(dir.exists(outDir) & file.exists(rstudioapi::getActiveDocumentContext()$path))
pacman::p_load(tidyverse, magrittr, rio, caret, doMC, parallel)
n_cores <- as.numeric(future::availableCores())


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


# Load preoperative machine learning RNA-seq results ------------------------------------------

comparisons <- c("Inf+ vs. SIRS-", "Sepsis+ vs. SIRS-", "Sepsis+ vs. SIRS+")
attach("rnaseq_ml_results.RData")


# Load Barcella's data (GSE110487) ------------------------------------------------------------

# Counts

barcella_counts <-
  list.files("Barcella_GSE110487", "counts", full.names = TRUE) %>%
  import() %>%
  column_to_rownames("Geneid") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  as_tibble()

# Metadata

barcella_meta <-
  list.files("Barcella_GSE110487", "sra", full.names = TRUE, ignore.case = TRUE) %>%
  import() %>%
  as_tibble() %>%
  select(
    sample_name = `Sample Name`,
    patient = PATIENT,
    class = clinical_classification,
    timepoint = Timepoint
  ) %>%
  mutate(sample_id = paste0(patient, timepoint)) %>%
  distinct(sample_name, .keep_all = TRUE)

# Combine Barcella's metadata with counts data 

barcella_meta_counts <-
  barcella_meta %>%
  select(sample_id, class, timepoint) %>%
  inner_join(barcella_counts, by = "sample_id") %>%
  distinct(sample_id, .keep_all = TRUE)


# Machine Learning ----------------------------------------------------------------------------

doML <- function(model, comparison, iter) {

  # Parallelization

  registerDoMC(cores = n_cores)

  # Define features to be used

  features <- str_subset(colnames(model$trainingData), ".outcome", negate = TRUE)

  # Define case and control as defined in the RNA-seq preoperative model

  case <- str_subset(levels(model), "^inf_pos$|^sepsis_pos$")
  control <- str_subset(levels(model), "^inf_pos$|^sepsis_pos$", negate = TRUE)

  # Adapt data for prediction

  data_to_predict <-
    barcella_meta_counts %>%
    filter(timepoint == "T1") %>%
    mutate(contrast = factor(case, c(case, control)), .before = class) %>%
    as.data.frame()

  # Prediction

  preds <-
    caret::predict.train(
      object = model,
      newdata = data_to_predict %>% select(all_of(features)),
      type = "prob"
    )

  # Optimal cutoff

  opt_cut <-
    InformationValue::optimalCutoff(
      actuals = ifelse(model$pred$obs == case, 1, 0),
      predictedScores = model$pred[[case]],
      optimiseFor = "Both"
    )

  # Predictions based on optimal cutoff

  opt_pred <-
    factor(
      x = ifelse(preds[[case]] > opt_cut, case, control),
      levels = c(case, control)
    )

  # Confusion matrix based on optimal predictions

  cm <-
    caret::confusionMatrix(
      data = opt_pred,
      reference = data_to_predict$contrast,
      positive = case
    )

  cm$table %>%
    as_tibble() %>%
    rename_with(str_to_lower) %>%
    mutate(class = case_when(
      prediction == case & reference == case ~ "TP",
      prediction == case & reference == control ~ "FP",
      prediction == control & reference == case ~ "FN",
      prediction == control & reference == control ~ "TN"
    )) %>%
    select(class, n) %>%
    pivot_wider(names_from = class, values_from = n) %>%
    mutate(
      Accuracy = cm$overall[["Accuracy"]],
      Levels = paste(case, control, sep = ","),
      Case = case,
      Control = control,
      Comparison = comparison,
      Model = iter,
      Optimal_Cutoff = opt_cut,
      Sample_Size = sum(c_across(c(TP, FN, FP, TN))),
      `TPR (%)` = (TP/Sample_Size) * 100,
      `FNR (%)` = (FN/Sample_Size) * 100,
      DEGs_n = length(features),
      DEGs = paste(features, collapse = ",")
    ) %>%
    select(
      Model, Comparison, Levels, Case, Control,
      Optimal_Cutoff, Accuracy, Sample_Size, `TPR (%)`, `FNR (%)`, TP, FN,
      DEGs_n, DEGs
    )
}

results_pred_barcella <-
  map(
    .x = comparisons,
    .f = function(comp) {
      imap(
        .x = all_models[[comp]],
        .f = function(model, .y) {
          doML(model = model, comparison = comp, iter = .y)
        },
        .progress = list(clear = FALSE, name = comp)
      )
    }
  ) %>%
  set_names(comparisons)


# Analyze results -----------------------------------------------------------------------------

# Mean performance metrics over 100 models for each comparison

results_pred_barcella_mean <-
  map_dfr(
    .x = comparisons,
    .f = function(comp) {
      results_pred_barcella[[comp]] %>%
        bind_rows() %>%
        mutate(across(Optimal_Cutoff:FN, ~ mean(.x, na.rm = TRUE))) %>%
        select(-Model)
    }
  )

# Reference models

reference_models_ids <-
  map(
    .x = comparisons,
    .f = function(comp) {
      map_lgl(
        .x = all_models[[comp]],
        .f = ~ identical(.x, reference_models[[comp]])
      ) %>%
        which()
    }
  ) %>%
  set_names(comparisons)

results_pred_barcella_reference_models <-
  map_dfr(
    .x = comparisons,
    .f = function(comp) {
      results_pred_barcella[[comp]] %>% chuck(reference_models_ids[[comp]])
    }
  )

