# Title: Preoperative models validation on COVID19 dataset
# Project: Risk assessment with gene expression markers in sepsis development
# Maintainer: Albert Garcia-Lopez
# Last update: October 20, 2023


# Define output directory and libraries -------------------------------------------------------

outDir <- dirname(rstudioapi::getActiveDocumentContext()$path)
stopifnot(dir.exists(outDir) & file.exists(rstudioapi::getActiveDocumentContext()$path))
pacman::p_load(tidyverse, rio, caret, InformationValue, magrittr, doMC)
n_cores <- as.numeric(future::availableCores())


# Load custom functions -----------------------------------------------------------------------

usethis::use_zip(
  url = "https://github.com/thisisalbert/AlbertML/archive/refs/heads/main.zip",
  destdir = outDir,
  cleanup = TRUE
)

outDir %>%
  dir("\\.R$", full.names = TRUE, recursive = TRUE) %>%
  str_subset("AlbertML") %>%
  sapply(source) %>%
  invisible()


# Load preoperative machine learning RNA-seq results ------------------------------------------

comparisons <- c("Inf+ vs. Inf-", "Sepsis+ vs. SIRS-", "Sepsis+ vs. SIRS+")
attach("rnaseq_ml_results.RData")


# Define DEGs to be used ----------------------------------------------------------------------

degs_all_comparisons <-
  map(
    .x = comparisons,
    .f = function(comp) {
      all_models[[comp]] %>%
        map(~ str_subset(colnames(.x$trainingData), ".outcome", negate = TRUE)) %>%
        reduce(intersect)
    }
  ) %>%
  set_names(comparisons)


# Load COVID-19 data --------------------------------------------------------------------------

covid_data <-
  import("tpm_covid19.csv") %>%
  as_tibble() %>%
  filter(V1 %in% reduce(degs_all_comparisons, union)) %>%
  column_to_rownames("V1") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_name") %>%
  as_tibble() %>%
  mutate(type = case_when(
    str_detect(sample_name, "_06[0-9]$") ~ "Mild",
    TRUE ~ "Severe"
  ), .after = sample_name)


# Apply ML function ---------------------------------------------------------------------------

doML <- function(model, comparison, iter) {

  # Define case and control

  case <- str_subset(levels(model), "^inf_pos$|^sepsis_pos$")
  control <- str_subset(levels(model), "^inf_pos$|^sepsis_pos$", negate = TRUE)

  # Adapt data for prediction

  if (comparison %in% c("Inf+ vs. Inf-")) {
    data_to_predict <-
      covid_data %>%
      mutate(contrast = factor(case, c(case, control)), .after = type) %>%
      as.data.frame()
  }

  if (comparison %in% c("Sepsis+ vs. SIRS-", "Sepsis+ vs. SIRS+")) {
    data_to_predict <-
      covid_data %>%
      filter(type == "Severe") %>%
      mutate(contrast = factor(case, c(case, control)), .after = type) %>%
      as.data.frame()
  }

  # Prediction

  preds <-
    caret::predict.train(
      object = model,
      newdata = data_to_predict %>% select(starts_with("ENSG")),
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

  # Summarize results

  tibble(
    type = data_to_predict$type,
    reference = data_to_predict$contrast,
    prediction = opt_pred
  ) %>%
    group_by(type) %>%
    count(reference, prediction, name = "freq") %>%
    ungroup() %>%
    mutate(status = case_when(
      prediction == case & reference == case ~ "TP",
      prediction == control & reference == case ~ "FN",
      prediction == case & reference == control ~ "FP",
      prediction == control & reference == control  ~ "TN"
    )) %>%
    select(-reference, -prediction) %>%
    pivot_wider(names_from = status, values_from = freq) %>%
    rename(Subtype = type) %>%
    rowwise() %>%
    mutate(Group_Size = sum(c_across(-Subtype), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(Total_Size = sum(Group_Size), Comparison = comparison, Model = iter) %>%
    rowwise() %>%
    mutate(
      `TPR (%)` = ifelse(any(colnames(.) == "TP"), (TP/Group_Size) * 100, NA),
      `TNR (%)` = ifelse(any(colnames(.) == "TN"), (TN/Group_Size) * 100, NA),
      `FPR (%)` = ifelse(any(colnames(.) == "FP"), (FP/Group_Size) * 100, NA),
      `FNR (%)` = ifelse(any(colnames(.) == "FN"), (FN/Group_Size) * 100, NA)
    ) %>%
    ungroup() %>%
    select(
      Model, Comparison, Subtype, Group_Size, Total_Size,
      `TPR (%)`, `TNR (%)`, `FPR (%)`, `FNR (%)`,
      everything()
    )

}

registerDoMC(cores = n_cores)

results_pred_covid <-
  map(
    .x = comparisons,
    .f = function(comp) {
      imap_dfr(
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

# Mean predictions over 100 models for each comparison

results_pred_covid_mean <-
  map_dfr(
    .x = comparisons,
    .f = function(comp) {
      results_pred_covid[[comp]] %>%
        group_by(Comparison, Subtype) %>%
        summarise(across(Group_Size:last_col(), ~ mean(.x, na.rm = TRUE))) %>%
        ungroup() %>%
        mutate(across(Group_Size:last_col(), ~ na_if(.x, NaN)))
    }
  ) %>%
  select(-`FPR (%)`, -`FNR (%)`, -FP, -FN)

