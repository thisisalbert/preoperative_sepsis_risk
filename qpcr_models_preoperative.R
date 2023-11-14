# Title: qPCR models preoperative
# Project: Risk assessment with gene expression markers in sepsis development
# Maintainer: Albert Garcia-Lopez
# Last update: September 13, 2023


# Define output directory and libraries -------------------------------------------------------

outDir <- dirname(rstudioapi::getActiveDocumentContext()$path)
stopifnot(dir.exists(outDir) & file.exists(rstudioapi::getActiveDocumentContext()$path)
pacman::p_load(tidyverse, rio, caret, progress, doMC, ROSE, magrittr)
n_cores <- as.numeric(future::availableCores())


# Load custom functions -----------------------------------------------------------------------

usethis::use_zip(
  url = "https://github.com/thisisalbert/AlbertML/archive/refs/heads/main.zip",
  destdir = outDir,
  cleanup = TRUE
)

outDir %>%
  dir("\\.R", recursive = TRUE, full.names = TRUE) %>%
  str_subset("AlbertML") %>%
  lapply(source) %>%
  invisible()


# Load qPCR data ------------------------------------------------------------------------------

qpcr_data <- import("qpcr_data.rds")


# Load DEGs -----------------------------------------------------------------------------------

degs_rnaseq_ml <-
  import("rnaseq_reference_models.xlsx") %>%
  as_tibble() %>%
  select(Comparison, DEGs) %>%
  mutate(DEGs = str_split(DEGs, ",")) %>%
  deframe() %>%
  map(~ AnnotationDbi::mapIds(
    x = org.Hs.eg.db::org.Hs.eg.db,
    keys = .x,
    keytype = "ENSEMBL",
    column = "SYMBOL"
  ) %>% as.character())


# Generate random sequence of seeds -----------------------------------------------------------

set.seed(230913)
seq_seeds <- sample.int(n = 2^15, size = 1e4, replace = FALSE)


# Machine Learning ----------------------------------------------------------------------------

doML <- function(comparison, type, iterations = 100) {

  # Define settings

  res_ml <- list()
  i <- 1
  pb <- makePB(iterations)

  # Parallelization

  registerDoMC(cores = n_cores)

  # Loop

  message("\n[", comparison, "]\n")

  while (i < iterations + 1) {

    # Set progress bar

    if (pb$finished == FALSE) {
      pb$update(i/iterations)
    } else {
      break
    }

    # Set fixed seed

    fix_seed <- sample(seq_seeds, 1)
    seq_seeds <- seq_seeds[seq_seeds != fix_seed]

    # Select DEGs

    degs <- degs_rnaseq_ml[[comparison]]

    # Define case and control

    levels_comparison <- unique(qpcr_data[[comparison]][[type]]$contrast)
    case <- str_subset(levels_comparison, "^inf_pos$|^sepsis_pos$")
    control <- str_subset(levels_comparison, "^inf_pos$|^sepsis_pos$", negate = TRUE)

    # Define data selection

    set.seed(fix_seed)
    data_selection <-
      qpcr_data[[comparison]][[type]] %>%
      select(id, sex, contrast, all_of(degs)) %>%
      slice_sample(n = min(table(.$sex)), replace = FALSE, by = sex) %>%
      mutate(contrast = factor(contrast, c(case, control))) %>%
      as.data.frame()

    # Sanity checks

    stopifnot(
      anyNA(data_selection) == FALSE &
      as.logical(anyDuplicated(data_selection)) == FALSE &
      all(degs %in% colnames(data_selection)) &
      all(prop.table(table(data_selection$sex)) == 0.5) &
      all(levels(data_selection$contrast) == c(case, control)) &
      is.data.frame(data_selection)
    )

    # Training

    if (all(prop.table(table(data_selection$contrast)) == 0.5)) {
      subsampling <- NULL
    } else {
      subsampling <- "rose"
    }

    set.seed(fix_seed)
    model_train <-
      try(
        expr = suppressMessages(suppressWarnings(
          caret::train(
            x = data_selection %>% select(all_of(degs)),
            y = data_selection %>% pull(contrast),
            method = "ranger",
            preProcess = c("center", "scale", "spatialSign"),
            tuneLength = 250,
            metric = "AUROC",
            maximize = TRUE,
            trControl = trainControl(
              method = "cv",
              number = 10,
              search = "random",
              sampling = subsampling,
              returnData = TRUE,
              returnResamp = "final",
              savePredictions = "final",
              classProbs = TRUE,
              summaryFunction = improvedSummary,
              selectionFunction = "best",
              seeds = makeCVseeds(type = "tuneLength", tunes = 250, numCV = 10),
              allowParallel = TRUE
            )
          )
        )),
        silent = TRUE
      )

    if (inherits(model_train, "try-error")) next

    # Export results

    res_ml[[i]] <-
      list(
        i = i,
        fix_seed = fix_seed,
        subsampling = subsampling,
        model_train = model_train
      )

    i <- i + 1

  }

  pb$terminate()
  gc()
  return(res_ml)

}

all_combinations <- crossing(comparisons, data_types)

results_ml <-
  map2(
    .x = all_combinations$comparisons,
    .y = all_combinations$data_types,
    .f = ~ doML(comparison = .x, type = .y)
  ) %>%
  set_names(paste0(all_combinations$comparisons, " & ", all_combinations$data_types))


# Analyze results -----------------------------------------------------------------------------

getRes <- function(comp, data_type) {
  imap_dfr(
    .x = results_ml %>%
      keep(str_detect(names(.), data_type)) %>%
      chuck(str_subset(names(.), stringr::fixed(comp))),
    .f = function(res, .y) {
      features = str_subset(colnames(res$model_train$trainingData), ".outcome", negate = TRUE)
      res$model_train$resample %>%
        as_tibble() %>%
        summarise(across(-Resample, ~ mean(.x, na.rm = TRUE))) %>%
        mutate(Model = .y, Comparison = comp, Type = data_type, .before = everything()) %>%
        mutate(DEGs_n = length(features), DEGs_Symbol = paste(features, collapse = ","))
    }
  )
}

results_ml_combined <-
  map_dfr(
    .x = comparisons,
    .f = ~ getRes(comp = .x, data_type = "\\sdelta_delta_ct")
  ) %>%
  mutate(Type = str_remove_all(Type, stringr::fixed("\\s")))

results_ml_mean <-
  map2_dfr(
    .x = all_combinations$comparisons,
    .y = all_combinations$data_types,
    .f = ~ results_ml_combined %>%
      filter(Comparison == .x & Type == .y) %>%
      mutate(across(AUROC:F1, ~ calcConfInt(.x, perc = 0.95), .names = "{.col} 95% CI")) %>%
      mutate(PrevalencePPVandNPV(Sens = Sensitivity, Spec = Specificity)) %>%
      mutate(across(c(AUROC:F1, matches("^PPV"), matches("^NPV")), ~ mean(.x, na.rm = TRUE))) %>%
      distinct(AUROC, .keep_all = TRUE) %>%
      select(
        Comparison, Type, matches("^AUROC"), matches("^Acc"), matches("^Bal"), matches("^Kappa"),
        matches("^Sens"), matches("^Spec"), matches("^F1"), matches("^PPV"), matches("^NPV"),
        DEGs_n, DEGs_Symbol
      )
  ) %>%
  mutate(DEGs_Ensembl = map(
    DEGs_Symbol, function(genes) {
      genes %>%
        str_split(",") %>%
        unlist() %>%
        AnnotationDbi::mapIds(
          x = org.Hs.eg.db::org.Hs.eg.db,
          keys = .,
          keytype = "SYMBOL",
          column = "ENSEMBL"
        ) %>%
        as.character() %>%
        paste(collapse = ",")
    }
  ))

