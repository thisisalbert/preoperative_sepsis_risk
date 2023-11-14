# Title: Feature selection
# Project: Risk assessment with gene expression markers in sepsis development
# Maintainer: Albert Garcia-Lopez
# Last update: June 09, 2023


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


# Load TPM counts -----------------------------------------------------------------------------

tpm_counts <- import("tpm_counts.rds")


# Load metadata -------------------------------------------------------------------------------

metadata <- import("metadata.xlsx")


# Create metadata files for each comparison ---------------------------------------------------

comparisons <- c(
  "Inf+ vs. Inf-", "Inf+ vs. SIRS-", "Inf+ vs. SIRS+",
  "Sepsis+ vs. Inf-", "Sepsis+ vs. SIRS-", "Sepsis+ vs. SIRS+",
  "Sepsis+ vs. UInf+"
)

metadata_extended <-
  metadata %>%
  select(id, sex, class, organ_dysfunction) %>%
  mutate(
    `Inf+ vs. Inf-` = case_when(
      class %in% c("Sepsis+", "UInf+") ~ "inf_pos",
      class %in% c("SIRS+", "SIRS-") ~ "inf_neg",
      TRUE ~ NA
    ),
    `Inf+ vs. SIRS-` = case_when(
      class %in% c("Sepsis+", "UInf+") ~ "inf_pos",
      class == "SIRS-" ~ "sirs_neg",
      TRUE ~ NA
    ),
    `Inf+ vs. SIRS+` = case_when(
      class %in% c("Sepsis+", "UInf+") ~ "inf_pos",
      class == "SIRS+" ~ "sirs_pos",
      TRUE ~ NA
    ),
    `Sepsis+ vs. Inf-` = case_when(
      class == "Sepsis+" ~ "sepsis_pos",
      class %in% c("SIRS+", "SIRS-") ~ "inf_neg",
      TRUE ~ NA
    ),
    `Sepsis+ vs. SIRS-` = case_when(
      class == "Sepsis+" ~ "sepsis_pos",
      class == "SIRS-" ~ "sirs_neg",
      TRUE ~ NA
    ),
    `Sepsis+ vs. SIRS+` = case_when(
      class == "Sepsis+" ~ "sepsis_pos",
      class == "SIRS+" ~ "sirs_pos",
      TRUE ~ NA
    ),
    `Sepsis+ vs. UInf+` = case_when(
      class == "Sepsis+" ~ "sepsis_pos",
      class == "UInf+" ~ "uinf_pos",
      TRUE ~ NA
    )
  )


# Load DEGs -----------------------------------------------------------------------------------

degs_all <-
  bind_rows(
    import("df_degs_general.xlsx"),
    import("df_degs_male.xlsx"),
    import("df_degs_female.xlsx")
  ) %>%
  as_tibble() %>%
  pull(id) %>%
  unique()


# Select only protein-coding DEGs -------------------------------------------------------------

pacman::p_load(biomaRt)
set.seed(123)
degs_all_protcod <-
  getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id", "gene_biotype"),
    filters = "ensembl_gene_id",
    values = degs_all,
    mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  ) %>%
  as_tibble() %>%
  filter(gene_biotype == "protein_coding") %>%
  mutate(across(everything(), ~ ifelse(.x == "", NA, .x)))
pacman::p_unload("biomaRt")


# Generate sequence of seeds ------------------------------------------------------------------

seq_seeds <- sample.int(n = 2^15, size = 2^15, replace = FALSE)


# Boruta + Train Machine Learning -------------------------------------------------------------

doBoruta <- function(comparison) {

  # Define iteration

  res_boruta <- list()
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
      filter(rownames(.) %in% degs_all_protcod$ensembl_gene_id) %>%
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
      all(str_subset(colnames(metadata_balanced_with_counts), "^ENSG") %in% degs_all_protcod$ensembl_gene_id) &
      all(as.character(unique(metadata_balanced_with_counts$contrast)) %in% c(case, control)) &
      is.data.frame(metadata_balanced_with_counts)
    )

    # Boruta

    set.seed(fix_seed)
    boruta_obj <-
      try(
        expr = Boruta(
          x = metadata_balanced_with_counts %>% select(all_of(degs_all_protcod$ensembl_gene_id)),
          y = metadata_balanced_with_counts %>% pull(contrast)
        ),
        silent = TRUE
      )

    if (inherits(boruta_obj, "try-error")) next
    if (length(getSelectedAttributes(boruta_obj)) == 0) next

    boruta_df <-
      boruta_obj %>%
      chuck("ImpHistory") %>%
      as_tibble() %>%
      select(getSelectedAttributes(boruta_obj)) %>%
      filter(if_all(everything(), ~ is.finite(.x))) %>%
      summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
      pivot_longer(everything(), names_to = "Genes", values_to = "Boruta_Score") %>%
      arrange(desc(Boruta_Score)) %>%
      slice_max(order_by = Boruta_Score, n = 25, with_ties = FALSE)

    # Training with Boruta-selected genes

    set.seed(fix_seed)
    model_train <-
      try(
        expr = suppressWarnings(
          caret::train(
            x = metadata_balanced_with_counts %>% select(all_of(boruta_df$Genes)),
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

    if (inherits(model_train, "try-error")) next

    # Export results

    res_boruta[[i]] <-
      list(
        i = i,
        fix_seed = fix_seed,
        metadata_balanced_with_counts = metadata_balanced_with_counts,
        case = case,
        control = control,
        boruta_df = boruta_df,
        model_train = model_train
      )

    i <- i + 1

  }

  close(pb)
  stopImplicitCluster()
  return(res_boruta)

}

results_boruta <-
  map(
    .x = comparisons,
    .f = ~ doBoruta(comparison = .x)
  ) %>%
  set_names(comparisons)


# Analyze results -----------------------------------------------------------------------------

results_boruta_metrics <-
  map(
    .x = comparisons,
    .f = function(comp) {
      map_dfr(
        .x = results_boruta[[comp]],
        .f = function(res) {
          res %>%
            chuck("model_train") %>%
            chuck("resample") %>%
            summarise(across(-Resample, ~ mean(.x, na.rm = TRUE))) %>%
            as_tibble() %>%
            mutate(
              Comparison = comp,
              Iteration = res$i,
              Seed = res$fix_seed,
              Boruta_DEGs_n = length(predictors(res$model_train)),
              Boruta_DEGs = paste(predictors(res$model_train), collapse = ","),
              .before = everything()
            ) %>%
            relocate(Boruta_DEGs_n, Boruta_DEGs, .after = last_col())
        }
      )
    }
  ) %>%
  set_names(comparisons)

# Select the 10 iterations with highest AUC on training cross-validation and get DEGs for 
# downstream use.

set.seed(123)
selection_boruta <-
  map(
    .x = comparisons,
    .f = function(comp) {
      results_boruta_metrics[[comp]] %>%
        slice_max(order_by = AUROC, n = 10, with_ties = FALSE) %>%
        pull(Boruta_DEGs) %>%
        str_split(",") %>%
        unlist() %>%
        unique()
    }
  ) %>%
  set_names(comparisons)

map_dfr(
  .x = results_boruta_metrics,
  .f = ~ .x %>% slice_max(order_by = AUROC, n = 10, with_ties = FALSE)
)


# RFE Machine Learning ------------------------------------------------------------------------

doRFE <- function(comparison) {

  res_rfe <- list()
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

    # Load DEGs

    degs <- selection_boruta[[comparison]]

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

    # RFE

    set.seed(fix_seed)
    model_rfe <-
      try(
        expr = suppressWarnings(
          caret::rfe(
            x = metadata_balanced_with_counts %>% select(all_of(degs)),
            y = metadata_balanced_with_counts %>% pull(contrast),
            sizes = seq(from = 1, to = length(degs) - 1, by = 1),
            metric = "AUROC",
            maximize = TRUE,
            rfeControl = rfeControl(
              functions = rangerFuncs,
              rerank = FALSE,
              method = "cv",
              number = 10,
              saveDetails = TRUE,
              returnResamp = "final",
              seeds = makeRFEseeds(features = length(degs), numCV = 10),
              allowParallel = TRUE
            )
          )
        ),
        silent = TRUE
      )

    if (inherits(model_rfe, "try-error")) next
    if (max(model_rfe$results$AUROC) < 0.8) next

    # Export

    res_rfe[[i]] <-
      list(
        i = i,
        fix_seed = fix_seed,
        metadata_balanced_with_counts = metadata_balanced_with_counts,
        case = case,
        control = control,
        degs = degs,
        model_rfe = model_rfe
      )

    i <- i + 1

  }

  close(pb)
  stopImplicitCluster()
  return(res_rfe)

}

results_rfe <-
  map(
    .x = comparisons,
    .f = ~ doRFE(comparison = .x)
  ) %>%
  set_names(comparisons)

