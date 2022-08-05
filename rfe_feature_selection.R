# Title: Machine Learning Recursive Feature Elimination (RFE)
# Project: Risk assessment with gene expression markers in sepsis development
# Mantainer: Albert Garcia-Lopez
# Last update: June 9, 2022

outDir <- "" # <-- fill in
stopifnot(dir.exists(outDir))
pacman::p_load(tidyverse, devtools, caret, parallel, doParallel, readxl)
options(readr.show_col_types = FALSE)


# Load counts reads ############################################################

load("after_pcas_counting_done.RData") # # output from pre-processing pipeline
rm(list = setdiff(ls(), c("counts", "outDir")))
stopifnot(ncol(counts) == 267)


# Load custom functions ########################################################

source_url("https://github.com/thisisalbert/machine_learning/raw/main/customSummary.R")
source_url("https://github.com/thisisalbert/machine_learning/raw/main/makeRFEseeds.R")
source_url("https://github.com/thisisalbert/bioinformatics/raw/main/saveWorkspace.R")


# Load 64 DEG subset ###########################################################

degs <- read_excel("Final_List_Genes_64_Annotated_20_11_26.xlsx") %>%
  pull(Ensembl_Gene_ID)


# Get counts from the 64 DEGs only #############################################

counts_degs <- counts %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  filter(gene_id %in% degs) %>%
  column_to_rownames("gene_id") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  as_tibble()

stopifnot(all(degs %in% colnames(counts_degs)))


# Load metadata files ##########################################################

meta_infpos_infneg <- read_tsv("infpos_infneg.tsv")
meta_infpos_sirsneg <- read_tsv("infpos_sirsneg.tsv")
meta_infpos_sirspos <- read_tsv("infpos_sirspos.tsv")
meta_sepsispos_infneg <- read_tsv("sepsispos_infneg.tsv")
meta_sepsispos_sirsneg <- read_tsv("sepsispos_sirsneg.tsv")
meta_sepsispos_sirspos <- read_tsv("sepsispos_sirspos.tsv")
meta_sepsispos_uinfpos <- read_tsv("sepsispos_uinfpos.tsv")
meta_sepsispos_sepsisneg <- read_tsv("sepsispos_sepsisneg.tsv")


# Sequence of seeds ############################################################

seq_seeds <- sample.int(n = 2^15, size = 2^15, replace = FALSE)


# RFE function #################################################################

RFEprocedure <- function(meta) {

  rfe_res <- list()
  i <- 1

  while (i <= 10) {

    # Parallelization

    registerDoParallel(cores = detectCores() - 5)

    # Start iteration

    message(paste0("Model ", i))
    fix_seed <- sample(seq_seeds, 1)
    seq_seeds <- seq_seeds[seq_seeds != fix_seed]

    # Define case and control

    case <- meta$Comparison %>%
      unique() %>%
      str_subset(pattern = "inf_pos|sepsis_pos")

    control <- meta$Comparison %>%
      unique() %>%
      str_subset(pattern = "inf_pos|sepsis_pos", negate = TRUE)

    # Balance genders

    num_females <- meta %>%
      filter(Gender == "Female") %>%
      nrow()

    set.seed(fix_seed)
    balanced_metadata <- meta %>%
      group_by(Gender) %>%
      slice_sample(n = num_females, replace = FALSE) %>%
      ungroup()

    stopifnot(all(duplicated(balanced_metadata$ID)) == FALSE)
    stopifnot(identical(
      balanced_metadata %>% filter(Gender == "Female") %>% nrow(),
      balanced_metadata %>% filter(Gender == "Male") %>% nrow()
    ))

    # Mix balanced metadata with counts

    mix_counts_meta <- left_join(
      x = balanced_metadata,
      y = counts_degs,
      by = "ID"
    ) %>%
      mutate(Comparison = factor(Comparison, c(case, control))) %>%
      as.data.frame()

    # Define function for model training

    customRFE <- caretFuncs
    customRFE$summary <- customSummary

    # RFE

    set.seed(fix_seed)
    model_rfe <- tryCatch(suppressWarnings(
      caret::rfe(
        x = mix_counts_meta %>% select(all_of(degs)),
        y = mix_counts_meta %>% chuck("Comparison"),
        sizes = seq(1, length(degs) - 1),
        maximize = TRUE,
        method = "svmLinear",
        metric = "ROC",
        preProcess = c("center", "scale", "nzv"),
        rfeControl = rfeControl(
          functions = customRFE,
          rerank = FALSE,
          method = "cv",
          number = 10,
          returnResamp = "final",
          allowParallel = TRUE,
          seeds = NULL
        ),
        trControl = trainControl(
          search = "random",
          sampling = "down",
          classProbs = TRUE,
          summaryFunction = customSummary
        )
      )
    ), error = function(e) {e})

    if (inherits(model_rfe, "error")) {next}
    if (all(model_rfe$results$ROC < 0.8)) {next}

    # Export result

    rfe_res[[i]] <- tibble(
      model = i,
      fix_seed = fix_seed,
      mix_counts_meta = list(mix_counts_meta),
      model_rfe = list(model_rfe)
    )

    i <- i + 1

  }

  return(rfe_res)

}


# Perform RFE on each contrast #################################################

rfe_infpos_infneg <- RFEprocedure(meta_infpos_infneg)
rfe_infpos_sirsneg <- RFEprocedure(meta_infpos_sirsneg)
rfe_infpos_sirspos <- RFEprocedure(meta_infpos_sirspos)
rfe_sepsispos_infneg <- RFEprocedure(meta_sepsispos_infneg)
rfe_sepsispos_sirsneg <- RFEprocedure(meta_sepsispos_sirsneg)
rfe_sepsispos_sirspos <- RFEprocedure(meta_sepsispos_sirspos)
rfe_sepsispos_uinfpos <- RFEprocedure(meta_sepsispos_uinfpos)
rfe_sepsispos_sepsisneg <- RFEprocedure(meta_sepsispos_sepsisneg)


# Export RFE objects ###########################################################

write_rds(rfe_infpos_infneg, paste0(outDir, "rfe_infpos_infneg.rds"))
write_rds(rfe_infpos_sirsneg, paste0(outDir, "rfe_infpos_sirsneg.rds"))
write_rds(rfe_infpos_sirspos, paste0(outDir, "rfe_infpos_sirspos.rds"))
write_rds(rfe_sepsispos_infneg, paste0(outDir, "rfe_sepsispos_infneg.rds"))
write_rds(rfe_sepsispos_sirsneg, paste0(outDir, "rfe_sepsispos_sirsneg.rds"))
write_rds(rfe_sepsispos_sirspos, paste0(outDir, "rfe_sepsispos_sirspos.rds"))
write_rds(rfe_sepsispos_uinfpos, paste0(outDir, "rfe_sepsispos_uinfpos.rds"))
write_rds(rfe_sepsispos_sepsisneg, paste0(outDir, "rfe_sepsispos_sepsisneg.rds"))


# Save workspace ###############################################################

saveWorkspace()
