# Title: DEG quantification
# Project: Risk assessment with gene expression markers in sepsis development
# Maintainer: Albert Garcia-Lopez
# Last update: June 01, 2023


# Define output directory and libraries -------------------------------------------------------

outDir <- dirname(rstudioapi::getActiveDocumentContext()$path)
stopifnot(dir.exists(outDir) & file.exists(rstudioapi::getActiveDocumentContext()$path))
pacman::p_load(tidyverse, rio, Geo2RNAseq)
devtools::source_url("https://github.com/thisisalbert/AlbertML/raw/main/R/saveWorkspace.R")


# Make design matrix --------------------------------------------------------------------------

dm_general <-
  import("metadata.xlsx") %>%
  as_tibble() %>%
  select(id, sex, class, organ_dysfunction) %>%
  mutate(
    `Inf+ vs. Inf-` = case_when(
      class %in% c("Sepsis+", "UInf+") ~ "treatment",
      class %in% c("SIRS+", "SIRS-") ~ "control",
    ),
    `Inf+ vs. SIRS-` = case_when(
      class %in% c("Sepsis+", "UInf+") ~ "treatment",
      class == "SIRS-" ~ "control",
      TRUE ~ "none"
    ),
    `Inf+ vs. SIRS+` = case_when(
      class %in% c("Sepsis+", "UInf+") ~ "treatment",
      class == "SIRS+" ~ "control",
      TRUE ~ "none"
    ),
    `Sepsis+ vs. Inf-` = case_when(
      class == "Sepsis+" ~ "treatment",
      class %in% c("SIRS+", "SIRS-") ~ "control",
      TRUE ~ "none"
    ),
    `Sepsis+ vs. SIRS-` = case_when(
      class == "Sepsis+" ~ "treatment",
      class == "SIRS-" ~ "control",
      TRUE ~ "none"
    ),
    `Sepsis+ vs. SIRS+` = case_when(
      class == "Sepsis+" ~ "treatment",
      class == "SIRS+" ~ "control",
      TRUE ~ "none"
    ),
    `Sepsis+ vs. UInf+` = case_when(
      class == "Sepsis+" ~ "treatment",
      class == "UInf+" ~ "control",
      TRUE ~ "none"
    )
  )


# Make male and female design matrices --------------------------------------------------------

dm_male <-
  dm_general %>%
  filter(sex == "Male") %>%
  select(-sex, -class, -organ_dysfunction) %>%
  column_to_rownames("id") %>%
  as.matrix()

dm_female <-
  dm_general %>%
  filter(sex == "Female") %>%
  select(-sex, -class, -organ_dysfunction) %>%
  column_to_rownames("id") %>%
  as.matrix()

dm_general <-
  dm_general %>%
  select(-sex, -class, -organ_dysfunction) %>%
  column_to_rownames("id") %>%
  as.matrix()


# Sanity check --------------------------------------------------------------------------------

stopifnot(
  is.matrix(dm_general) & 
  is.matrix(dm_female) & 
  is.matrix(dm_male) & 
  is.matrix(counts)
)

stopifnot(
  identical(rownames(dm_general), colnames(counts[, rownames(dm_general)])) &
  identical(rownames(dm_female), colnames(counts[, rownames(dm_female)])) &
  identical(rownames(dm_male), colnames(counts[, rownames(dm_male)]))
)


# DEG quantification --------------------------------------------------------------------------

# Settings

pValCut <- 0.05
logfcCut <- 0.25
tools <- c("DESeq2", "edgeR")
logfcnorm <- c("mrn", "tpm", "rpkm")

# Results directory

generalDir <- file.path(outDir, "generalDir")
maleDir <- file.path(outDir, "maleDir")
femaleDir <- file.path(outDir, "femaleDir")

map(
  .x = c(generalDir, maleDir, femaleDir),
  .f = function(dir) {
    if (dir.exists(dir)) {
      unlink(dir, recursive = TRUE, force = TRUE)
      dir.create(dir)
    } else {
      dir.create(dir)
    }
  }
)

# General (both sexes) results

degs_general <-
  calculate_DEGs(
    counts = counts,
    geneLengths = gene_lengths,
    libSizes = lib_sizes,
    designMatrix = dm_general,
    pValCut = pValCut,
    logfcCut = logfcCut,
    logfcnorm = logfcnorm,
    tools = tools,
    outDir = generalDir,
    cpus = 1,
    prefix = ""
  )

# Male results

degs_male <-
  calculate_DEGs(
    counts = counts[, rownames(dm_male)],
    geneLengths = gene_lengths,
    libSizes = lib_sizes[rownames(dm_male)],
    designMatrix = dm_male,
    pValCut = pValCut,
    logfcCut = logfcCut,
    logfcnorm = logfcnorm,
    tools = tools,
    outDir = maleDir,
    cpus = 1,
    prefix = ""
  )

# Female results

degs_female <-
  calculate_DEGs(
    counts = counts[, rownames(dm_female)],
    geneLengths = gene_lengths,
    libSizes = lib_sizes[rownames(dm_female)],
    designMatrix = dm_female,
    pValCut = pValCut,
    logfcCut = logfcCut,
    logfcnorm = logfcnorm,
    tools = tools,
    outDir = femaleDir,
    cpus = 1,
    prefix = ""
  )

