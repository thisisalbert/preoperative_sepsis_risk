# Title: RNA-Seq analysis with GEO2RNA-Seq data processing on Sepsis pre-surgery samples
# Project: Risk assessment with gene expression markers in sepsis development
# Mantainer: Albert Garcia-Lopez
# Last update: September 3, 2020

workdir <- "" # <-- fill in
setwd(workdir)
stopifnot(dir.exists(workdir))


# Load libraries ###############################################################

pacman::p_load(R.utils, devtools, BiocManager, tidyverse, purrr, BiocParallel, DESeq2, edgeR, limma, vroom)


# Load workspace and keep just necessary data ##################################

load("after_pcas_counting_done.RData") # output from pre-processing pipeline
rm(list = setdiff(ls(), c("counts", "gene_lengths", "lib_sizes")))
stopifnot(ncol(counts) == 267)


# Load GEO2RNASeq functions and resources ######################################

sourceDirectory("geo2rnaseq_sourcefiles")


# Load design matrix updated with the new contrast #############################

designmat <- read.csv("design_matrix.csv") %>% 
  rename(ID = X) %>% 
  mutate(ID = str_replace_all(ID, "R0+", "R")) %>% 
  column_to_rownames("ID")
stopifnot(designmat$ID == colnames(counts))


# DEG analysis #################################################################

pvalcut <- 0.05
logfcCut <- 0.25
tools <- c("DESeq2", "edgeR")
logfcnorm <- c("mrn", "tpm", "rpkm")

# 1) Inf+ vs. Inf-

dm1 <- designmat %>% select(Infection_vs_Control) %>% as.matrix()

deg_result1 <- calculate_DEGs(
  counts        = counts,
  geneLengths   = gene_lengths,
  libSizes      = lib_sizes,
  designMatrix  = dm1,
  pValCut       = pvalcut,
  logfcCut      = logfcCut,
  tools         = tools,
  logfcnorm     = logfcnorm,
  outDir        = "./Inf+_vs_Inf-/",
  prefix        = "",
  no.csv        = FALSE,
  cpus          = 1
)

deg_overview1 <- make_deg_overview_plot_UNION(
  outDir = "./Inf+_vs_Inf-/",
  degs = deg_result1$DEGs,
  tools = tools,
  fileName = paste("DEG_Overview_Inf+_vs_Inf-_deseq2_edger", Sys.Date(), ".pdf")
)

print(deg_overview1)
write.csv(deg_overview1, "./Inf+_vs_Inf-/DEG_Overview_Inf+_vs_Inf-_deseq2_edger.csv", quote = FALSE)

# 2) Inf+ vs. SIRS-
 
dm2 <- designmat %>% select(Infection_vs_Control_ODneg) %>% as.matrix()

deg_result2 <- calculate_DEGs(
  counts        = counts,
  geneLengths   = gene_lengths,
  libSizes      = lib_sizes,
  designMatrix  = dm2,
  pValCut       = pvalcut,
  logfcCut      = logfcCut,
  tools         = tools,
  logfcnorm     = logfcnorm,
  outDir        = "./Inf+_vs_SIRS-/",
  prefix        = "",
  no.csv        = FALSE,
  cpus          = 1
)

deg_overview2 <- make_deg_overview_plot_UNION(
  outDir = "./Inf+_vs_SIRS-/",
  degs = deg_result2$DEGs,
  tools = tools,
  fileName = paste("DEG_Overview_Inf+_vs_SIRS-_ODneg_deseq2_edger", Sys.Date(), ".pdf")
)

print(deg_overview2)
write.csv(deg_overview2, "./Inf+_vs_SIRS-/DEG_Overview_Inf+_vs_SIRS-_deseq2_edger.csv", quote = FALSE)

# 3) Sepsis_vs_Inf-

dm3 <- designmat %>% select(Infection_ODpos_vs_Control) %>% as.matrix()

deg_result3 <- calculate_DEGs(
  counts        = counts,
  geneLengths   = gene_lengths,
  libSizes      = lib_sizes,
  designMatrix  = dm3,
  pValCut       = pvalcut,
  logfcCut      = logfcCut,
  tools         = tools,
  logfcnorm     = logfcnorm,
  outDir        = "./Sepsis_vs_Inf-/",
  prefix        = "",
  no.csv        = FALSE,
  cpus          = 1
)

deg_overview3 <- make_deg_overview_plot_UNION(
  outDir = "./Sepsis_vs_Inf-/",
  degs = deg_result3$DEGs,
  tools = tools,
  fileName = paste("DEG_Overview_Sepsis_vs_Inf-_deseq2_edger", Sys.Date(), ".pdf")
)

print(deg_overview3)
write.csv(deg_overview3, "./Sepsis_vs_Inf-/DEG_Overview_Sepsis_vs_Inf-_deseq2_edger.csv", quote = FALSE)

# 4) Sepsis_vs_SIRS-

dm4 <- designmat %>% select(Infection_ODpos_vs_Control_ODneg) %>% as.matrix()

deg_result4 <- calculate_DEGs(
  counts        = counts,
  geneLengths   = gene_lengths,
  libSizes      = lib_sizes,
  designMatrix  = dm4,
  pValCut       = pvalcut,
  logfcCut      = logfcCut,
  tools         = tools,
  logfcnorm     = logfcnorm,
  outDir        = "./Sepsis_vs_SIRS-/",
  prefix        = "",
  no.csv        = FALSE,
  cpus          = 1
)

deg_overview4 <- make_deg_overview_plot_UNION(
  outDir = "./Sepsis_vs_SIRS-/",
  degs = deg_result4$DEGs,
  tools = tools,
  fileName = paste("DEG_Overview_Sepsis_vs_SIRS-_deseq2_edger", Sys.Date(), ".pdf")
)

print(deg_overview4)
write.csv(deg_overview4, "./Sepsis_vs_SIRS-/DEG_Overview_Sepsis_vs_SIRS-_deseq2_edger.csv", quote = FALSE)

# 5) Sepsis_vs_SIRS+

dm5 <- designmat %>% select(Infection_ODpos_vs_Control_ODpos) %>% as.matrix()

deg_result5 <- calculate_DEGs(
  counts        = counts,
  geneLengths   = gene_lengths,
  libSizes      = lib_sizes,
  designMatrix  = dm5,
  pValCut       = pvalcut,
  logfcCut      = logfcCut,
  tools         = tools,
  logfcnorm     = logfcnorm,
  outDir        = "./Sepsis_vs_SIRS+/",
  prefix        = "",
  no.csv        = FALSE,
  cpus          = 1
)

deg_overview5 <- make_deg_overview_plot_UNION(
  outDir = "./Sepsis_vs_SIRS+/",
  degs = deg_result5$DEGs,
  tools = tools,
  fileName = paste("DEG_Overview_Sepsis_vs_SIRS+_deseq2_edger", Sys.Date(), ".pdf")
)

print(deg_overview5)
write.csv(deg_overview5, "./Sepsis_vs_SIRS+/DEG_Overview_Sepsis_vs_SIRS+_deseq2_edger.csv", quote = FALSE)

# 6) Sepsis_vs_UInf+

dm6 <- designmat %>% select(Infection_ODpos_vs_Infection_ODneg) %>% as.matrix()

deg_result6 <- calculate_DEGs(
  counts        = counts,
  geneLengths   = gene_lengths,
  libSizes      = lib_sizes,
  designMatrix  = dm6,
  pValCut       = pvalcut,
  logfcCut      = logfcCut,
  tools         = tools,
  logfcnorm     = logfcnorm,
  outDir        = "./Sepsis_vs_UInf+/",
  prefix        = "",
  no.csv        = FALSE,
  cpus          = 1
)

deg_overview6 <- make_deg_overview_plot_UNION(
  outDir = "./Sepsis_vs_UInf+/",
  degs = deg_result6$DEGs,
  tools = tools,
  fileName = paste("DEG_Overview_Sepsis_vs_UInf+_deseq2_edger", Sys.Date(), ".pdf")
)

print(deg_overview6)
write.csv(deg_overview6, "./Sepsis_vs_UInf+/DEG_Overview_Sepsis_vs_UInf+_deseq2_edger.csv", quote = FALSE)

# 7) Sepsis+_vs_Sepsis-

dm7 <- designmat %>% select(Infection_ODpos_vs_Infection_ODneg_and_Control) %>% as.matrix()

deg_result7 <- calculate_DEGs(
  counts        = counts,
  geneLengths   = gene_lengths,
  libSizes      = lib_sizes,
  designMatrix  = dm7,
  pValCut       = pvalcut,
  logfcCut      = logfcCut,
  tools         = tools,
  logfcnorm     = logfcnorm,
  outDir        = "./Sepsis+_vs_Sepsis-/",
  prefix        = "",
  no.csv        = FALSE,
  cpus          = 1
)

deg_overview7 <- make_deg_overview_plot_UNION(
  outDir = "./Sepsis+_vs_Sepsis-/",
  degs = deg_result7$DEGs,
  tools = tools,
  fileName = paste("DEG_Overview_Sepsis+_vs_Sepsis-_deseq2_edger", Sys.Date(), ".pdf")
)

print(deg_overview7)
write.csv(deg_overview7, "./Sepsis+_vs_Sepsis-/DEG_Overview_Sepsis+_vs_Sepsis-_deseq2_edger.csv", quote = FALSE)

# 8) Inf+_vs_SIRS+

dm8 <- designmat %>% select(Infection_vs_Control_ODpos) %>% as.matrix()

deg_result8 <- calculate_DEGs(
  counts        = counts,
  geneLengths   = gene_lengths,
  libSizes      = lib_sizes,
  designMatrix  = dm8,
  pValCut       = pvalcut,
  logfcCut      = logfcCut,
  tools         = tools,
  logfcnorm     = logfcnorm,
  outDir        = "./Inf+_vs_SIRS+/",
  prefix        = "",
  no.csv        = FALSE,
  cpus          = 1
)

deg_overview8 <- make_deg_overview_plot_UNION(
  outDir = "./Inf+_vs_SIRS+/",
  degs = deg_result8$DEGs,
  tools = tools,
  fileName = paste("DEG_Overview_Inf+_vs_SIRS+_deseq2_edger", Sys.Date(), ".pdf")
)

print(deg_overview8)
write.csv(deg_overview8, "./Inf+_vs_SIRS+/DEG_Overview_Inf+_vs_SIRS+_deseq2_edger.csv", quote = FALSE)

# 9) OD+_vs_OD-

dm9 <- designmat %>% select(ODpos_vs_ODneg) %>% as.matrix()

deg_result9 <- calculate_DEGs(
  counts        = counts,
  geneLengths   = gene_lengths,
  libSizes      = lib_sizes,
  designMatrix  = dm9,
  pValCut       = pvalcut,
  logfcCut      = logfcCut,
  tools         = tools,
  logfcnorm     = logfcnorm,
  outDir        = "./OD+_vs_OD-/",
  prefix        = "",
  no.csv        = FALSE,
  cpus          = 1
)

deg_overview9 <- make_deg_overview_plot_UNION(
  outDir = "./OD+_vs_OD-/",
  degs = deg_result9$DEGs,
  tools = tools,
  fileName = paste("DEG_Overview_OD+_vs_OD-_deseq2_edger", Sys.Date(), ".pdf")
)

print(deg_overview9)
write.csv(deg_overview9, "./OD+_vs_OD-/DEG_Overview_OD+_vs_OD-_deseq2_edger.csv", quote = FALSE)


# Save table with all DEGs for all contrasts ###################################

deg_overview_all <- rbind(
  deg_overview1,
  deg_overview2,
  deg_overview3,
  deg_overview4,
  deg_overview5,
  deg_overview6,
  deg_overview7,
  deg_overview8,
  deg_overview9
)

write.csv(deg_overview_all, "./table_degs_all_contrasts_counts100.csv")


# Save results #################################################################

save.image(paste0(Sys.Date(), "_pval005_logFC025_deseq2_edger.RData"))
