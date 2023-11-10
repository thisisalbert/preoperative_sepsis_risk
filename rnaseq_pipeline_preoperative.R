# Title: Preoperative RNA-seq preprocessing (Script based on GEO2RNAseq pipeline version R320)
# Project: Risk assessment with gene expression markers in sepsis development
# Maintainer: Albert Garcia-Lopez
# Last update: July 22, 2020

# Load libraries ###############################################################

library("Geo2RNAseq")
library("R.utils")
library("tidyverse")
library("ggplot2")
library("data.table")
library("gtools")
library("stringi")
# R.utils::sourceDirectory("R/") # only if not installed via install.packages()


# Set WD #######################################################################

workdir <- "" # <-- fill in
setwd(workdir)
stopifnot(dir.exists(getwd()))


# Global settings ##############################################################

# Output directories

fastqDir <- "./fastq"
qualDir  <- "./quality"
mapDir   <- "./mapping"
countDir <- "./counting"
tabDir   <- "./result_tables"
plotDir  <- "./result_plots"
degDir   <- "./diff_exp_genes"

# If TRUE, existing files will be overwritten without questioning.
# If FALSE, most methods will skip step for existing output files.
# In these cases, the method will return 'call[x] = "NOT USED"'

FORCE_OVERWRITE <- FALSE
paired <- TRUE
MAX_CPUS <- 30


# Files ########################################################################

genome <- "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
anno   <- "Homo_sapiens.GRCh38.99.gtf"
stopifnot(file.exists(genome))
stopifnot(file.exists(anno))

# Build index with make_HiSat2_index() or make_Tophat_index() ##################

index  <- make_HiSat2_index("Homo_sapiens.GRCh38.dna.primary_assembly.fa")

# Checkpoint

save.image("./index_done.RData")


# Raw FASTQ files ##############################################################

raw_fastq_files <- list.files("raw_data", full.names = TRUE, pattern = "\\.gz")
raw_fastq_files <- raw_fastq_files[grep("md5", raw_fastq_files, invert = TRUE)]

if (paired == TRUE) {
  print(length(raw_fastq_files)/2)
} else {
  print(length(raw_fastq_files))
}

if (!file.exists(raw_fastq_files[1])){
  raw_fastq_files <- file.path(fastqDir, raw_fastq_files)
}

writeLines("Working with files:")
print(raw_fastq_files)


# Quality control 1 (before trimming) ##########################################

qualityRawDir <- file.path(qualDir, "raw")
if (length(dir(qualityRawDir)) > 0) {
    warning(paste0("Directory \"", qualityRawDir, "\" not empty! Overwrite? ", FORCE_OVERWRITE), immediate. = TRUE)
}

if (FORCE_OVERWRITE || length(dir(qualityRawDir)) == 0) {
    writeLines("FastQC - raw ...")
    fq_res <- run_FastQC(raw_fastq_files, outDir = qualityRawDir, cpus = MAX_CPUS, extend = TRUE)
}


# Trimming #####################################################################

trimmedDir <- fastqDir
windowsizetrimming <- 15
qualcuttrimming <- 25
phred <- "-phred33"
leading <- 3
trailing <- 3
minlen <- 30

if (length(list.files(trimmedDir, pattern = "\\.trimo(\\.pe)*\\.fastq$")) > 0) {
  warning(paste0("Directory \"", trimmedDir, "\" not empty! Overwrite? ", FORCE_OVERWRITE), immediate. = TRUE)
}

trimming_res <- run_Trimmomatic(
    files      = raw_fastq_files,
    is.paired  = paired,
    outDir     = trimmedDir,
    cpus       = MAX_CPUS,
    windowsize = windowsizetrimming,
    qualcut    = qualcuttrimming,
    phred      = phred,
    leading    = leading,
    trailing   = trailing,
    minlen     = minlen
)

# Number of reads after trimming based on Trimmomatic output
number_raw     <- trimming_res$input
number_trimmed <- trimming_res$surviving

if (paired) {
    trimmed_fastq_files <- asPairVector(trimming_res$files)
} else {
    trimmed_fastq_files <- trimming_res$files
}

names(number_raw)     <- basename(raw_fastq_files)
names(number_trimmed) <- basename(trimmed_fastq_files)
fastq_files <- trimmed_fastq_files


# SortMeRNA (removes rRNA reads) ###############################################

filterrRNA <- TRUE
sortmeDir <- file.path(fastqDir, "sortmerna")

if (filterrRNA && length(dir(sortmeDir)) > 0) {
  warning(paste0("Directory \"", sortmeDir, "\" not empty! Overwrite? ", FORCE_OVERWRITE), immediate. = TRUE)
}

if (filterrRNA) {
  sortmerna_res <- run_SortMeRNA(
    fastq_files,
    outDir    = fastqDir,
    mode      = "fast",
    paired    = paired,
    cpus      = MAX_CPUS,
    overwrite = TRUE
  )
  non_rrna_files <- sortmerna_res$files
  number_nonrRNA <- sortmerna_res$nonrrna
}

# If sample/file order changed:

stopifnot(unique(order(trimmed_fastq_files) == order(non_rrna_files)) == TRUE)
fastq_files <- if (filterrRNA) non_rrna_files else trimmed_fastq_files


# Quality control 2 (after trimming) ###########################################

qualitySubDir <- if (filterrRNA) file.path(qualDir, "non_rrna") else file.path(qualDir, "trimmed")
if (length(dir(qualitySubDir)) > 0) {
    warning(paste0("Directory \"" , qualitySubDir, "\" not empty! overwrite? ", FORCE_OVERWRITE))
}

if (FORCE_OVERWRITE || length(dir(qualitySubDir)) == 0) {
    writeLines("FastQC - trimmed ...")
    fq_res <- run_FastQC(fastq_files, outDir = qualitySubDir, cpus = MAX_CPUS, extend = FALSE)
}


# Mapping ######################################################################

index_map <- index
anno_map  <- anno
mapper <- "hisat2"
bamDir <- file.path(mapDir, "bamfiles")
samDir <- file.path(mapDir, "samfiles")
convert_sam <- TRUE # Should SAM files be converted to BAM files if BAM files are not found?

if (mapper == "tophat2") {
    mapping_res <- run_Tophat(
        fastq_files,
        index        = index_map,
        outDir       = mapDir,
        is.paired    = paired,
        anno         = anno_map,
        addArgs      = NA,      # "-g 2 --b2-very-sensitive --no-coverage-search",
        cpus         = MAX_CPUS,
        worker       = if (paired) 5 else 10,
        use.existing = T
    )
    bam_files <- mapping_res$files

} else {

    if (FORCE_OVERWRITE || length(dir(samDir)) == 0) {
        mapping_res <- run_Hisat2(
            files     = fastq_files,
            index     = index_map,
            outDir    = mapDir,
            is.paired = paired,
            addArgs   = "",
            splice    = T,
            cpus      = MAX_CPUS
        )

    ## NOTE: use that only if the run_Hisat2() parameter 'as.bam' is FALSE.
    # bam_files <- sam_to_bam(mapping_res$mapFiles, sort=TRUE, bamDir = bamDir)
    bam_files <- mapping_res$files

    } else {
        if (convert_sam) {
            sam_files <- list.files(samDir, pattern = "\\.sam$", full.names = TRUE)
            if (length(sam_files) != number_samples)
                stop(paste("Invalid number of SAM files. Expected", number_samples, "but got", length(sam_files)))
            bam_files <- sam_to_bam(mapping_res$mapFiles, sort=TRUE, bamDir = bamDir)
        } else
            bam_files <- list.files(bamDir, pattern = "\\.bam$", full.names = TRUE)

        if (length(bam_files) != number_samples)
            warning(paste("Found", length(bam_files),
            "BAM files, but expected", number_samples,
            ". Maybe SAMtools was interrupted. Try to convert again? ",
            convert_sam), immediate. = TRUE
        )
    }
}


# Counting #####################################################################

anno_count <- anno

if (length(dir(countDir)) > 0)
    warning(paste("Directory", countDir, "not empty! Overwrite?", FORCE_OVERWRITE), immediate. = TRUE)

if (FORCE_OVERWRITE || length(dir(countDir)) == 0) {
    res_counting <- run_featureCounts(
        files             = bam_files,
        annotation        = anno_count,
        isGTF             = TRUE,
        IDtype            = "gene_id",
        featureType       = "exon",
        outDir            = countDir,
        isPairedEnd       = paired,
        allowMultiOverlap = FALSE,
        cpus              = MAX_CPUS
    )

} else {
    res_counting <- read.featureCounts.files(countDir)
}

counts <- res_counting$counts
countfile <- res_counting$countFile
sumfile <- res_counting$sumFile

# Change the names of samples names in counts and lib_sizes:

colnames(counts) <- gsub("_.*", "", colnames(counts))
info <- res_counting$anno
lib_sizes <- res_counting$summary[,1]
gene_lengths <- info$Length
names(gene_lengths) <- info$GeneID
names(lib_sizes) <- gsub("_.*", "", names(lib_sizes))
samples <- colnames(counts)

# Write counts to CSV file

write_count_table(
    file   = file.path(tabDir, "counts"),
    counts = counts
)

# Write counts to XLS file

write_count_table(
    file       = file.path(tabDir, "counts"),
    counts     = counts,
    as.xls     = TRUE,
    sheetNames = basename(getwd())
)


# Check for rRNA contamination #################################################

detect_high_coverage(counts, lib_sizes)

rrna_contaminated_transcripts <- read.table(
  "./ordered_transcripts_rRNA_contaminated.txt"
) %>%
  unlist() %>%
  as.character() %>%
  unique()

filtered_counts <- counts[-which(rownames(counts) %in% rrna_contaminated_transcripts),]
detect_high_coverage(filtered_counts, lib_sizes)
rm(filtered_counts)

counts <- counts[-which(rownames(counts) %in% rrna_contaminated_transcripts),]
detect_high_coverage(counts, lib_sizes) # Works!


# SAMtools #####################################################################

flagstatDir <- file.path(mapDir, "flagstats")
flag_files <- make_flagstats(bam_files, flagstatDir, MAX_CPUS)


# MultiQC ######################################################################

muliqc_config_yaml <- file.path(
  '/sbidata/seelbind/Geo2RNAseq/hki_rna_seq/inst/extdata/multiqc_config.yaml'
)

run_MultiQC(
    tools = c(
        "fastqc",
        "trimmomatic",
        "sortmerna",
        "tophat",
        "hisat2",
        "samtools",
        "featureCounts"
    ),
    config = muliqc_config_yaml,
    force = FORCE_OVERWRITE
)


# Mapping stats ################################################################

if (!exists("number_trimmed")) {
  warning("Variable 'number_trimmed' undefined. Setting it to NA.")
  number_trimmed <- NA
}
if (!exists("number_nonrRNA")) {
  warning("Variable 'number_nonrRNA' undefined. Setting it to NA.")
  number_nonrRNA <- NA
}

# "precise" mapping stats requires BAM sorting.
# If TRUE, Bioconductor functions are used to determine mapping stats.

precise <- FALSE

mapping_stats_df <- calc_mapping_stats(
    bamFiles   = bam_files,
    fqFiles    = raw_fastq_files,
    anno       = anno_map,       # annotation used for mapping ~ can be different from counting!
    numReads   = number_raw,     # number of raw reads - mandatory!
    numTrimmed = number_trimmed, # number of remaining reads after trimming - or NA
    numNonrRNA = number_nonrRNA, # number of remaining reads after SortMeRNA - or NA
    libSizes   = lib_sizes,
    paired     = paired,
    precise    = precise,
    remove.na  = FALSE,
    cpus       = MAX_CPUS
)

# Save mapping stats

write_count_table(
    file       = file.path(tabDir, "mapping_stats"),
    counts     = mapping_stats_df,
    as.xls     = TRUE,
    sheetNames = basename(getwd()),
    rnames     = TRUE
)
write_count_table(
    file       = file.path(tabDir, "mapping_stats"),
    counts     = mapping_stats_df,
    as.xls     = FALSE,
    sheetNames = basename(getwd()),
    rnames     = TRUE
)

# Sanity check

if (paired) {
    if (FALSE %in% (asPaired(number_raw) >= lib_sizes)) {
        stop("Error: number of reads mapping in exons should not exceed original number of reads!")
    }
} else {
    if (FALSE %in% (number_raw >= lib_sizes)) {
        stop("Error: number of reads mapping in exons should not exceed original number of reads!")
    }
}


# Normalized count values ######################################################

rpkm <- get_rpkm(counts, gene_lengths, lib_sizes)
tpm  <- get_tpm(counts, gene_lengths, lib_sizes)
mrn  <- get_mrn(counts)

# Write values to CSV file

write_count_table(
    file   = file.path(tabDir, "rpkm"),
    counts = rpkm
)
write_count_table(
    file   = file.path(tabDir, "tpm"),
    counts = tpm
)
write_count_table(
    file   = file.path(tabDir, "mrn"),
    counts = mrn
)

# Write to XLS file

write_count_table(
    file       = file.path(tabDir, "rpkm"),
    counts     = rpkm,
    as.xls     = TRUE,
    sheetNames = basename(getwd())
)
write_count_table(
    file       = file.path(tabDir, "tpm"),
    counts     = tpm,
    as.xls     = TRUE,
    sheetNames = basename(getwd())
)
write_count_table(
    file       = file.path(tabDir, "mrn"),
    counts     = mrn,
    as.xls     = TRUE,
    sheetNames = basename(getwd())
)


# Design Matrix ################################################################

# Note: the design matrix needs to be a "matrix" object for the downstream
# analysis to work!

design_matrix <- read.csv("./design_matrix_267_samples_final_edited.csv") %>%
  column_to_rownames('X') %>%
  as.matrix()
stopifnot(class(design_matrix) == "matrix")

# Make sure only values present in the dm are "treatment", "control", and "none"
if (unique((design_matrix %>% sapply(., unique) %>% unique) %in% c("treatment", "control", "none")) == TRUE) {
  message("The values of the Design Matrix are correct. You may proceed!")
} else {
  message("The values of the Design Matrix are NOT correct. Fix them!")
}

# Order counts sample names as in the design matrix

counts <- counts[, rownames(design_matrix)]

# Fix gene_lengths transcripts as in the counts dataset

gene_lengths <- gene_lengths[-which(!names(gene_lengths) %in% rownames(counts))]

# Create data frame for each Contrast

# dm_1 = Inf+ vs. Inf-
# dm_2 = Inf+ vs. SIRS-
# dm_3 = Sepsis vs. Inf-
# dm_4 = Sepsis vs. SIRS-
# dm_5 = Sepsis vs. SIRS+
# dm_6 = Sepsis vs. UInf+
# dm_7 = Sepsis+ vs. Sepsis-

dm_1 <- design_matrix %>% as.data.frame() %>% select(Infection_vs_Control) %>% as.matrix()
dm_2 <- design_matrix %>% as.data.frame() %>% select(Infection_vs_Control_ODneg) %>% as.matrix()
dm_3 <- design_matrix %>% as.data.frame() %>% select(Infection_ODpos_vs_Control) %>% as.matrix()
dm_4 <- design_matrix %>% as.data.frame() %>% select(Infection_ODpos_vs_Control_ODneg) %>% as.matrix()
dm_5 <- design_matrix %>% as.data.frame() %>% select(Infection_ODpos_vs_Control_ODpos) %>% as.matrix()
dm_6 <- design_matrix %>% as.data.frame() %>% select(Infection_ODpos_vs_Infection_ODneg) %>% as.matrix()
dm_7 <- design_matrix %>% as.data.frame() %>% select(Infection_ODpos_vs_Infection_ODneg_and_Control) %>% as.matrix()

# Sanity check

stopifnot(colnames(dm_1) == "Infection_vs_Control")
stopifnot(colnames(dm_2) == "Infection_vs_Control_ODneg")
stopifnot(colnames(dm_3) == "Infection_ODpos_vs_Control")
stopifnot(colnames(dm_4) == "Infection_ODpos_vs_Control_ODneg")
stopifnot(colnames(dm_5) == "Infection_ODpos_vs_Control_ODpos")
stopifnot(colnames(dm_6) == "Infection_ODpos_vs_Infection_ODneg")
stopifnot(colnames(dm_7) == "Infection_ODpos_vs_Infection_ODneg_and_Control")


# Clustering ###################################################################

# Create conditions for each design matrix created before

conds1 <- conditions_from_design(dm_1)
conds2 <- conditions_from_design(dm_2)
conds3 <- conditions_from_design(dm_3)
conds4 <- conditions_from_design(dm_4)
conds5 <- conditions_from_design(dm_5)
conds6 <- conditions_from_design(dm_6)
conds7 <- conditions_from_design(dm_7)

# Hierarchical clustering heat map version

make_heat_clustering_plot(
    file.path(plotDir, "heat_hierarchical_clustering_1"),
    counts    = counts[, conds1 != "none"],
    conds     = conds1[conds1 != "none"]
)
make_heat_clustering_plot(
  file.path(plotDir, "heat_hierarchical_clustering_2"),
  counts    = counts[, conds2 != "none"],
  conds     = conds2[conds2 != "none"]
)
make_heat_clustering_plot(
  file.path(plotDir, "heat_hierarchical_clustering_3"),
  counts    = counts[, conds3 != "none"],
  conds     = conds3[conds3 != "none"]
)
make_heat_clustering_plot(
  file.path(plotDir, "heat_hierarchical_clustering_4"),
  counts    = counts[, conds4 != "none"],
  conds     = conds4[conds4 != "none"]
)
make_heat_clustering_plot(
  file.path(plotDir, "heat_hierarchical_clustering_5"),
  counts    = counts[, conds5 != "none"],
  conds     = conds5[conds5 != "none"]
)
make_heat_clustering_plot(
  file.path(plotDir, "heat_hierarchical_clustering_6"),
  counts    = counts[, conds6 != "none"],
  conds     = conds6[conds6 != "none"]
)
make_heat_clustering_plot(
  file.path(plotDir, "heat_hierarchical_clustering_7"),
  counts    = counts[, conds7 != "none"],
  conds     = conds7[conds7 != "none"]
)

# Hierarchical clustering

make_hclust_plot(
  file.path(plotDir, "hierarchical_clustering_1"),
  counts    = counts[, conds1 != "none"],
  conds     = conds1[conds1 != "none"]
)
make_hclust_plot(
  file.path(plotDir, "hierarchical_clustering_2"),
  counts    = counts[, conds2 != "none"],
  conds     = conds2[conds2 != "none"]
)
make_hclust_plot(
  file.path(plotDir, "hierarchical_clustering_3"),
  counts    = counts[, conds3 != "none"],
  conds     = conds3[conds3 != "none"]
)
make_hclust_plot(
  file.path(plotDir, "hierarchical_clustering_4"),
  counts    = counts[, conds4 != "none"],
  conds     = conds4[conds4 != "none"]
)
make_hclust_plot(
  file.path(plotDir, "hierarchical_clustering_5"),
  counts    = counts[, conds5 != "none"],
  conds     = conds5[conds5 != "none"]
)
make_hclust_plot(
  file.path(plotDir, "hierarchical_clustering_6"),
  counts    = counts[, conds6 != "none"],
  conds     = conds6[conds6 != "none"]
)
make_hclust_plot(
  file.path(plotDir, "hierarchical_clustering_7"),
  counts    = counts[, conds7 != "none"],
  conds     = conds7[conds7 != "none"]
)


# Pearson correlation ##########################################################

# Note: Regression is performed. This takes a considerable amount of time
# for large gene sets and/or samples.

make_correlation_plots(
  dat          = mrn,
  outDir       = plotDir,
  prefix       = "corr_1",
  designMatrix = dm_1
)
make_correlation_plots(
  dat          = mrn,
  outDir       = plotDir,
  prefix       = "corr_2",
  designMatrix = dm_2
)
make_correlation_plots(
  dat          = mrn,
  outDir       = plotDir,
  prefix       = "corr_3",
  designMatrix = dm_3
)
make_correlation_plots(
  dat          = mrn,
  outDir       = plotDir,
  prefix       = "corr_4",
  designMatrix = dm_4
)
make_correlation_plots(
  dat          = mrn,
  outDir       = plotDir,
  prefix       = "corr_5",
  designMatrix = dm_5
)
make_correlation_plots(
  dat          = mrn,
  outDir       = plotDir,
  prefix       = "corr_6",
  designMatrix = dm_6
)
make_correlation_plots(
  dat          = mrn,
  outDir       = plotDir,
  prefix       = "corr_7",
  designMatrix = dm_7
)


# PCA plots ####################################################################

# Note: if designMatrix is supplied, 'conds' is ignored. In that case, use
# 'designMatrix = NA'

shapes <- NULL            # <-- optional. fill in!
shapeName <- "Sample"     # <-- optional. fill in!
colName <- "Contrasts"    # <-- fill in!

# Inf+ vs. Inf-

unique(conds1)

make_PCA_plot(
  file         = file.path(plotDir, "pca_Infection_vs_Control_MRN"),
  counts       = counts[, conds1 != "none"],
  conds        = conds1[conds1 != "none"],
  main         = "PCA Infection_vs_Control",
  colName      = colName,
  shapeName    = shapeName,
  norm         = "mrn",
  overwrite    = TRUE,
  add_eclipse  = TRUE
)

# Inf+ vs. SIRS-

unique(conds2)

make_PCA_plot(
  file         = file.path(plotDir, "pca_Infection_vs_Control_ODneg_MRN"),
  counts       = counts[, conds2 != "none"],
  conds        = conds2[conds2 != "none"],
  main         = "PCA Infection_vs_Control_ODneg",
  colName      = colName,
  shapeName    = shapeName,
  norm         = "mrn",
  overwrite    = TRUE,
  add_eclipse  = TRUE
)

# Sepsis vs. Inf-

unique(conds3)

make_PCA_plot(
  file         = file.path(plotDir, "pca_Infection_ODpos_vs_Control_MRN"),
  counts       = counts[, conds3 != "none"],
  conds        = conds3[conds3 != "none"],
  main         = "PCA Infection_ODpos_vs_Control",
  colName      = colName,
  shapeName    = shapeName,
  norm         = "mrn",
  overwrite    = TRUE,
  add_eclipse  = TRUE
)

# Sepsis vs. SIRS-

unique(conds4)

make_PCA_plot(
  file         = file.path(plotDir, "pca_Infection_ODpos_vs_Control_ODneg_MRN"),
  counts       = counts[, conds4 != "none"],
  conds        = conds4[conds4 != "none"],
  main         = "PCA Infection_ODpos_vs_Control_ODneg",
  colName      = colName,
  shapeName    = shapeName,
  norm         = "mrn",
  overwrite    = TRUE,
  add_eclipse  = TRUE
)

# Sepsis vs. SIRS+

unique(conds5)

make_PCA_plot(
  file         = file.path(plotDir, "pca_Infection_ODpos_vs_Control_ODpos_MRN"),
  counts       = counts[, conds5 != "none"],
  conds        = conds5[conds5 != "none"],
  main         = "PCA Infection_ODpos_vs_Control_ODpos",
  colName      = colName,
  shapeName    = shapeName,
  norm         = "mrn",
  overwrite    = TRUE,
  add_eclipse  = TRUE
)

# Sepsis vs. UInf+

unique(conds6)

make_PCA_plot(
  file         = file.path(plotDir, "pca_Infection_ODpos_vs_Infection_ODneg_MRN"),
  counts       = counts[, conds6 != "none"],
  conds        = conds6[conds6 != "none"],
  main         = "PCA Infection_ODpos_vs_Infection_ODneg",
  colName      = colName,
  shapeName    = shapeName,
  norm         = "mrn",
  overwrite    = TRUE,
  add_eclipse  = TRUE
)

# Sepsis+ vs. Sepsis-

unique(conds7)

make_PCA_plot(
  file         = file.path(plotDir, "pca_Infection_ODpos_vs_Infection_ODneg_and_Control_MRN"),
  counts       = counts[, conds7 != "none"],
  conds        = conds7[conds7 != "none"],
  main         = "PCA Infection_ODpos_vs_Infection_ODneg_and_Control",
  colName      = colName,
  shapeName    = shapeName,
  norm         = "mrn",
  overwrite    = TRUE,
  add_eclipse  = TRUE
)

# Save workspace ###############################################################

save.image("./after_pcas_counting_done.RData")
