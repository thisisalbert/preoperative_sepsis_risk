# Title: COVID19 RNA-seq preprocessing (Script based on GEO2RNAseq pipeline version R320)
# Project: Risk assessment with gene expression markers in sepsis development
# Maintainer: Albert Garcia-Lopez
# Last update: December 6, 2021


# Load libraries ------------------------------------------------------------------------------

library("Geo2RNAseq")
library("R.utils")
library("tidyverse")
library("data.table")
library("gtools")
library("stringi")
# R.utils::sourceDirectory("R/") # only if "Geo2RNAseq" is not installed via `install.packages()`


# Set working directory -----------------------------------------------------------------------

workdir <- "" # <-- fill in
setwd(workdir)
stopifnot(dir.exists(getwd()))


# Global settings -----------------------------------------------------------------------------

# Output directories

indexDir <- paste0(getwd(), "/index")
fastqDir <- paste0(getwd(), "/fastq")
qualDir  <- paste0(getwd(), "/quality")
mapDir   <- paste0(getwd(), "/mapping")
countDir <- paste0(getwd(), "/counting")
tabDir   <- paste0(getwd(), "/result_tables")
plotDir  <- paste0(getwd(), "/result_plots")
degDir   <- paste0(getwd(), "/diff_exp_genes")

dirList <- c(
  indexDir, fastqDir, qualDir, mapDir, countDir, tabDir, plotDir, degDir
)

for (i in dirList) {
  dir.create(i)
}

# If TRUE, existing files will be overwritten without questioning.
# If FALSE, most methods will skip step for existing output files.
# In these cases, the method will return 'call[x] = "NOT USED"'

FORCE_OVERWRITE <- FALSE
paired <- TRUE
MAX_CPUS <- 30


# Genome and annotation files -----------------------------------------------------------------

genome <- "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
anno   <- "Homo_sapiens.GRCh38.99.gtf"
stopifnot(file.exists(genome) & file.exists(anno))


# Build index ---------------------------------------------------------------------------------

index <- make_HiSat2_index(
  genomeFile = genome,
  customOut = paste0(indexDir, "/index"),
  cpus = MAX_CPUS
)


# Load raw FASTQ files ------------------------------------------------------------------------

raw_fastq_files <- list.files(
  path = "/sbiarchiv/AGarcia/COVID19_61Samples/data/files/raw_data/",
  full.names = TRUE, 
  pattern = "\\.gz",
  recursive = TRUE
)

if (paired == TRUE) {
  writeLines(paste0(
    "Total number of files = ", length(raw_fastq_files), "\n",
    "Total number of samples = ", length(raw_fastq_files)/2
  ))
} else {
  writeLines(paste0(
    "Total number of files = ", length(raw_fastq_files), "\n",
    "Total number of samples = ", length(raw_fastq_files)
  ))
}


# Quality control 1 (before trimming) ---------------------------------------------------------

qualityRawDir <- file.path(qualDir, "raw")
if (length(dir(qualityRawDir)) > 0) {
  warning(paste0(
    "Directory \"", qualityRawDir, "\" not empty! Overwrite? ", FORCE_OVERWRITE
  ), immediate. = TRUE)
}

if (FORCE_OVERWRITE || length(dir(qualityRawDir)) == 0) {
  writeLines("FastQC - raw ...")
  fq_res <- run_FastQC(
    files = raw_fastq_files, 
    outDir = qualityRawDir, 
    cpus = MAX_CPUS, 
    extend = TRUE
  )
}


# Trimming ------------------------------------------------------------------------------------

trimmedDir <- fastqDir
windowsizetrimming <- 15
qualcuttrimming <- 25
phred <- "-phred33"
leading <- 3
trailing <- 3
minlen <- 30

if (length(list.files(trimmedDir, pattern = "\\.trimo(\\.pe)*\\.fastq$")) > 0) {
  warning(paste0(
    "Directory \"", trimmedDir, "\" not empty! Overwrite? ", FORCE_OVERWRITE
  ), immediate. = TRUE
  )
}

trimming_res <- run_Trimmomatic(
  files = raw_fastq_files,
  is.paired = paired,
  outDir = trimmedDir,
  cpus = MAX_CPUS,
  windowsize = windowsizetrimming,
  qualcut = qualcuttrimming,
  phred = phred,
  leading = leading,
  trailing = trailing,
  minlen = minlen
)

# Number of reads after trimming based on Trimmomatic output

number_raw <- trimming_res$input
number_trimmed <- trimming_res$surviving

if (paired) {
  trimmed_fastq_files <- asPairVector(trimming_res$files)
} else {
  trimmed_fastq_files <- trimming_res$files
}

names(number_raw) <- basename(raw_fastq_files)
names(number_trimmed) <- basename(trimmed_fastq_files)

fastq_files <- trimmed_fastq_files


# SortMeRNA -----------------------------------------------------------------------------------

filterrRNA <- TRUE
sortmeDir <- file.path(fastqDir, "sortmerna")

if (filterrRNA && length(dir(sortmeDir)) > 0) {
  warning(paste0(
    "Directory \"", sortmeDir, "\" not empty! Overwrite? ", FORCE_OVERWRITE
  ), 
  immediate. = TRUE
  )
}

if (filterrRNA) {
  sortmerna_res <- run_SortMeRNA(
    files = fastq_files,
    outDir = fastqDir,
    mode = "fast",
    paired = paired,
    cpus = MAX_CPUS,
    overwrite = TRUE
  )
  non_rrna_files <- sortmerna_res$files
  number_nonrRNA <- sortmerna_res$nonrrna
}

# If sample/file order changed:

stopifnot(unique(order(trimmed_fastq_files) == order(non_rrna_files)) == TRUE)

fastq_files <- if (filterrRNA) {
  non_rrna_files
} else {
  trimmed_fastq_files
}


# Quality control 2 (after trimming) ----------------------------------------------------------

qualitySubDir <- if (filterrRNA) {
  file.path(qualDir, "non_rrna")
} else {
  file.path(qualDir, "trimmed")
}

if (length(dir(qualitySubDir)) > 0) {
  warning(paste0("Directory \"" , qualitySubDir, "\" not empty! overwrite? ", 
    FORCE_OVERWRITE))
}

if (FORCE_OVERWRITE || length(dir(qualitySubDir)) == 0) {
  writeLines("FastQC - trimmed ...")
  fq_res <- run_FastQC(
    files = fastq_files, 
    outDir = qualitySubDir, 
    cpus = MAX_CPUS, 
    extend = FALSE)
}


# Mapping -------------------------------------------------------------------------------------

index_map <- index
anno_map  <- anno
mapper <- "hisat2"
bamDir <- file.path(mapDir, "bamfiles")
samDir <- file.path(mapDir, "samfiles")
convert_sam <- TRUE # Should SAM files be converted to BAM files if BAM files are not found?

if (mapper == "tophat2") {
  mapping_res <- run_Tophat(
    files = fastq_files,
    index = index_map,
    outDir = mapDir,
    is.paired = paired,
    anno = anno_map,
    addArgs = NA,      # "-g 2 --b2-very-sensitive --no-coverage-search",
    cpus = MAX_CPUS,
    worker = if (paired) 5 else 10,
    use.existing = TRUE
  )
  bam_files <- mapping_res$files
} else {
  if (FORCE_OVERWRITE || length(dir(samDir)) == 0) {
    mapping_res <- run_Hisat2(
      files = fastq_files,
      index = index_map,
      outDir = mapDir,
      is.paired = paired,
      addArgs = "",
      splice = TRUE,
      cpus = MAX_CPUS
    )
    ## NOTE: use that only if the run_Hisat2() parameter 'as.bam' is FALSE.
    # bam_files <- sam_to_bam(mapping_res$mapFiles, sort=TRUE, bamDir = bamDir)
    bam_files <- mapping_res$files
  } else {
    if (convert_sam) {
      sam_files <- list.files(
        samDir, 
        pattern = "\\.sam$", 
        full.names = TRUE
      )
      if (length(sam_files) != number_samples)
        stop(paste(
          "Invalid number of SAM files. Expected", 
          number_samples, 
          "but got", 
          length(sam_files)
        ))
      bam_files <- sam_to_bam(
        mapping_res$mapFiles, 
        sort = TRUE, 
        bamDir = bamDir
      )
    } else
      bam_files <- list.files(
        bamDir, 
        pattern = "\\.bam$", 
        full.names = TRUE
      )
    if (length(bam_files) != number_samples)
      warning(paste(
        "Found", 
        length(bam_files), 
        "BAM files, but expected", 
        number_samples, 
        ". Maybe SAMtools was interrupted. Try to convert again? ", 
        convert_sam
      ), 
      immediate. = TRUE
      )
  }
}


# Counting ------------------------------------------------------------------------------------

anno_count <- anno

if (length(dir(countDir)) > 0) {
  warning(paste(
    "Directory", 
    countDir, 
    "not empty! Overwrite?", 
    FORCE_OVERWRITE
  ), immediate. = TRUE
  )
}

if (FORCE_OVERWRITE || length(dir(countDir)) == 0) {
  res_counting <- run_featureCounts(
    files = bam_files,
    annotation = anno_count,
    isGTF = TRUE,
    IDtype = "gene_id",
    featureType = "exon",
    outDir = countDir,
    isPairedEnd = paired,
    allowMultiOverlap = FALSE,
    cpus = MAX_CPUS
  )
} else {
  res_counting <- read.featureCounts.files(countDir)
}

counts <- res_counting$counts
countfile <- res_counting$countFile
sumfile <- res_counting$sumFile


# Change the names of samples names in counts and lib_sizes

colnames(counts) <- gsub("_1.trimo.pe.non_rrna.bam", "", colnames(counts))
info <- res_counting$anno
lib_sizes <- res_counting$summary[,1]
gene_lengths <- info$Length
names(gene_lengths) <- info$GeneID
names(lib_sizes) <- gsub("_1.trimo.pe.non_rrna.bam", "", names(lib_sizes))
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


# Check for rRNA contamination ----------------------------------------------------------------

detect_high_coverage(counts, lib_sizes)

cont_ids <- c(
  "ENSG00000038427", "ENSG00000101596", "ENSG00000105835", "ENSG00000124942",
  "ENSG00000156508", "ENSG00000163220", "ENSG00000166710", "ENSG00000170315",
  "ENSG00000173821", "ENSG00000245532", "ENSG00000251562", "ENSG00000266412",
  "ENSG00000276168", "ENSG00000278771"
)

counts_clean <- counts[!rownames(counts) %in% cont_ids,]

# Perform detection after cleanup

detect_high_coverage(counts_clean, lib_sizes)
counts <- counts_clean

# Remove also those reads from gene_lengths

gene_lengths <- gene_lengths[!names(gene_lengths) %in% cont_ids]


# SAMtools ------------------------------------------------------------------------------------

flagstatDir <- file.path(mapDir, "flagstats")
flag_files <- make_flagstats(bam_files, flagstatDir, MAX_CPUS)


# MultiQC -------------------------------------------------------------------------------------

muliqc_config_yaml <- file.path(
  "/sbidata/AGarcia/RNA-Seq/multiqc_config_bastian.yaml"
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


# Mapping stats -------------------------------------------------------------------------------

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

precise <- TRUE

mapping_stats_df <- calc_mapping_stats(
  bamFiles = bam_files,
  fqFiles = raw_fastq_files,
  anno = anno_map, # annotation used for mapping ~ can be different from counting!
  numReads = number_raw, # number of raw reads - mandatory!
  numTrimmed = number_trimmed, # number of remaining reads after trimming - or NA
  numNonrRNA = number_nonrRNA, # number of remaining reads after SortMeRNA - or NA
  libSizes = lib_sizes,
  paired = paired,
  precise = precise,
  remove.na = FALSE,
  cpus = MAX_CPUS
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
  if (FALSE %in% (asPaired(number_raw) >= lib_sizes))
    stop("Error: number of reads mapping in exons should not exceed original number of reads!")
} else {
  if (FALSE %in% (number_raw >= lib_sizes))
    stop("Error: number of reads mapping in exons should not exceed original number of reads!")
}


# Normalized count values ---------------------------------------------------------------------

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

