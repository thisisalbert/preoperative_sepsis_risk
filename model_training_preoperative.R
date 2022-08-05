# Title: Machine Learning Training with selected features
# Project: Risk assessment with gene expression markers in sepsis development
# Mantainer: Albert Garcia-Lopez
# Last update: June 23, 2022

outDir <- "" # <-- fill in
stopifnot(dir.exists(outDir))
pacman::p_load(
  tidyverse, devtools, caret, ROSE, parallel, doParallel, DescTools, paletteer
)
options(readr.show_col_types = FALSE)


# Aesthetics ###################################################################

font_size <- 25
col_pal <- paletteer_d("ggthemes::Tableau_10")
name_contrasts <- c(
  "Inf+ vs. Inf-", "Inf+ vs. SIRS-", "Inf+ vs. SIRS+", "Sepsis+ vs. Inf-",
  "Sepsis+ vs. SIRS-", "Sepsis+ vs. SIRS+", "Sepsis+ vs. UInf+",
  "Sepsis+ vs. Sepsis-"
)
y_scale <- seq.int(0.4, 0, length.out = 8)


# Load counts reads ############################################################

load("after_pcas_counting_done.RData") # # output from pre-processing pipeline
rm(list = setdiff(ls(), c("counts", "outDir")))
stopifnot(ncol(counts) == 267)


# Load custom functions ########################################################

source_url("https://github.com/thisisalbert/machine_learning/raw/main/customSummary.R")
source_url("https://github.com/thisisalbert/machine_learning/raw/main/makeCVseeds.R")
source_url("https://github.com/thisisalbert/machine_learning/raw/main/getCM.R")
source_url("https://github.com/thisisalbert/machine_learning/raw/main/plotCM.R")
source_url("https://github.com/thisisalbert/bioinformatics/raw/main/saveWorkspace.R")


# Load RFE results and get the best DEGs #######################################

getDEGs <- function(rds_path) {
  read_rds(rds_path) %>%
    map(~ .x %>% chuck("coefnames")) %>%
    unlist() %>%
    unique()
}

degs_infpos_infneg <- getDEGs("01_Infection+_vs_Infection-/models_iteration.rds")
degs_infpos_sirsneg <- getDEGs("03_Infection+_vs_SIRS-/models_iteration.rds")
degs_infpos_sirspos <- getDEGs("02_Infection+_vs_SIRS+/models_iteration.rds")
degs_sepsispos_infneg <- getDEGs("06_Sepsis_vs_Infection-/models_iteration.rds")
degs_sepsispos_sirsneg <- getDEGs("08_Sepsis_vs_SIRS-/models_iteration.rds")
degs_sepsispos_sirspos <- getDEGs("07_Sepsis_vs_SIRS+/models_iteration.rds")
degs_sepsispos_uinfpos <- getDEGs("09_Sepsis_vs_Uncomplicated_Infection/models_iteration.rds")
degs_sepsispos_sepsisneg <- getDEGs("05_Sepsis_vs_All/models_iteration.rds")


# Load metadata files ##########################################################

meta_infpos_infneg <- read_tsv("infpos_infneg.tsv")
meta_infpos_sirspos <- read_tsv("infpos_sirspos.tsv")
meta_infpos_sirsneg <- read_tsv("infpos_sirsneg.tsv")
meta_sepsispos_infneg <- read_tsv("sepsispos_infneg.tsv")
meta_sepsispos_sirspos <- read_tsv("sepsispos_sirspos.tsv")
meta_sepsispos_sirsneg <- read_tsv("sepsispos_sirsneg.tsv")
meta_sepsispos_uinfpos <- read_tsv("sepsispos_uinfpos.tsv")
meta_sepsispos_sepsisneg <- read_tsv("sepsispos_sepsisneg.tsv")


# Sequence of seeds ############################################################

seq_seeds <- sample.int(n = 2^15, size = 2^15, replace = FALSE)


# Model training ###############################################################

MLProcedure <- function(
    degs,
    meta,
    comparison,
    n_iter = 100,
    n_cores = 20,
    tunes = 1000,
    metric_to_optimize = "Bal_Accuracy",
    cm_optimization = "misclasserror"
) {

  res_iteration <- list()
  i <- 1

  while (i <= n_iter) {

    # Parallelization

    registerDoParallel(cores = n_cores)

    # Start iteration

    message(paste0("Model ", i))
    fix_seed <- sample(seq_seeds, 1)
    seq_seeds <- seq_seeds[seq_seeds != fix_seed]

    # Selected counts for the DEGs to be used as predictors

    counts_degs <- counts %>%
      as.data.frame() %>%
      rownames_to_column("gene_id") %>%
      as_tibble() %>%
      filter(gene_id %in% degs) %>%
      column_to_rownames("gene_id") %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("ID") %>%
      as_tibble()

    stopifnot(all(degs %in% colnames(counts_degs)) == TRUE)

    # Define case and control

    case <- meta$Comparison %>%
      unique() %>%
      str_subset("^inf_pos$|^sepsis_pos$")

    control <- meta$Comparison %>%
      unique() %>%
      str_subset("^inf_pos$|^sepsis_pos$", negate = TRUE)

    stopifnot(all(unique(meta$Comparison) %in% c(case, control)) == TRUE)

    # Balance genders

    num_females <- meta %>%
      filter(Gender == "Female") %>%
      nrow()

    set.seed(fix_seed)
    balanced_metadata <- meta %>%
      group_by(Gender) %>%
      slice_sample(n = num_females, replace = FALSE) %>%
      ungroup()

    stopifnot(all(
      duplicated(balanced_metadata$ID)) == FALSE
    )
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
      as.data.frame() %>%
      mutate(Comparison = factor(Comparison, c(case, control)))

    # Training

    set.seed(fix_seed)
    model_train <- tryCatch(suppressWarnings(
      caret::train(
        x = mix_counts_meta %>% select(all_of(degs)),
        y = mix_counts_meta %>% chuck("Comparison"),
        maximize = TRUE,
        method = "svmLinear",
        metric = metric_to_optimize,
        preProcess = c("center", "scale", "nzv"),
        tuneLength = tunes,
        trControl = trainControl(
          method = "cv",
          number =  10,
          search = "random",
          sampling = "rose",
          returnData = TRUE,
          returnResamp = "final",
          savePredictions = "final",
          classProbs = TRUE,
          summaryFunction = customSummary,
          selectionFunction = "best",
          allowParallel = TRUE,
          seeds = makeCVseeds(type = "tuneLength", tunes = tunes, numCV = 10)
        )
      )), error = function(e) {e}
    )

    if (inherits(model_train, "error")) {next}

    # Variable Importance

    set.seed(fix_seed)
    var_imp <- varImp(model_train) %>%
      chuck("importance") %>%
      rownames_to_column("Ensembl_Gene_ID") %>%
      mutate(Importance = rowMeans(.[,-1])) %>%
      mutate(Levels = paste(case, control, sep = ",")) %>%
      arrange(desc(Importance)) %>%
      select(Ensembl_Gene_ID, Levels, Importance)

    # Confusion Matrix

    set.seed(fix_seed)
    cm_obj <- tryCatch(
      getCM(
        model = model_train,
        prob_preds = model_train$pred,
        real_labels = model_train$pred$obs,
        case = case,
        control = control,
        type_opt = cm_optimization
      ), error = function(e) {e}
    )

    if (inherits(cm_obj, "error")) {next}

    # AUROC

    set.seed(fix_seed)
    roc_df <- tryCatch(
      PRROC::roc.curve(
        scores.class0 = model_train$pred[[case]],
        weights.class0 = ifelse(model_train$pred$obs == case, 1, 0),
        curve = TRUE,
        max.compute = TRUE,
        min.compute = TRUE
      ), error = function(e) {e}
    )

    if (inherits(roc_df, "error")) {next}

    roc_curve <- roc_df %>%
      chuck("curve") %>%
      as.data.frame() %>%
      rename(FPR = V1, TPR = V2, Threshold = V3) %>%
      as_tibble()

    # Compile results

    res_df <- model_train %>%
      chuck("resample") %>%
      summarise(across(Accuracy:F1, ~ mean(.x, na.rm = TRUE))) %>%
      mutate(
        Model = i,
        Seed = fix_seed,
        Levels = paste(case, control, sep = ","),
        Selected_DEGs_n = length(colnames(model_train$finalModel@xmatrix[[1]])),
        Selected_DEGs = paste(colnames(model_train$finalModel@xmatrix[[1]]), collapse = ","),
      ) %>%
      select(Model, Seed, Levels, everything()) %>%
      as_tibble()

    # Export to list

    res_iteration[[i]] <- tibble(
      Comparison = comparison,
      model = i,
      fix_seed = fix_seed,
      degs = list(degs),
      model_train = list(model_train),
      mix_counts_meta = list(mix_counts_meta),
      res_df = list(res_df),
      var_imp = list(var_imp),
      cm_obj = list(cm_obj),
      roc_curve = list(roc_curve)
    )

    i <- i + 1

  }

  return(res_iteration)

}


# Apply ML custom function #####################################################

res_infpos_infneg <- MLProcedure(
  degs = degs_infpos_infneg,
  meta = meta_infpos_infneg,
  comparison = "Inf+ vs. Inf-"
)

res_infpos_sirspos <- MLProcedure(
  degs = degs_infpos_sirspos,
  meta = meta_infpos_sirspos,
  comparison = "Inf+ vs. SIRS+"
)

res_infpos_sirsneg <- MLProcedure(
  degs = degs_infpos_sirsneg,
  meta = meta_infpos_sirsneg,
  comparison = "Inf+ vs. SIRS-"
)

res_sepsispos_infneg <- MLProcedure(
  degs = degs_sepsispos_infneg,
  meta = meta_sepsispos_infneg,
  comparison = "Sepsis+ vs. Inf-"
)

res_sepsispos_sirspos <- MLProcedure(
  degs = degs_sepsispos_sirspos,
  meta = meta_sepsispos_sirspos,
  comparison = "Sepsis+ vs. SIRS+"
)

res_sepsispos_sirsneg <- MLProcedure(
  degs = degs_sepsispos_sirsneg,
  meta = meta_sepsispos_sirsneg,
  comparison = "Sepsis+ vs. SIRS-"
)

res_sepsispos_uinfpos <- MLProcedure(
  degs = degs_sepsispos_uinfpos,
  meta = meta_sepsispos_uinfpos,
  comparison = "Sepsis+ vs. UInf+"
)

res_sepsispos_sepsisneg <- MLProcedure(
  degs = degs_sepsispos_sepsisneg,
  meta = meta_sepsispos_sepsisneg,
  comparison = "Sepsis+ vs. Sepsis-"
)

saveWorkspace()


# Export objects ###############################################################

write_rds(res_infpos_infneg, paste0(outDir, "res_infpos_infneg.rds"))
write_rds(res_infpos_sirsneg, paste0(outDir, "res_infpos_sirsneg.rds"))
write_rds(res_infpos_sirspos, paste0(outDir, "res_infpos_sirspos.rds"))
write_rds(res_sepsispos_infneg, paste0(outDir, "res_sepsispos_infneg.rds"))
write_rds(res_sepsispos_sirsneg, paste0(outDir, "res_sepsispos_sirsneg.rds"))
write_rds(res_sepsispos_sirspos, paste0(outDir, "res_sepsispos_sirspos.rds"))
write_rds(res_sepsispos_uinfpos, paste0(outDir, "res_sepsispos_uinfpos.rds"))
write_rds(res_sepsispos_sepsisneg, paste0(outDir, "res_sepsispos_sepsisneg.rds"))


# Mean metrics #################################################################

getMeanMetrics <- function(res_list, comparison) {
  res_list %>%
    map_dfr(~ .x %>% chuck("res_df")) %>%
    mutate(across(Accuracy:F1, ~ mean(.x, na.rm = TRUE))) %>%
    select(Levels:ncol(.)) %>%
    distinct() %>%
    mutate(Comparison = comparison, .before = everything())
}

res_mean <- bind_rows(
  getMeanMetrics(res_infpos_infneg, "Inf+ vs. Inf-"),
  getMeanMetrics(res_infpos_sirspos, "Inf+ vs. SIRS+"),
  getMeanMetrics(res_infpos_sirsneg, "Inf+ vs. SIRS-"),
  getMeanMetrics(res_sepsispos_infneg, "Sepsis+ vs. Inf-"),
  getMeanMetrics(res_sepsispos_sirspos, "Sepsis+ vs. SIRS+"),
  getMeanMetrics(res_sepsispos_sirsneg, "Sepsis+ vs. SIRS-"),
  getMeanMetrics(res_sepsispos_uinfpos, "Sepsis+ vs. UInf+"),
  getMeanMetrics(res_sepsispos_sepsisneg, "Sepsis+ vs. Sepsis-")
) %>%
  mutate(Comparison = factor(Comparison, c(
    "Inf+ vs. Inf-", "Inf+ vs. SIRS-", "Inf+ vs. SIRS+", "Sepsis+ vs. Inf-",
    "Sepsis+ vs. SIRS-", "Sepsis+ vs. SIRS+", "Sepsis+ vs. UInf+",
    "Sepsis+ vs. Sepsis-"
  ))) %>%
  arrange(Comparison)

write_tsv(
  x = res_mean,
  file = paste0(outDir, "res_mean.tsv")
)


# Select model representing mean metrics #######################################

getMeanModel <- function(res_list, comparison) {

  mean_values <- res_list %>%
    map_dfr(~ .x %>% chuck("res_df")) %>%
    summarise(across(Accuracy:F1, ~ mean(.x, na.rm = TRUE)))

  mean_model_id <- res_list %>%
    map_dfr(~ .x %>% chuck("res_df")) %>%
    filter(ROC == Closest(x = ROC, a = mean_values$ROC)) %>%
    filter(Bal_Accuracy == Closest(x = Bal_Accuracy, a = mean_values$Bal_Accuracy)) %>%
    filter(Accuracy == Closest(x = Accuracy, a = mean_values$Accuracy)) %>%
    pull(Model)

  res_list[[mean_model_id]] %>%
    mutate(Comparison = comparison, .before = everything())

}

selected_models <- bind_rows(
  getMeanModel(res_infpos_infneg, "Inf+ vs. Inf-"),
  getMeanModel(res_infpos_sirsneg, "Inf+ vs. SIRS-"),
  getMeanModel(res_infpos_sirspos, "Inf+ vs. SIRS+"),
  getMeanModel(res_sepsispos_infneg, "Sepsis+ vs. Inf-"),
  getMeanModel(res_sepsispos_sirsneg, "Sepsis+ vs. SIRS-"),
  getMeanModel(res_sepsispos_sirspos, "Sepsis+ vs. SIRS+"),
  getMeanModel(res_sepsispos_uinfpos, "Sepsis+ vs. UInf+"),
  getMeanModel(res_sepsispos_sepsisneg, "Sepsis+ vs. Sepsis-")
)

write_rds(
  x = selected_models,
  file = paste0(outDir, "selected_models.rds")
)


# Get CM (Training) ############################################################

calculateCM <- function(comparison, color) {

  model = selected_models %>%
    filter(Comparison == comparison) %>%
    pull(model_train) %>%
    chuck(1)

  case = model$levels %>% str_subset("^inf_pos$|^sepsis_pos$")
  control = model$levels %>% str_subset("^inf_pos$|^sepsis_pos$", negate = TRUE)

  cm_obj <- getCM(
    model = model,
    prob_preds = model$pred,
    real_labels = model$pred$obs,
    case = case,
    control = control,
    type_opt = "misclasserror"
  )

  new_case <- ifelse(case == "inf_pos", "Inf+", "Sepsis+")
  new_control <- case_when(
    control == "inf_neg" ~ "Inf-",
    control == "sirs_neg" ~ "SIRS-",
    control == "sirs_pos" ~ "SIRS+",
    control == "uinf_pos" ~ "UInf+",
    control == "sepsis_neg" ~ "Sepsis-"
  )

  n_case <- cm_obj$table %>%
    as.data.frame() %>%
    filter(Reference == case) %>%
    pull(Freq) %>%
    sum() %>%
    paste0(new_case, "\n(n = ", ., ")")

  n_control <- cm_obj$table %>%
    as.data.frame() %>%
    filter(Reference == control) %>%
    pull(Freq) %>%
    sum() %>%
    paste0(new_control, "\n(n = ", ., ")")

  cm_plot <- cm_obj$table %>%
    as.data.frame() %>%
    mutate(Status = ifelse(Reference == Prediction, "Hit", "Miss")) %>%
    mutate(Prediction = ifelse(Prediction == case, new_case, new_control)) %>%
    mutate(Prediction = factor(Prediction, c(new_case, new_control))) %>%
    mutate(Reference = ifelse(Reference == case, n_case, n_control)) %>%
    mutate(Reference = factor(Reference, c(n_case, n_control))) %>%
    plotCM(
      x = "Reference",
      y = "Prediction",
      label = "Freq",
      fill = "Status",
      font_size = font_size,
      color = color
    )

  return(list("cm_obj" = cm_obj, "cm_plot" = cm_plot))

}

cm_infpos_infneg <- calculateCM("Inf+ vs. Inf-", col_pal[1])
cm_infpos_sirsneg <- calculateCM("Inf+ vs. SIRS-", col_pal[2])
cm_infpos_sirspos <- calculateCM("Inf+ vs. SIRS+", col_pal[3])
cm_sepsispos_infneg <- calculateCM("Sepsis+ vs. Inf-", col_pal[4])
cm_sepsispos_sirsneg <- calculateCM("Sepsis+ vs. SIRS-", col_pal[5])
cm_sepsispos_sirspos <- calculateCM("Sepsis+ vs. SIRS+", col_pal[6])
cm_sepsispos_uinfpos <- calculateCM("Sepsis+ vs. UInf+", col_pal[7])
cm_sepsispos_sepsisneg <- calculateCM("Sepsis+ vs. Sepsis-", col_pal[8])

pdf(
  file = paste0(outDir, format(Sys.time(), "%y_%m_%d_cm_training.pdf")),
  width = 21,
  height = 10,
  onefile = TRUE
)
cowplot::plot_grid(
  plotlist = list(
    cm_infpos_infneg$cm_plot,
    cm_infpos_sirsneg$cm_plot,
    cm_infpos_sirspos$cm_plot,
    cm_sepsispos_infneg$cm_plot
  ),
  align = "hv",
  axis = "tblr"
)
cowplot::plot_grid(
  plotlist = list(
    cm_sepsispos_sirsneg$cm_plot,
    cm_sepsispos_sirspos$cm_plot,
    cm_sepsispos_uinfpos$cm_plot,
    cm_sepsispos_sepsisneg$cm_plot
  ),
  align = "hv",
  axis = "tblr"
)
dev.off()


# Plot ROC curves ##############################################################

mean_aurocs <- res_mean %>%
  select(Comparison, ROC) %>%
  mutate(ROC = formatC(ROC, 3, format = "f")) %>%
  deframe()

pdf(
  file = paste0(outDir, format(Sys.time(), "%y_%m_%d_roc_plot_training.pdf")),
  width = 12,
  height = 12
)
map_dfr(
  name_contrasts,
  ~ selected_models %>%
    filter(Comparison == .x) %>%
    pull(roc_curve) %>%
    chuck(1) %>%
    mutate(Comparison = .x, .before = everything())
) %>%
  mutate(Comparison = factor(Comparison, name_contrasts)) %>%
  ggplot(aes(x = FPR, y = TPR, color = Comparison)) +
  geom_path(size = 1.5) +
  coord_fixed() +
  scale_color_paletteer_d("ggthemes::Tableau_10") +
  geom_abline(slope = 1, intercept = 0, color = "darkgrey") +
  annotate(
    geom = "text", x = 1, y = y_scale[1], color = col_pal[1],
    size = 14, fontface = "bold", hjust = 1,
    label = paste0("Inf+ vs. Inf- = ", mean_aurocs[["Inf+ vs. Inf-"]])
  ) +
  annotate(
    geom = "text", x = 1, y = y_scale[2], color = col_pal[2],
    size = 14, fontface = "bold", hjust = 1,
    label = paste0("Inf+ vs. SIRS- = ", mean_aurocs[["Inf+ vs. SIRS-"]])
  ) +
  annotate(
    geom = "text", x = 1, y = y_scale[3], color = col_pal[3],
    size = 14, fontface = "bold", hjust = 1,
    label = paste0("Inf+ vs. SIRS+ = ", mean_aurocs[["Inf+ vs. SIRS+"]])
  ) +
  annotate(
    geom = "text", x = 1, y = y_scale[4], color = col_pal[4],
    size = 14, fontface = "bold", hjust = 1,
    label = paste0("Sepsis+ vs. Inf- = ", mean_aurocs[["Sepsis+ vs. Inf-"]])
  ) +
  annotate(
    geom = "text", x = 1, y = y_scale[5], color = col_pal[5],
    size = 14, fontface = "bold", hjust = 1,
    label = paste0("Sepsis+ vs. SIRS- = ", mean_aurocs[["Sepsis+ vs. SIRS-"]])
  ) +
  annotate(
    geom = "text", x = 1, y = y_scale[6], color = col_pal[6],
    size = 14, fontface = "bold", hjust = 1,
    label = paste0("Sepsis+ vs. SIRS+ = ", mean_aurocs[["Sepsis+ vs. SIRS+"]])
  ) +
  annotate(
    geom = "text", x = 1, y = y_scale[7], color = col_pal[7],
    size = 14, fontface = "bold", hjust = 1,
    label = paste0("Sepsis+ vs. UInf+ = ", mean_aurocs[["Sepsis+ vs. UInf+"]])
  ) +
  annotate(
    geom = "text", x = 1, y = y_scale[8], color = col_pal[8],
    size = 14, fontface = "bold", hjust = 1,
    label = paste0("Sepsis+ vs. Sepsis- = ", mean_aurocs[["Sepsis+ vs. Sepsis-"]])
  ) +
  theme(
    text = element_text(size = font_size),
    axis.text = element_text(size = font_size),
    axis.title = element_text(face = "bold", size = font_size),
    legend.position = "none",
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 2),
    panel.grid = element_blank()
  )
dev.off()


# Save workspace ###############################################################

saveWorkspace()
