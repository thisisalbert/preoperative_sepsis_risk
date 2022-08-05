# Title: Machine Learning with qPCR data
# Project: Risk assessment with gene expression markers in sepsis development
# Mantainer: Albert Garcia-Lopez
# Last update: June 27, 2022

outDir <- "" # <-- fill in
stopifnot(dir.exists(outDir))
pacman::p_load(
  tidyverse, devtools, caret, ROSE, parallel, doParallel, paletteer, DescTools
)
options(readr.show_col_types = FALSE)
font_size <- 25


# Load custom functions ########################################################

source_url("https://github.com/thisisalbert/machine_learning/raw/main/customSummary.R")
source_url("https://github.com/thisisalbert/machine_learning/raw/main/makeCVseeds.R")
source_url("https://github.com/thisisalbert/machine_learning/raw/main/getCM.R")
source_url("https://github.com/thisisalbert/machine_learning/raw/main/plotCM.R")
source_url("https://github.com/thisisalbert/bioinformatics/raw/main/saveWorkspace.R")


# Load DEGs for each comparison ################################################

degs_symbols <- readxl::read_excel("list_64_genes_annotated.xlsx") %>%
  select(ensembl = Ensembl_Gene_ID, symbol = External_Gene_Name) %>%
  mutate(symbol = case_when(
    symbol == "MT-ND4L" ~ "MT_ND4L",
    symbol == "Z82206.1" ~ "Z82206_1",
    TRUE ~ symbol
  ))

degs_list <- read_rds("selected_models.rds") %>%
  select(Comparison, degs) %>%
  deframe() %>%
  map(
    ~ tibble(ensembl = .) %>%
      left_join(degs_symbols, by = "ensembl") %>%
      filter(!symbol %in% c("LY6G6E", "F2RL3")) %>%
      pull(symbol)
  )


# Load metadata + qPCR data with inputed NAs ###################################

meta_qpcr <- read_tsv("qPCR/Input_NA_Delta_Ct/meta_qpcr_na_inputed.tsv") %>%
  rename(HROB = C17ORF53)


# Gender info ##################################################################

gender_info <- read_csv("metadata.csv") %>%
  select(ID, Gender) %>%
  mutate(ID = str_replace_all(ID, "R0+", "R"))


# Add gender info to qPCR data #################################################

meta_qpcr <- meta_qpcr %>%
  left_join(gender_info, by = "ID") %>%
  relocate(Gender, .after = ID)


# Make metadata for each comparison ############################################

meta_qpcr_sepsispos_infneg <- meta_qpcr %>%
  filter(GROUP %in% c("Infection_ODpos", "Control_ODneg", "Control_ODpos")) %>%
  mutate(Comparison = ifelse(GROUP == "Infection_ODpos", "sepsis_pos", "inf_neg")) %>%
  mutate(Comparison = factor(Comparison, c("sepsis_pos", "inf_neg"))) %>%
  select(-GROUP) %>%
  select(ID, Comparison, everything())

meta_qpcr_sepsispos_sirsneg <- meta_qpcr %>%
  filter(GROUP %in% c("Infection_ODpos", "Control_ODneg")) %>%
  mutate(Comparison = ifelse(GROUP == "Infection_ODpos", "sepsis_pos", "sirs_neg")) %>%
  mutate(Comparison = factor(Comparison, c("sepsis_pos", "sirs_neg"))) %>%
  select(-GROUP) %>%
  select(ID, Comparison, everything())

meta_qpcr_sepsispos_sirspos <- meta_qpcr %>%
  filter(GROUP %in% c("Infection_ODpos", "Control_ODpos")) %>%
  mutate(Comparison = ifelse(GROUP == "Infection_ODpos", "sepsis_pos", "sirs_pos")) %>%
  mutate(Comparison = factor(Comparison, c("sepsis_pos", "sirs_pos"))) %>%
  select(-GROUP) %>%
  select(ID, Comparison, everything())

meta_qpcr_sepsispos_uinfpos <- meta_qpcr %>%
  filter(GROUP %in% c("Infection_ODpos", "Infection_ODneg")) %>%
  mutate(Comparison = ifelse(GROUP == "Infection_ODpos", "sepsis_pos", "uinf_pos")) %>%
  mutate(Comparison = factor(Comparison, c("sepsis_pos", "uinf_pos"))) %>%
  select(-GROUP) %>%
  select(ID, Comparison, everything())


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

    # Filter data for the DEGs to be used as predictors

    meta_filtered <- meta %>%
      select(ID, Gender, Comparison, all_of(degs))

    stopifnot(all(degs %in% colnames(meta)) == TRUE)

    # Define case and control

    case <- meta_filtered$Comparison %>%
      unique() %>%
      str_subset("^inf_pos$|^sepsis_pos$")

    control <- meta_filtered$Comparison %>%
      unique() %>%
      str_subset(case, negate = TRUE)

    stopifnot(all(unique(meta_filtered$Comparison) %in% c(case, control)) == TRUE)

    # Define subsampling method to correct class imbalances (if necessary)

    if (identical(
      length(which(meta_filtered$Comparison == case)),
      length(which(meta_filtered$Comparison == control))
    )) {
      sampling <- NULL
    } else {
      sampling <- "rose"
    }

    # Balance genders (only if necessary)

    if (!identical(
      meta_filtered %>% filter(Gender == "Male") %>% nrow(),
      meta_filtered %>% filter(Gender == "Female") %>% nrow()
    )) {

      message("Correcting gender imbalance")

      num_females <- meta_filtered %>%
        filter(Gender == "Female") %>%
        nrow()

      set.seed(fix_seed)
      meta_filtered <- meta_filtered %>%
        group_by(Gender) %>%
        slice_sample(n = num_females, replace = FALSE) %>%
        ungroup()

      stopifnot(
        all(duplicated(meta_filtered$ID) == FALSE) &
          identical(
            meta_filtered %>% filter(Gender == "Female") %>% nrow(),
            meta_filtered %>% filter(Gender == "Male") %>% nrow()
          )
      )
    }

    # Training

    set.seed(fix_seed)
    model_train <- tryCatch(suppressWarnings(
      caret::train(
        x = meta_filtered %>% select(all_of(degs)),
        y = meta_filtered %>% chuck("Comparison"),
        maximize = TRUE,
        method = "svmLinear",
        metric = metric_to_optimize,
        preProcess = c("center", "scale", "nzv"),
        tuneLength = tunes,
        trControl = trainControl(
          method = "cv",
          number =  10,
          search = "random",
          sampling = sampling,
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
    var_imp <- filterVarImp(
      x = meta_filtered %>% select(all_of(degs)),
      y = meta_filtered %>% chuck("Comparison"),
    ) %>%
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
    roc_curve <- tryCatch(
      PRROC::roc.curve(
        scores.class0 = model_train$pred[[case]],
        weights.class0 = ifelse(model_train$pred$obs == case, 1, 0),
        curve = TRUE,
        max.compute = TRUE,
        min.compute = TRUE
      ) %>%
        chuck("curve") %>%
        as.data.frame() %>%
        rename(FPR = V1, TPR = V2, Threshold = V3) %>%
        as_tibble()
      , error = function(e) {e}
    )

    if (inherits(roc_curve, "error")) {next}

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
      comparison = comparison,
      model = i,
      fix_seed = fix_seed,
      degs = list(degs),
      model_train = list(model_train),
      meta_filtered = list(meta_filtered),
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

res_ml <- list(
  `Sepsis+ vs. Inf-` = MLProcedure(
    degs = degs_list$`Sepsis+ vs. Inf-`,
    meta = meta_qpcr_sepsispos_infneg,
    comparison = "Sepsis+ vs. Inf-"
  ),
  `Sepsis+ vs. SIRS-` = MLProcedure(
    degs = degs_list$`Sepsis+ vs. SIRS-`,
    meta = meta_qpcr_sepsispos_sirsneg,
    comparison = "Sepsis+ vs. SIRS-"
  ),
  `Sepsis+ vs. SIRS+` = MLProcedure(
    degs = degs_list$`Sepsis+ vs. SIRS+`,
    meta = meta_qpcr_sepsispos_sirspos,
    comparison = "Sepsis+ vs. SIRS+"
  ),
  `Sepsis+ vs. UInf+` = MLProcedure(
    degs = degs_list$`Sepsis+ vs. UInf+`,
    meta = meta_qpcr_sepsispos_uinfpos,
    comparison = "Sepsis+ vs. UInf+"
  )
)

saveWorkspace()

write_rds(
  x = res_ml,
  file = paste0(outDir, "res_ml.rds")
)


# Analyze results ##############################################################

res_mean <- map_dfr(
  res_ml,
  ~ .x %>%
    map_dfr(~ .x %>% chuck("res_df")) %>%
    mutate(across(Accuracy:F1, ~ mean(.x, na.rm = TRUE))) %>%
    select(Levels:ncol(.)) %>%
    distinct()
) %>%
  mutate(Comparison = case_when(
    str_detect(Levels, "inf_neg") ~ "Sepsis+ vs. Inf-",
    str_detect(Levels, "sirs_neg") ~ "Sepsis+ vs. SIRS-",
    str_detect(Levels, "sirs_pos") ~ "Sepsis+ vs. SIRS+",
    str_detect(Levels, "uinf_pos") ~ "Sepsis+ vs. UInf+"
  ), .before = everything())

write_tsv(
  x = res_mean,
  file = paste0(outDir, "res_mean.tsv")
)


# Select average model #########################################################

selected_models_id <- map(
  res_ml,
  ~ .x %>%
    bind_rows() %>%
    pull(res_df) %>%
    bind_rows() %>%
    filter(ROC == Closest(ROC, mean(.$ROC))) %>%
    filter(Bal_Accuracy == Closest(Bal_Accuracy, mean(.$Bal_Accuracy))) %>%
    filter(Accuracy == Closest(Accuracy, mean(.$Accuracy))) %>%
    filter(F1 == Closest(F1, mean(.$F1))) %>%
    filter(Kappa == Closest(Kappa, mean(.$Kappa))) %>%
    slice_sample(n = 1) %>%
    pull(Model)
)

selected_models <- list(
  `Sepsis+ vs. Inf-` = res_ml$`Sepsis+ vs. Inf-` %>%
    chuck(selected_models_id$`Sepsis+ vs. Inf-`),
  `Sepsis+ vs. SIRS-` = res_ml$`Sepsis+ vs. SIRS-` %>%
    chuck(selected_models_id$`Sepsis+ vs. SIRS-`),
  `Sepsis+ vs. SIRS+` = res_ml$`Sepsis+ vs. SIRS+` %>%
    chuck(selected_models_id$`Sepsis+ vs. SIRS+`),
  `Sepsis+ vs. UInf+` = res_ml$`Sepsis+ vs. UInf+` %>%
    chuck(selected_models_id$`Sepsis+ vs. UInf+`)
)


# Plot AUROC ###################################################################

plotROC <- function(comparison) {

  selected_models[[comparison]] %>%
    pull(roc_curve) %>%
    chuck(1) %>%
    ggplot(aes(x = FPR, y = TPR, color = Threshold)) +
    geom_path(size = 2, linejoin = "round", lineend = "round") +
    geom_abline(slope = 1, intercept = 0, color = "darkgrey", size = 1) +
    geom_label(aes(
      x = 0.5,
      y = 0.5,
      label = paste0(
        "AUROC = ",
        res_mean %>%
          filter(Comparison == comparison) %>%
          pull(ROC) %>%
          round(3)
      )
    ),
    size = 14,
    label.size = 1,
    color = "black"
    ) +
    coord_fixed() +
    scale_color_gradientn(colors = rainbow(5), limits = c(0, 1)) +
    labs(title = paste0("RT-qPCR: ", comparison)) +
    theme(
      text = element_text(size = font_size, color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.border = element_rect(fill = NA, size = 2, color = "black"),
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold"),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.width = unit(2, "cm")
    )

}

pdf(
  file = paste0(outDir, format(Sys.time(), "%y_%m_%d_auroc_plot.pdf")),
  width = 12,
  height = 12,
  onefile = TRUE
)
plotROC("Sepsis+ vs. Inf-")
plotROC("Sepsis+ vs. SIRS-")
plotROC("Sepsis+ vs. SIRS+")
plotROC("Sepsis+ vs. UInf+")
dev.off()


# Plot CM ######################################################################

getCMplot <- function(comparison, color_id) {

  df = selected_models[[comparison]] %>%
    pull(cm_obj) %>%
    chuck(1) %>%
    chuck("table") %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(Status = ifelse(Reference == Prediction, "Hit", "Miss")) %>%
    mutate(Status = factor(Status, c("Hit", "Miss")))

  case = "sepsis_pos"

  control = df$Reference %>%
    unique() %>%
    str_subset(case, negate = TRUE)

  new_control = case_when(
    control == "inf_neg" ~ "Inf-",
    control == "sirs_neg" ~ "SIRS-",
    control == "sirs_pos" ~ "SIRS+",
    control == "uinf_pos" ~ "UInf+"
  )

  n_case = df %>%
    filter(Reference == case) %>%
    pull(Freq) %>%
    sum() %>%
    paste0("Sepsis+\n(n = ", ., ")")

  n_control = df %>%
    filter(Reference == control) %>%
    pull(Freq) %>%
    sum() %>%
    paste0(new_control, "\n(n = ", ., ")")

  df %>%
    mutate(Prediction = ifelse(Prediction == case, "Sepsis+", new_control)) %>%
    mutate(Prediction = factor(Prediction, c("Sepsis+", new_control))) %>%
    mutate(Reference = ifelse(Reference == case, n_case, n_control)) %>%
    mutate(Reference = factor(Reference, c(n_case, n_control))) %>%
    plotCM(
      df = .,
      x = "Reference",
      y = "Prediction",
      label = "Freq",
      fill = "Status",
      font_size = font_size,
      color = color_id
    )

}

pdf(
  file = paste0(outDir, format(Sys.time(), "%y_%m_%d_cm_plots.pdf")),
  width = 21,
  height = 12
)
cowplot::plot_grid(
  getCMplot(comparison = "Sepsis+ vs. Inf-", color_id = "#76B7B2FF"),
  getCMplot(comparison = "Sepsis+ vs. SIRS-", color_id = "#59A14FFF"),
  getCMplot(comparison = "Sepsis+ vs. SIRS+", color_id = "#EDC948FF"),
  getCMplot(comparison = "Sepsis+ vs. UInf+", color_id = "#B07AA1FF"),
  align = "hv",
  axis = "tblr"
)
dev.off()


# Save workspace ###############################################################

saveWorkspace()
