# Title: Machine Learning preoperative model validation on COVID-19 data
# Project: Risk assessment with gene expression markers in sepsis development
# Mantainer: Albert Garcia-Lopez
# Last update: July 6, 2022

outDir <- "" # <-- fill in
dataDir <- "" # <-- fill in
stopifnot(dir.exists(outDir) & dir.exists(dataDir))
pacman::p_load(tidyverse, DescTools, devtools, caret)
options(readr.show_col_types = FALSE)
source_url("https://github.com/thisisalbert/machine_learning/raw/main/plotCM.R")
source_url("https://github.com/thisisalbert/machine_learning/raw/main/saveWorkspace.R")
source_url("https://github.com/thisisalbert/machine_learning/raw/main/getTestResults.R")
font_size <- 25


# Load COVID-19 data ###########################################################

counts_covid <- read_tsv("counts.tsv") %>%
  select(!matches(match = "_06[0-9]"))


# Load selected models #########################################################

selected_models <- read_rds(paste0(dataDir, "selected_models.rds"))


# ML Procedure #################################################################

MLProcedure <- function(comparison) {

  defined_model <- selected_models %>% filter(Comparison == comparison)

  # Select random seed

  fix_seed <- defined_model$fix_seed

  # Select model

  presurgery_model <- defined_model$model_train[[1]]

  # Define DEGs

  degs <- colnames(presurgery_model$finalModel@xmatrix[[1]])

  # Define case and control

  model_levels <- defined_model %>%
    chuck("res_df") %>%
    chuck(1) %>%
    pull(Levels) %>%
    str_split(",") %>%
    unlist()

  case <- model_levels %>%
    str_subset("^inf_pos$|^sepsis_pos$")

  control <- model_levels %>%
    str_subset("^inf_pos$|^sepsis_pos$", negate = TRUE)

  # Adapt counts data with class labels

  meta_counts_covid <- counts_covid %>%
    filter(ID %in% degs) %>%
    column_to_rownames("ID") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Sample_ID") %>%
    as_tibble() %>%
    mutate(Comparison = case, .before = everything()) %>%
    mutate(Comparison = factor(Comparison, c(case, control)))

  # Prediction on unlabeled data

  set.seed(fix_seed)
  preds <- list(
    raws = predict.train(
      object = presurgery_model,
      newdata = meta_counts_covid %>% select(all_of(degs)),
      type = "raw"
    ),
    probs = predict.train(
      object = presurgery_model,
      newdata = meta_counts_covid %>% select(all_of(degs)),
      type = "prob"
    )
  )

  InformationValue::confusionMatrix(
    actuals = ifelse(meta_counts_covid$Comparison == case, 1, 0),
    predictedScores = preds$probs[[case]],
    threshold = InformationValue::optimalCutoff(
      actuals = ifelse(presurgery_model$pred$obs == case, 1, 0),
      predictedScores = presurgery_model$pred[[case]],
      optimiseFor = "misclasserror"
    )
  ) %>%
    rownames_to_column("Predicted") %>%
    pivot_longer(cols = "1", names_to = "Reference", values_to = "Freq") %>%
    mutate(across(c(Predicted, Reference), ~ ifelse(.x == "1", case, control))) %>%
    mutate(Comparison = comparison, .before = everything())

}


# Generate predictions #########################################################

res_all <- list()

for (i in selected_models$Comparison) {
  res_all[[i]] <- map_dfr(
    c("Both", "misclasserror", "Ones", "Zeros"),
    ~ MLProcedure(comparison = i, type_cutoff = .x)
  )
}

# Barplot to visualize differences between optimization cutoff types

pdf(
  file = paste0(outDir, format(Sys.time(), "%y_%m_%d_barplot_tp_misclasserror.pdf")),
  width = 11,
  height = 8
)
map_dfr(
  res_all,
  ~ .x %>% filter(str_detect(Predicted, "^inf_pos$|^sepsis_pos$"))
) %>%
  filter(Type_Cutoff == "misclasserror") %>%
  mutate(Comparison = factor(Comparison, c(
    "Inf+ vs. Inf-", "Inf+ vs. SIRS-", "Inf+ vs. SIRS+", "Sepsis+ vs. Inf-",
    "Sepsis+ vs. SIRS-", "Sepsis+ vs. SIRS+", "Sepsis+ vs. UInf+",
    "Sepsis+ vs. Sepsis-"
  ))) %>%
  ggplot(aes(x = Comparison, y = Freq, fill = Comparison, label = Freq)) +
  geom_col(color = "black") +
  geom_text(nudge_y = 1.5, fontface = "bold", size = 8) +
  paletteer::scale_fill_paletteer_d("ggthemes::Tableau_10") +
  scale_y_continuous(n.breaks = 10) +
  labs(
    x = NULL,
    y = "True Positives",
    fill = NULL,
    title = "True Positives 51 COVID-19",
    subtitle = "Prediction using RNA-seq pre-surgery models optimized for 'misclasserror'"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    legend.direction = "horizontal",
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    panel.border = element_rect(size = 1)
  )
dev.off()


# Plot selected model CMs -------------------------------------------------

generateCMplot <- function(res_df, type_cutoff = "misclasserror", col_id = "firebrick") {

  cm_df <- res_df %>%
    filter(Type_Cutoff == type_cutoff) %>%
    select(Comparison, Type_Cutoff, Reference, Predicted, Freq) %>%
    mutate(across(c(Predicted, Reference), ~ case_when(
      .x == "inf_pos" ~ "Inf+",
      .x == "inf_neg" ~ "Inf-",
      .x == "sirs_pos" ~ "SIRS+",
      .x == "sirs_neg" ~ "SIRS-",
      .x == "sepsis_pos" ~ "Sepsis+",
      .x == "sepsis_neg" ~ "Sepsis-",
      .x == "uinf_pos" ~ "UInf+"
    ))) %>%
    mutate(Status = ifelse(Reference == Predicted, "Hit", "Miss"))

  classes <- cm_df %>%
    pull(Predicted) %>%
    unique()

  case <- classes %>% str_subset("^Inf\\+$|^Sepsis\\+$")
  control <- classes %>% str_subset("^Inf\\+$|^Sepsis\\+$", negate = TRUE)

  if (is_empty(case)) {case <- NA}
  if (is_empty(control)) {control <- NA}

  n_case <- cm_df %>%
    filter(Reference == case) %>%
    pull(Freq) %>%
    sum() %>%
    paste0(case, "\n(n = ", ., ")")

  n_control <- cm_df %>%
    filter(Reference == control) %>%
    pull(Freq) %>%
    sum() %>%
    paste0(control, "\n(n = ", ., ")")

  cm_df %>%
    mutate(Reference = case_when(
      Reference == case ~ n_case,
      Reference == control ~ n_control
    )) %>%
    mutate(Reference = factor(Reference, c(n_case, n_control))) %>%
    mutate(Predicted = factor(Predicted, c(case, control))) %>%
    plotCM(
      x = "Reference",
      y = "Predicted",
      label = "Freq",
      fill = "Status",
      color = col_id,
      font_size = 16,
      tile_font_size = 16,
      title = type_cutoff
    )

}

plots_misclasserror <- map(
  1:length(res_all),
  ~ generateCMplot(
    res_df = res_all[[.x]],
    type_cutoff = "misclasserror",
    col_id = paletteer::paletteer_d("ggthemes::Tableau_10")[.x]
  )
)

plots_both <- map(
  1:length(res_all),
  ~ generateCMplot(
    res_df = res_all[[.x]],
    type_cutoff = "Both",
    col_id = paletteer::paletteer_d("ggthemes::Tableau_10")[.x]
  )
)

plots_ones <- map(
  1:length(res_all),
  ~ generateCMplot(
    res_df = res_all[[.x]],
    type_cutoff = "Ones",
    col_id = paletteer::paletteer_d("ggthemes::Tableau_10")[.x]
  )
)

plots_zeros <- map(
  1:length(res_all),
  ~ generateCMplot(
    res_df = res_all[[.x]],
    type_cutoff = "Zeros",
    col_id = paletteer::paletteer_d("ggthemes::Tableau_10")[.x]
  )
)

pdf(
  file = paste0(outDir, format(Sys.time(), "%y_%m_%d_cm_just_covid19.pdf")),
  width = 21,
  height = 10,
  onefile = TRUE
)
cowplot::plot_grid(
  plotlist = plots_misclasserror,
  align = "hv",
  axis = "tblr",
  nrow = 2,
  ncol = 4
)
cowplot::plot_grid(
  plotlist = plots_both,
  align = "hv",
  axis = "tblr",
  nrow = 2,
  ncol = 4
)
cowplot::plot_grid(
  plotlist = plots_ones,
  align = "hv",
  axis = "tblr",
  nrow = 2,
  ncol = 4
)
cowplot::plot_grid(
  plotlist = plots_zeros,
  align = "hv",
  axis = "tblr",
  nrow = 2,
  ncol = 4
)
dev.off()


# Save workspace ###############################################################

saveWorkspace()
