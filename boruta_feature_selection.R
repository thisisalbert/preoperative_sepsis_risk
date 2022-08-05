# Title: Machine Learning feature selection with Boruta
# Project: Risk assessment with gene expression markers in sepsis development
# Mantainer: Albert Garcia-Lopez
# Last update: September 6, 2020

workdir <- "" # <-- fill in
setwd(workdir)
stopifnot(dir.exists(workdir))


# Load packages ################################################################

pacman::p_load(
  tidyverse, stringr, ggfortify, gplots, data.table, devtools, ggpubr, vegan, 
  caret, gmodels, ROCR, pROC, mlbench, plyr, lattice, ranger, randomForest, 
  glmnet, xgboost, e1071, Boruta, biomaRt, shiny, shinythemes, modelgrid, 
  magrittr, foreach, recipes, purrr, BradleyTerry2, corrplot, DataExplorer, 
  vroom
)


# Load counts reads ############################################################

load("after_pcas_counting_done.RData") # output from pre-processing pipeline
rm(list = setdiff(ls(), "counts"))
stopifnot(ncol(counts) == 267)


# Load DEGs ####################################################################

degs <- read.csv("all_unique_degs_in_all_contrasts.csv")
degs <- degs %>% select(x) %>% unlist() %>% as.character()


# Get count reads from DEGs ####################################################

counts_degs <- counts[degs,]
stopifnot(rownames(counts_degs) == degs)


# Create a random sequence of seeds ############################################

seq_seeds_1 <- round(runif(100, min = 0, max = 2^15), 0)
seq_seeds_2 <- round(runif(100, min = 0, max = 2^15), 0)


# Find the best performing DEGs ################################################

results_iteration_1 <- list()

for (i in seq_seeds_1) {
  
  print(paste0("Starting iteration: ", i))
  
  # Import metadata and fix it 
  
  metadata <- read.csv("metadata.csv") %>%
    select(-X) %>%
    mutate(Type = str_replace_all(string = Type, pattern = "Comparator", replacement = "Noninfective")) %>%
    mutate(Type = str_replace_all(string = Type, pattern = "SIRS", replacement = "Noninfective")) %>%
    mutate(Type = str_replace_all(string = Type, pattern = "Sepsis", replacement = "Infective")) %>%
    mutate(Type = factor(Type)) %>%
    mutate(ID = str_replace_all(string = ID, pattern = "R0+", replacement = "R"))
  
  stopifnot(class(metadata$Type) == "factor")
  
  # Equalize males and females in the dataset
  # Note: use prop = 1 on slice_sample() to randomly sample all rows
  
  set.seed(i)
  males_df <- metadata %>% filter(Gender == "Male") %>% slice_sample(prop = 1)
  females_df <- metadata %>% filter(Gender == "Female") %>% slice_sample(prop = 1)
  males_df <- males_df[1:nrow(females_df),]
  metadata <- rbind(males_df, females_df)
  
  # Merge metadata with transcriptomic data
  
  meta_counts <- counts_degs %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("ID") %>% 
    left_join(metadata, ., by = "ID") %>% 
    select(Type, degs) %>% 
    na.omit() %>% 
    remove_rownames()
  
  stopifnot(degs == colnames(meta_counts)[-1])
  stopifnot(meta_counts %>% is.na() %>% sum == 0) 
  
  # Boruta
  
  set.seed(i)
  message("Starting Boruta")
  boruta_degs <- Boruta(
    formula = Type ~ ., 
    data = meta_counts, 
    doTrace = 0, 
    holdHistory = TRUE, 
    pValue = 0.001, 
    maxRuns = 1000
  )
  boruta_degs <- TentativeRoughFix(boruta_degs)
  boruta_degs_names <- getSelectedAttributes(boruta_degs, withTentative = FALSE)
  if (length(boruta_degs_names) <= 25){
    boruta_vars <- boruta_degs_names
  } else {
    boruta_degs <- boruta_degs$ImpHistory[is.finite(rowSums(boruta_degs$ImpHistory)), ]
    boruta_degs <- boruta_degs %>% as.data.frame() %>% dplyr::select(boruta_degs_names)
    boruta_vars <- colSums(boruta_degs) %>% 
    sort(., decreasing = T) %>% 
    .[1:25] %>% 
    names() %>% 
    na.omit()
  }
  
  if (length(boruta_vars) == 0){
    message("No variables were found significant in the Boruta procedure! Jumping to next iteration")
    next
  } else {
    message(
      paste0("Boruta variables selected (", length(boruta_vars), 
        "). Proceeding to train the model with these variables.")
    )
  }
  
  # Training
  
  message("Starting training")
  
  set.seed(i)
  model_train <- 
    train(
      form = formula(paste("Type ~ ", paste(boruta_vars, collapse = "+"))),
      data = meta_counts,
      maximize = TRUE,
      method = "svmLinear",
      metric = "ROC",
      preProcess = c("center", "scale", "nzv"),
      tuneLength = 25,
      trControl = trainControl(
        method = "cv", 
        number =  10,
        search = "random", 
        returnData = TRUE, 
        returnResamp = "final", 
        savePredictions = "final", 
        classProbs = TRUE, 
        summaryFunction = twoClassSummary, 
        selectionFunction = "best",
        allowParallel = TRUE, 
        seeds = NULL)
    )
  
  # Save result of model on CV and go to next iteration
  
  print(paste0("Finishing iteration: ", i))
  
  results_iteration_1[[i]] <- model_train$resample %>%
    select(ROC) %>%
    unlist() %>%
    mean() %>%
    data.frame(
      Iteration = i, 
      Model = "svmLinear", 
      AUC_CV = ., 
      Boruta_n = length(model_train$coefnames),
      Boruta_DEGs = paste(model_train$coefnames, collapse = ",")
    )
}

results_iteration_1_df <- results_iteration_1 %>% bind_rows()
results_iteration_1_export <- results_iteration_1_df %>% 
  mutate(Iteration = as.character(Iteration)) %>% 
  add_row(Iteration = "Mean", 
          Model = unique(as.character(.$Model)),
          AUC_CV = mean(.$AUC_CV),
          Boruta_n = NULL,
          Boruta_DEGs = NULL
  ) %>% 
  mutate(AUC_CV = round(AUC_CV, 3))
write.csv(results_iteration_1_export, "results_iteration_1_df.csv")


# Select the 10 iterations with the greatest AUC ###############################

selected_degs <- results_iteration_1_df %>% 
  slice_max(order_by = AUC_CV, n = 10) %>% 
  select(Boruta_DEGs)
selected_degs <- paste(selected_degs$Boruta_DEGs, collapse = ",")
selected_degs <- str_split(string = selected_degs, pattern = ",") %>% 
  unlist() %>% 
  as.character() %>% 
  unique()
write.csv(selected_degs, "./best_selected_degs.csv")
message(paste0("Amount of best DEGs from top 10 iterations: ", length(selected_degs)))


# Selected DEGs on Test ########################################################

results_iteration_2 <- list()

for (i in seq_seeds_2) {
  
  print(paste0("Starting iteration: ", i))
  
  # Import metadata and fix it 
  
  metadata <- read.csv("metadata.csv") %>%
    select(-X) %>%
    mutate(Type = str_replace_all(string = Type, pattern = "Comparator", replacement = "Noninfective")) %>%
    mutate(Type = str_replace_all(string = Type, pattern = "SIRS", replacement = "Noninfective")) %>%
    mutate(Type = str_replace_all(string = Type, pattern = "Sepsis", replacement = "Infective")) %>%
    mutate(Type = factor(Type)) %>%
    mutate(ID = str_replace_all(string = ID, pattern = "R0+", replacement = "R"))
  
  stopifnot(class(metadata$Type) == "factor")
  
  # Equalize males and females in the dataset.
  
  set.seed(i)
  males_df <- metadata %>% filter(Gender == "Male") %>% slice_sample(prop = 1)
  females_df <- metadata %>% filter(Gender == "Female") %>% slice_sample(prop = 1)
  males_df <- males_df[1:nrow(females_df),]
  metadata <- rbind(males_df, females_df)
  
  # Merge metadata with transcriptomic data
  
  meta_counts <- counts_degs %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("ID") %>% 
    left_join(metadata, ., by = "ID") %>% 
    select(Type, degs) %>% 
    na.omit() %>% 
    remove_rownames()
  
  stopifnot(degs == colnames(meta_counts)[-1])
  stopifnot(meta_counts %>% is.na() %>% sum == 0) 
  
  # Training
  
  message("Starting training")
  
  set.seed(i)
  model_train <- 
    train(
      form = formula(paste("Type ~ ", paste(selected_degs, collapse = "+"))),
      data = meta_counts,
      maximize = TRUE,
      method = "svmLinear",
      metric = "ROC",
      preProcess = c("center", "scale", "nzv"),
      tuneLength = 25,
      trControl = trainControl(
        method = "cv", 
        number =  10,
        search = "random", 
        returnData = TRUE, 
        returnResamp = "final", 
        savePredictions = "final", 
        classProbs = TRUE, 
        summaryFunction = twoClassSummary, 
        selectionFunction = "best",
        allowParallel = TRUE, 
        seeds = NULL)
    )
  
  # Save result of model on CV and go to next iteration
  
  print(paste0("Finishing iteration: ", i))
  
  results_iteration_2[[i]] <- model_train$resample %>%
    select(ROC) %>%
    unlist() %>%
    mean() %>%
    data.frame(
      Iteration = i, 
      Model = "svmLinear", 
      AUC_CV = ., 
      Selected_DEGs_n = length(model_train$coefnames),
      Selected_DEGs = paste(model_train$coefnames, collapse = ",")
    )
}

results_iteration_2_df <- results_iteration_2 %>% 
  bind_rows() %>% 
  mutate(Iteration = as.character(Iteration)) %>% 
  add_row(Iteration = "Mean", 
          Model = "svmLinear",
          AUC_CV = mean(.$AUC_CV),
          Selected_DEGs_n = length(selected_degs),
          Selected_DEGs = paste(selected_degs, collapse = ",")
          ) %>% 
  mutate(AUC_CV = round(AUC_CV, 3))

write.csv(results_iteration_2_df, "results_iteration_2_df.csv")


# Save workspace ###############################################################

save.image(paste0(Sys.Date(), "_workspace.RData"))