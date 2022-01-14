## !/usr/bin/Rscript
#
# Perform inference on the ELISA data from Dopico et al., 2021.
# We use five chains for the MVN and MVT mixture models and look at the inferred
#   * datasets
#   * classes
#   * parameters
# We plot the complete log-likelihood to check within/across chain convergence.
#
# === Setup ====================================================================

library(ggplot2)
library(magrittr)
library(tidyr)
library(dplyr)
library(tibble)
library(BatchMixtureModel)
library(randomForest)
library(kernlab)
library(stringr)

calcSeroprevalence <- function(predictions, time_vec, pos_value = 1) {
  timepoints <- sort(unique(time_vec))
  out_df <- data.frame(Seroprevalence = 0, Week = timepoints)
  for(t in timepoints) {
    .inds <- which(time_vec == t)
    N_t <- length(.inds)
    N_t_pos <- sum(predictions[.inds] == pos_value)
    sero_est <- (N_t_pos / N_t) * 100
    out_df$Seroprevalence[out_df$Week == t] <- sero_est
  }
  out_df
}

set.seed(1)
my_wd <- "./"

# "https://github.com/chr1swallace/seroprevalence-paper/blob/master/adjusted-data.RData"
# Please download before proceeding if it is not already in the Data directory.
# Read in the ELISA data
data_file <- paste0(my_wd, "Data/Sweden/adjusted-data.RData")
load(data_file)

save_dir <- paste0(my_wd, "Analysis/Sweden/Outputs/Plots/")

# These samples are poorly behaved, drop
drop_sample_in_12 <- which((m$Sample.ID) %in% 1:2)
m <- m[-drop_sample_in_12, ]

# Drop the "patient 4" samples
patient_4 <- which(m$type == "Patient 4")
m <- m[-patient_4, ]

# Find the controls
negative_controls <- which(m$type == "Historical controls")
positive_controls <- which(m$type == "COVID")
controls <- c(negative_controls, positive_controls)

N <- nrow(m)
N_positive <- length(positive_controls)
N_negative <- length(negative_controls)
N_controls <- length(controls)

# Find which observations have a known label
fixed <- rep(0, N)
fixed[c(positive_controls, negative_controls)] <- 1

# Label
truth <- m$type %>%
  factor(levels = unique(m$type)) %>%
  as.numeric()

truth <- rep(2, N)
truth[negative_controls] <- 0
truth[positive_controls] <- 1

# Our batch variable
batches_present <- m$group %>%
  unique()

batch_vec <- m$group %>%
  as.factor() %>%
  as.numeric()

# The numeric data to model is contained in the features ``SPIKE`` and ``RBD``.
cols_of_interest <- c("SPIKE", "RBD")
cols_used <- which(colnames(m) %in% cols_of_interest)

# Log transform the data when in matrix form
X <- log(as.matrix(m[, ..cols_used])) %>%
  set_rownames(paste0("Person_", 1:N)) %>%
  set_colnames(cols_of_interest)

# Model inputs
K_max <- 2
B <- length(unique(batch_vec))
P <- ncol(X)


train_ind <- which(fixed == 1)
test_ind <- which(fixed == 0)
train <- train_data <- data.frame(
  SPIKE = X[train_ind, 1],
  RBD = X[train_ind, 2],
  Labels = factor(truth[train_ind])
  # Batch = factor(batch_vec[train_ind])
)

test <- test_data <- data.frame(
  SPIKE = X[-train_ind, 1],
  RBD = X[-train_ind, 2],
  Labels = factor(truth[-train_ind])
  # Batch = factor(batch_vec[-train_ind])
)

# The non-control data is used when estimating seroprevalence
non_control_data <- m[-controls, ]

# Find the allocations across time
# The non-control data has the week of the sample collection included in the
# sample ID
weeks_str <- non_control_data$Sample.ID

# Extract the week
weeks <- weeks_str %>%
  str_match_all("Wk(\\d{1,2})") %>%
  do.call(rbind, .) %>%
  `[`(, 2) %>%
  as.numeric()

N_by_week <- data.frame(N = table(weeks))


# Other model types to run via ``parsnip``
# engines <- c("randomForest", "kernlab", "glm")
model_names <- c("RF", "SVM", "LR", "LR_with_batch")
n_models <- length(model_names)
models <- vector("list", n_models) %>%
  set_names(model_names)


# This tibble will hold the outputs we're interested in to compare models.
output_df <- tibble(
  "Seroprevalence" = numeric(),
  "Time" = numeric(),
  "Model" = character()
)

for (.mod in model_names) {
  set.seed(1)
  
  if (.mod == "LR_with_batch") {
    new_train <- data.frame(train)
    new_test <- data.frame(test)
    B <- length(unique(batch_vec))
    
    scaled_data <- matrix(0, nrow = N, ncol = P)
    
    for (b in seq(1, B)) {
      batch_ind <- which(batch_vec == b)
      scaled_data[batch_ind, ] <- scale(X[batch_ind, ])
    }
    
    new_train[, 1:2] <- scaled_data[train_ind, ]
    new_test[, 1:2] <- scaled_data[test_ind, ]
    
    t_0 <- Sys.time()
    
    models[[.mod]] <- my_lr <- glm(Labels ~ .,
                                   family = binomial(link = "logit"),
                                   data = new_train
    )
    probs <- predict(my_lr, new_test, type = "response")
    predictions <- lr_results <- ifelse(probs > 0.5, 1, 0)
    
    probs_mat <- matrix(0, nrow = length(probs), ncol = 2)
    probs_mat[, 1] <- 1 - probs
    probs_mat[, 2] <- probs
    
    probs <- probs_mat
    
    # lr_acc <- calcAccuracy(lr_results, test$labels)
    t_1 <- Sys.time() - t_0
  }
  
  if (.mod == "LR") {
    t_0 <- Sys.time()
    
    models[[.mod]] <- my_lr <- glm(Labels ~ .,
                                   family = binomial(link = "logit"),
                                   data = train
    )
    probs <- predict(my_lr, test, type = "response")
    predictions <- lr_results <- ifelse(probs > 0.5, 1, 0)
    
    probs_mat <- matrix(0, nrow = length(probs), ncol = 2)
    probs_mat[, 1] <- 1 - probs
    probs_mat[, 2] <- probs
    
    probs <- probs_mat
    
    # lr_acc <- calcAccuracy(lr_results, test$labels)
    t_1 <- Sys.time() - t_0
  }
  
  if (.mod == "SVM") {
    t_0 <- Sys.time()
    
    # Radial Basis kernel "Gaussian"
    models[[.mod]] <- my_svm <- ksvm(Labels ~ ., data = train, prob.model = TRUE)
    probs <- predict(my_svm, test, type = "prob")
    predictions <- svm_pred <- predict(my_svm, test, type = "response")
    
    t_1 <- Sys.time() - t_0
  }
  
  if (.mod == "RF") {
    t_0 <- Sys.time()
    
    models[[.mod]] <- my_rf <- randomForest(Labels ~ ., data = train)
    
    probs <- rf_probs <- predict(my_rf, test, type = "prob")
    predictions <- rf_pred <- predict(my_rf, test, type = "class")
    
    t_1 <- Sys.time() - t_0
  }
  
  
  seroprevalence <- calcSeroprevalence(predictions, weeks, pos_value = 1)
  seroprevalence$Model <- .mod
  seroprevalence$N <- N_by_week$N.Freq
  
  output_df <- rbind(output_df, seroprevalence)
}

# output_df %>% 
#   ggplot(aes(x = Time, y = Seroprevalence, colour = Model)) +
#   geom_point() + 
#   geom_line() + 
#   ggthemes::scale_color_colorblind()


# === Write outputs ============================================================

# Save the various data frames that we use to check model behaviour
write.csv(output_df,
  paste0(my_wd, "/Analysis/Sweden/Outputs/otherModelsSeroprevalence.csv"),
  row.names = FALSE
)

