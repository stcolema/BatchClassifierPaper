# eliseLikeModelling.R
# Perform modelling of the data simulated from the model learnt on the Swedish 
# data and save various forms of information.
# 
# === Setup ====================================================================
library(BatchMixtureModel)
library(mdiHelpR)
library(ggplot2)
library(magrittr)
library(tidyr)
library(mcclust)
library(tibble)
library(dplyr)
library(data.table)
library(randomForest)
library(kernlab)
library(optparse)

# Collect arguments from the command line using the ``optparse`` package
inputArguments <- function() {
  option_list <- list(
    
    # Data to cluster
    optparse::make_option(c("-d", "--dir"),
                          type = "character",
                          default = "./",
                          help = "File path to BatchPaper directory [default= %default]",
                          metavar = "character"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}


# Various functions to compare methods in the simulation study
calcSensitivity <- function(pred, truth) {
  sum((pred == truth) & (truth == 1)) / sum(truth == 1)
}

calcSpecificity <- function(pred, truth) {
  sum((pred == truth) & (truth == 0)) / sum(truth == 0)
}

calcAccuracy <- function(pred, truth) {
  sum(pred == truth) / length(truth)
}

calcPrecision <- function(pred, truth) {
  sum((pred == truth) & (pred == 1)) / sum(pred == 1)
}

calcRecall <- function(pred, truth) {
  sum((pred == truth) & (pred == 1)) / sum(truth == 1)
}

# The F Beta score.
calcFBScore <- function(pred, truth, beta = 1) {
  precision <- calcPrecision(pred, truth)
  recall <- calcRecall(pred, truth)

  # The F beta score
  fb_score <- (1 + beta**2) * (precision * recall) / (beta**2 * precision + recall)
  fb_score
}

oneHotEncoding <- function(v) {
  v <- as.numeric(factor(v, levels = sort(unique(v))))

  N <- length(v)
  n_col <- length(unique(v))
  out <- matrix(0, ncol = n_col, nrow = N)
  for (n in 1:N) {
    out[n, v[n]] <- 1
  }
  out
}

calcDistance <- function(probs, truth) {
  truth <- oneHotEncoding(truth)
  sum((probs - truth)**2)
}


calcSeroprevalence <- function(pred, positive_label = 2) {
  (sum(pred == positive_label) / length(pred)) * 100
}


# Set the gglot2 theme
setMyTheme()

# Set the random seed
set.seed(1)

# Read in commandline arguments
args <- inputArguments()

# Directory to read from and to save Sparse matrices to
my_wd <- args$dir

# on_hpc <- FALSE
# my_wd <- "./"
# 
# if(on_hpc) {
#   my_wd <- "/home/sdc56/rds/hpc-work/BatchMixtureModel/"
# }

data_dir <- paste0(my_wd, "Simulations/ELISA_like/Data/")
data_files <- list.files(data_dir, pattern = "*.csv") %>%
  stringr::str_sort(numeric = T)

n_sims <- length(data_files) #


save_dir <- paste0(my_wd, "Simulations/ELISA_like/Output/")

# === Modelling ================================================================

# This tibble will hold the outputs we're interested in to compare models.
output_df <- tibble(
  "Distance" = numeric(),
  "Accuracy" = numeric(),
  "F1" = numeric(),
  "Sensitivity" = numeric(),
  "Specificity" = numeric(),
  "Seroprevalence" = numeric(),
  "Seroprevalence_diff" = numeric(),
  "BICmcmc" = numeric(),
  "Time" = numeric(),
  "Model" = character(),
  "Simulation" = numeric(),
  "Seed" = numeric()
)

# Types of mixture densities
# types <- c("MVN", "MSN", "MVT")
types <- c("MVN", "MVT")

# MCMC parameters
mu_proposal_window <- 0.10**2
cov_proposal_window <- 250
m_proposal_window <- 0.14**2
S_proposal_window <- 125
t_df_proposal_window <- 15

# Batch hyperparameters
rho <- 11
theta <- 5
m_scale <- 0.1

R <- 50000
thin <- 100
burn <- 20000
eff_burn <- burn / thin

# Other model types to run via ``parsnip``
# engines <- c("randomForest", "kernlab", "glm")
model_names <- c("RF", "SVM", "LR", "LR_with_batch")
n_models <- length(model_names)
models <- vector("list", n_models) %>%
  set_names(model_names)

acceptance_full <- tibble(
  "Mu_1" = numeric(),
  "Mu_2" = numeric(),
  "Sigma_1" = numeric(),
  "Sigma_2" = numeric(),
  "nu_1" = numeric(),
  "nu_2" = numeric(),
  "m_1" = numeric(),
  "m_2" = numeric(),
  "m_3" = numeric(),
  "m_4" = numeric(),
  "m_5" = numeric(),
  "m_6" = numeric(),
  "m_7" = numeric(),
  "S_1" = numeric(),
  "S_2" = numeric(),
  "S_3" = numeric(),
  "S_4" = numeric(),
  "S_5" = numeric(),
  "S_6" = numeric(),
  "S_7" = numeric(),
  "Chain" = numeric(),
  "Model" = character()
)

# === Main loop ================================================================

convergence_lst <- list()
lst_entry <- 0

# Number of chains run
n_chains <- 10 

for (ii in 1:n_sims) {
  my_data <- fread(paste0(data_dir, data_files[ii]))
  
  N <- nrow(my_data)
  P <- 2
  
  X <- cbind(my_data$X_observed, my_data$Y_observed)
  X_corr <- cbind(my_data$X_true, my_data$Y_true)

  batch_vec <- my_data$Batch
  B <- length(unique(batch_vec))
  
  truth <- my_data$Labels
  fixed <- my_data$Fixed
  N_fixed <- sum(fixed)

  true_seroprevalence <- (table(my_data$Labels)[2] / N) * 100

  train_ind <- which(my_data$Fixed == 1)
  test_ind <- which(my_data$Fixed == 0)

  for (type in types) {

    # This object will hold the sampled values that allows us to judge
    # convergence.

    for (jj in 1:n_chains) {
      set.seed(jj)

      # Initial labels
      initial_labels <- truth
      initial_labels[fixed == 0] <- sample(c(1, 2),
        size = N - N_fixed,
        prob = c(1, 1),
        replace = T
      )

      t0 <- Sys.time()
      samples <- batchSemiSupervisedMixtureModel(X,
        R,
        thin,
        initial_labels,
        fixed,
        batch_vec,
        type,
        mu_proposal_window = mu_proposal_window,
        cov_proposal_window = cov_proposal_window,
        m_proposal_window = m_proposal_window,
        S_proposal_window = S_proposal_window,
        t_df_proposal_window = t_df_proposal_window,
        rho = rho,
        theta = theta,
        m_scale = m_scale
      )
      t_1 <- Sys.time() - t0

      # For saving numbers with a decimal place
      m_scale_sci_notation <- formatC(m_scale, format = "e", digits = 0)

      # Save the samples for the current chain
      saveRDS(samples, file = paste0(
        my_wd,
        "/Simulations/ELISA_like/Output/Chains/Sim", 
        ii,
        "/",
        type, 
        "_chain_",
        jj,
        "_m_scale_",
        m_scale_sci_notation,
        "_rho_",
        rho,
        "_theta_",
        theta,
        ".rds"
      ))


      samples_df <- samplesToDF(samples, type,
        R = R,
        thin = thin,
        keep_allocation = F
      )

      acceptance_df <- collectAcceptanceRates(samples, type)
      acceptance_df$Model <- type
      acceptance_df$Chain <- jj
      acceptance_df$Simulation <- ii

      if (type != "MVT") {
        acceptance_df$nu_1 <- NA
        acceptance_df$nu_2 <- NA
      }

      acceptance_full <- rbind(acceptance_full, acceptance_df)


      samples_df$Model <- type
      samples_df$Chain <- jj
      samples_df$Simulation <- ii

      if (type == "MVN") {
        samples_df$t_df_1 <- samples_df$t_df_2 <- NA
      }

      lst_entry <- lst_entry + 1
      convergence_lst[[lst_entry]] <- samples_df

      probs <- calcAllocProb(samples$alloc, eff_burn)
      predictions <- predictClass(probs)

      recall <- calcRecall(predictions[test_ind], truth[test_ind])
      accuracy <- calcAccuracy(predictions[test_ind], truth[test_ind])
      precision <- calcPrecision(predictions[test_ind], truth[test_ind])
      f1Score <- calcFBScore(predictions[test_ind], truth[test_ind])
      sensitivity <- calcSensitivity(predictions[test_ind], truth[test_ind])
      specificity <- calcSpecificity(predictions[test_ind], truth[test_ind])
      distance <- calcDistance(probs[test_ind, ], truth[test_ind])
      bicmcmc <- max(samples$BIC[-c(1:eff_burn)])

      seroprevalence <- calcSeroprevalence(predictions)
      difference_in_seroprevalence <- seroprevalence - true_seroprevalence

      new_out <- tibble(
        "Distance" = distance,
        "Accuracy" = accuracy,
        "F1" = f1Score,
        "Sensitivity" = sensitivity,
        "Specificity" = specificity,
	"Seroprevalence" = seroprevalence,
        "Seroprevalence_diff" = difference_in_seroprevalence,
        "BICmcmc" = bicmcmc,
        "Time" = t_1,
        "Model" = type,
        "Simulation" = ii,
        "Seed" = jj
      )

      output_df <- rbind(output_df, new_out)
    }
  }

  my_data$Labels <- factor(my_data$Labels)
  train <- my_data[train_ind, c(1, 4, 5)]
  test <- my_data[test_ind, c(1, 4, 5)]

  for (.mod in model_names) {
    set.seed(1)

    if (.mod == "LR_with_batch") {
      new_train <- data.frame(train)
      new_test <- data.frame(test)
      B <- length(unique(batch_vec))

      scaled_data <- matrix(0, nrow = N, ncol = P)

      for (b in seq(1, B)) {
        batch_ind <- which(batch_vec == b)
        scaled_data[batch_ind, ] <- scale(my_data[batch_ind, c(4, 5)])
      }

      new_train[, 2:3] <- scaled_data[train_ind, ]
      new_test[, 2:3] <- scaled_data[test_ind, ]
      
      t_0 <- Sys.time()

      models[[.mod]] <- my_lr <- glm(Labels ~ .,
        family = binomial(link = "logit"),
        data = new_train
      )
      probs <- predict(my_lr, new_test, type = "response")
      predictions <- lr_results <- ifelse(probs > 0.5, 2, 1)

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
      predictions <- lr_results <- ifelse(probs > 0.5, 2, 1)

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

    recall <- calcRecall(predictions, as.numeric(test$Labels))
    accuracy <- calcAccuracy(predictions, as.numeric(test$Labels))
    precision <- calcPrecision(predictions, as.numeric(test$Labels))
    f1Score <- calcFBScore(predictions, as.numeric(test$Labels))
    sensitivity <- calcSensitivity(predictions, as.numeric(test$Labels))
    specificity <- calcSpecificity(predictions, as.numeric(test$Labels))
    distance <- calcDistance(as.numeric(probs), as.numeric(test$Labels))

    seroprevalence <- calcSeroprevalence(predictions)
    difference_in_seroprevalence <- seroprevalence - true_seroprevalence
    
    new_out <- tibble(
      "Distance" = distance,
      "Accuracy" = accuracy,
      "F1" = f1Score,
      "Sensitivity" = sensitivity,
      "Specificity" = specificity,
      "Seroprevalence" = seroprevalence,
      "Seroprevalence_diff" = difference_in_seroprevalence,
      "BICmcmc" = NA,
      "Time" = t_1,
      "Model" = .mod,
      "Simulation" = ii,
      "Seed" = ii
    )

    output_df <- rbind(output_df, new_out)
  }
}


# === Acceptance ===============================================================

data.table::fwrite(acceptance_full, 
  paste0(my_wd, 
    "/Simulations/ELISA_like/Output/SimAcceptanceRates.csv"
  )
)

# === Convergence ==============================================================

saveRDS(convergence_lst, 
  file = paste0(my_wd, "/Simulations/ELISA_like/Output/ConvergenceList.rds")
)

# Save the sampled parameters in data frame format for ease of use
param_df <- do.call(rbind, convergence_lst)
data.table::fwrite(param_df, 
  paste0(my_wd, "/Simulations/ELISA_like/Output/SampledParameters.csv")
)

# === Model comparison =========================================================

data.table::fwrite(output_df, 
  paste0(my_wd, "/Simulations/ELISA_like/Output/SimModelPerformance.csv")
)

