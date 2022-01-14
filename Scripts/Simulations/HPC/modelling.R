# modelling.R
# Perform modelling of the simulated data and save various forms of information
# to be used by ``visualisation.R``.
# 
# === Swetup ===================================================================

# install and call various packages
cran_packages <- c(
#  "foreach",
  "tidyr",
  "ggplot2",
  "dplyr",
  "magrittr",
  "devtools",
  "tibble",
  "randomForest",
  "kernlab",
  "data.table"
)

install.packages(setdiff(cran_packages, rownames(installed.packages())))

# library(foreach)
library(tidyr)
library(ggplot2)
library(dplyr)
library(magrittr)
library(tibble)
library(coda)
library(mdiHelpR)
library(BatchMixtureModel)
library(randomForest)
library(kernlab)

wd <- "./"
on_hpc <- TRUE
if(on_hpc) {
  my_wd <- "/home/sdc56/rds/hpc-work/BatchMixtureModel/"
}

save_dir <- paste0(my_wd, "Simulations/Output/")

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

calcDistance <- function(probs, truth) {
  truth <- oneHotEncoding(truth)
  sum((probs - truth)**2)
}

# Function to predict a class and find allocation probability
allocProb <- function(samples, burn) {
  if (burn > 0) {
    samples <- samples[, , -c(1:burn)]
  }
  prob <- rowSums(samples, dims = 2) / dim(samples)[3]
  prob
}


predictedClass <- function(prob) {
  pred_cl <- apply(prob, 1, which.max)

  pred_cl
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

calcSeroprevalence <- function(pred, positive_label = 2) {
  (sum(pred == positive_label) / length(pred)) * 100
}


# generateInitialLabels <- function(labels, fixed) {
#   N <- length(labels)
#   N_fixed <- sum(fixed)
#   N_unfixed <- N - N_fixed
#   ratio <- table(labels[fixed == 1]) / N_fixed
#   labels[fixed == 0] <- sample(unique(labels), N_unfixed, replace = T, prob = ratio)
#   labels
# }


# Set the gglot2 theme
setMyTheme()

# Set the random seed
set.seed(1)

# Save the generated files to
save_path <- paste0(my_wd, "Simulations/Output/")
dir.create(save_path, showWarnings = F)

# The files we save
acceptance_file <- paste0(save_path, "SimAcceptanceRates.csv")
convergence_file <- paste0(save_path, "ConvergenceList.rds")
param_file <- paste0(save_path, "SampledParameters.csv")

# The generated data and the description tibble
data_path <- paste0(my_wd, "Simulations/Generated_data/")
scenario_descriptions <- readRDS(paste0(data_path, "scenario_descriptions.rds"))

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
  "Scenario" = character(),
  "Simulation" = numeric(),
  "Seed" = numeric()
)

# Types of mixture densities
# types <- c("MVN", "MSN", "MVT")
types <- c("MVN", "MVT")

# Number of chains run
n_chains <- 10

# MCMC parameters
mu_proposal_window <- 0.45**2
cov_proposal_window <- 75
m_proposal_window <- 0.5**2
S_proposal_window <- 12
t_df_proposal_window <- 4
phi_proposal_window <- 1.2**2
R <- 15000
thin <- 25
burn <- 7500
eff_burn <- burn / thin

# Prior distirbution parameters for the batch effects
rho <- list("MVN" = 21, "MVT" = 21)
theta <- list("MVN" = 10, "MVT" = 10)
lambda <- list("MVN" = 1, "MVT" = 1)

n_scn <- nrow(scenario_descriptions)
n_sim <- 10

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
  "phi_1" = numeric(),
  "phi_2" = numeric(),
  "m_1" = numeric(),
  "m_2" = numeric(),
  "m_3" = numeric(),
  "m_4" = numeric(),
  "m_5" = numeric(),
  "S_1" = numeric(),
  "S_2" = numeric(),
  "S_3" = numeric(),
  "S_4" = numeric(),
  "S_5" = numeric(),
  "Chain" = numeric(),
  "Model" = character(),
  "Scenario" = character()
)

# === Main loop ================================================================

convergence_lst <- list()
lst_entry <- 0

# n_scn <- n_sim <- 2

# n_sim <- 5
# n_chains <- 3
# n_scn <- 1

for (ii in 1 : n_scn) {
# for (ii in c(1, 5)) {
  N <- scenario_descriptions$N[ii]
  P <- scenario_descriptions$P[ii]
  scn <- scenario_descriptions$Scenario[ii]
  curr_dir <- paste0(data_path, scn)

  # for (jj in 11:12) {
  for (jj in 1:n_sim) {
    my_data <- read.csv(paste0(curr_dir, "/seed_", jj, ".csv"))
    X <- cbind(my_data$X_observed, my_data$Y_observed)
    X_corr <- cbind(my_data$X_true, my_data$Y_true)

    batch_vec <- my_data$batch
    truth <- my_data$labels
    fixed <- my_data$fixed

    true_seroprevalence <- 100 * sum(my_data$labels == 2) / N
    
    train_ind <- which(my_data$fixed == 1)
    test_ind <- which(my_data$fixed == 0)

    for (type in types) {

      # This object will hold the sampled values that allows us to judge
      # convergence.

      for (kk in 1:n_chains) {
        set.seed(kk)
        initial_labels <- generateInitialLabels(truth, fixed)

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
          phi_proposal_window = phi_proposal_window
        )
        t_1 <- Sys.time() - t0

        samples_df <- samplesToDF(samples, type,
          R = R,
          thin = thin,
          keep_allocation = F
        )

        acceptance_df <- collectAcceptanceRates(samples, type)
        acceptance_df$Model <- type
        acceptance_df$Chain <- kk
        acceptance_df$Simulation <- jj
        acceptance_df$Scenario <- scn

        if (type != "MSN") {
          acceptance_df$phi_1 <- NA
          acceptance_df$phi_2 <- NA
        }
        if (type != "MVT") {
          acceptance_df$nu_1 <- NA
          acceptance_df$nu_2 <- NA
        }

        # acceptance_full

        # if(kk == 1) {
        #   acceptance_full <- acceptance_df
        # } else {
        acceptance_full <- rbind(acceptance_full, acceptance_df)
        # }

        samples_df$Model <- type
        samples_df$Chain <- kk
        samples_df$Simulation <- jj
        samples_df$Scenario <- scn

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
          "Scenario" = scn,
          "Simulation" = jj,
          "Seed" = kk
        )

        output_df <- rbind(output_df, new_out)
      }
    }

    my_data$labels <- factor(my_data$labels)
    train <- my_data[train_ind, c(1, 4, 5)]
    test <- my_data[test_ind, c(1, 4, 5)]

    # .recipe <- recipe(labels ~ ., data = train) %>%
    #   step_center(all_numeric()) %>%
    #   step_scale(all_numeric()) %>%
    #   prep(training = train)
    #
    # .testing <- bake(.recipe, test)
    # .training <- juice(.recipe)

    for (.mod in model_names) {
      
      set.seed(1)
      
      if(.mod == "LR_with_batch") {
        
        new_train <- train
        new_test <- test
        B <- scenario_descriptions$B[ii]
        
        scaled_data <- matrix(0, nrow = N,ncol = P)
        
        for(b in seq(1, B)) {
          
          batch_ind <- which(batch_vec == b)
          scaled_data[batch_ind, ] <- scale(my_data[batch_ind, c(4, 5)])
          
          # batch_train_ind <- which(batch_vec[train_ind] == b)
          # batch_test_ind <- which(batch_vec[test_ind] == b)
          # n_bt <- length(batch_test_ind)
          # 
          # scaled_batch_test <- matrix(0, nrow = n_bt,ncol = P)
          # 
          # batch_train <- train[batch_train_ind, 2:3]
          # batch_test <- test[batch_test_ind, 2:3]
          # 
          # scaled_batch_train <- scale(batch_train)
          # 
          # scaled_centre <- attr(scaled_batch_train, "scaled:center")
          # scaled_scale <- attr(scaled_batch_train, "scaled:scale")
          # 
          # for(p in seq(1, P)) {
          #   scaled_batch_test[, p] <- (batch_test[, p] - scaled_centre[p]) / scaled_scale[p]
          # }
          # 
          # new_train[batch_train_ind, 2:3] <- scaled_batch_train
          # new_test[batch_test_ind, 2:3] <- scaled_batch_test
          
        }
        
        new_train[, 2:3] <- scaled_data[train_ind, ]
        new_test[, 2:3] <- scaled_data[test_ind, ]
        
        
        t_0 <- Sys.time()
        
        models[[.mod]] <- my_lr <- glm(labels ~ .,
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

        models[[.mod]] <- my_lr <- glm(labels ~ .,
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
        models[[.mod]] <- my_svm <- ksvm(labels ~ ., data = train,prob.model=TRUE)
        probs <- predict(my_svm, test, type = "prob")
        predictions <- svm_pred <- predict(my_svm, test, type = "response")
        # svm_acc <- calcAccuracy(svm_pred, test$labels)

        # models[[.mod]] <- svm_poly(mode = "classification") %>%
        #   set_engine("kernlab") %>%
        #   fit(labels ~ ., data = .training)
        t_1 <- Sys.time() - t_0
      }

      if (.mod == "RF") {
        t_0 <- Sys.time()

        models[[.mod]] <- my_rf <- randomForest(labels ~ ., data = train)

        probs <- rf_probs <- predict(my_rf, test, type = "prob")
        predictions <- rf_pred <- predict(my_rf, test, type = "class")
        # rf_acc <- calcAccuracy(rf_pred, test$labels)

        # models[[.mod]] <- rand_forest(trees = 100, mode = "classification") %>%
        #   set_engine("ranger") %>%
        #   fit(labels ~ ., data = .training)
        t_1 <- Sys.time() - t_0
      }

      # predictions <- predict(models[[.mod]], .testing, type = "class")
      # probs <- predict(models[[.mod]], .testing, type = "prob")

      recall <- calcRecall(predictions, as.numeric(test$labels))
      accuracy <- calcAccuracy(predictions, as.numeric(test$labels))
      precision <- calcPrecision(predictions, as.numeric(test$labels))
      f1Score <- calcFBScore(predictions, as.numeric(test$labels))
      sensitivity <- calcSensitivity(predictions, as.numeric(test$labels))
      specificity <- calcSpecificity(predictions, as.numeric(test$labels))
      distance <- calcDistance(as.numeric(probs), as.numeric(test$labels))

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
        "Scenario" = scn,
        "Simulation" = jj,
        "Seed" = jj
      )

      output_df <- rbind(output_df, new_out)
    }
  }
}

# === Acceptance ===============================================================

# long_acceptance_full <- acceptance_full %>%
#   select(-c(phi_1, phi_2)) %>%
#   pivot_longer(-c(Model, Chain, Scenario, Simulation), values_to = "Acceptance_rate", names_to = "Parameter")
# 
# # Should be the number of scenarios x number of saved samples x number of clusters
# # (i.e. )
# MVN_df <- which(is.na(long_acceptance_full$Acceptance_rate))

# if(length(MVN_df) > n_scn * n_chains * 2 *

# p_acceptance <- long_acceptance_full[-MVN_df, ] %>%
#   ggplot(aes(x = Parameter, y = Acceptance_rate)) +
#   geom_boxplot() +
#   facet_grid(Scenario ~ Model, scales = "free_x") + 
#   labs(title = "Parameter acceptance rates for Bayesian models",
#        subtitle = "All simulations and chains pooled"
#      )
# 
# ggsave("./Simulations/SimAcceptanceRates.png", p_acceptance)

data.table::fwrite(acceptance_full, acceptance_file)

# cols_ignored <- which(colnames(acceptance_full) %in% c("Model", "Chain", "Simulation", "Scenario"))
# 
# drop_seeds <- acceptance_full[, -cols_ignored] %>%
#   apply(1, function(x) {
#     any(x < 0.1)
#   }) %>%
#   which()
# 
# dropped_models <- acceptance_full[, cols_ignored]
# dropped_models$Dropped <- FALSE
# dropped_models$Dropped[drop_seeds] <- TRUE

# === Convergence ==============================================================

# used_lst <- convergence_lst
# if (length(drop_seeds) > 0) {
#   used_lst <- convergence_lst[-drop_seeds]
# }

#acceptance_file <- paste0(save_dir, "SimAcceptanceRates.csv")
#convergence_file <- paste0(save_dir, "ConvergenceList.rds")
#param_file <- paste0(save_dir, "SampledParameters.csv")

saveRDS(convergence_lst, file = convergence_file)

#   .curr_df <- convergence_lst[[i]]
#   # .curr_df$Simulation <- i
# 
#   if (i == 1) {
#     param_df <- .curr_df
#   } else {
#     param_df <- rbind(param_df, .curr_df)
#   }
# }

# Save the sampled parameters in data frame format for ease of use
param_df <- do.call(rbind, convergence_lst)
data.table::fwrite(param_df, param_file)

# === Likelihood test ==========================================================

# likelihood_bounds <- param_df %>% 
#   filter(Iteration > burn) %>% 
#   select(Iteration, Scenario, Simulation, Chain, Model, Complete_likelihood) %>% 
#   group_by(Scenario, Simulation, Model) %>% 
#   summarise(Mean_likelihood = mean(Complete_likelihood), 
#             SD_likelihood = sd(Complete_likelihood),
#             UB = Mean_likelihood + 2 * SD_likelihood,
#             LB = Mean_likelihood - 2 * SD_likelihood)
# 
# 
# 
# dropped_model_runs <- param_df %>%
#   filter(Iteration > burn) %>% 
#   select(Iteration, Scenario, Simulation, Chain, Model, Complete_likelihood) %>% 
#   group_by(Scenario, Simulation, Model, Chain) %>% 
#   summarise(Mean_chain_likelihood = median(Complete_likelihood)) %>% 
#   left_join(likelihood_bounds) %>% 
#   filter(Mean_chain_likelihood > UB | Mean_chain_likelihood < LB)
# 
# dropped_model_runs$Dropped_likelihood <- T
# 
# new_dropped <- left_join(dropped_models, dropped_model_runs, 
#           by = c("Scenario",
#                  "Simulation",
#                  "Model",
#                  "Chain")
#           ) %>% 
#   select(-c(Mean_chain_likelihood, Mean_likelihood, SD_likelihood, UB, LB)) 
# 
# new_dropped$Dropped_likelihood[is.na(new_dropped$Dropped_likelihood)] <- F
# new_dropped$Dropped <- new_dropped$Dropped | new_dropped$Dropped_likelihood
# 
# new_dropped <- new_dropped %>% select(-Dropped_likelihood)

# lkl2 <- param_df %>% 
#   # filter(Simulation == 1, Scenario == "MVT", Model == "MVT") %>% 
#   group_by(Chain, Simulation, Scenario, Model) %>% 
#   summarise(Likelihood_sum = sum(Complete_likelihood)) %>% 
#   ungroup() %>% 
#   group_by(Simulation, Scenario, Model) %>% 
#   summarise(Median_likelihood = median(Likelihood_sum), 
#             SD_likelihood = sd(Likelihood_sum),
#             LB = Median_likelihood - SD_likelihood,
#             UB = Median_likelihood+ SD_likelihood)
# 
# lkl2$Median_likelihood <-median(lkl2$Likelihood_sum)
# lkl2$SD_likelihood <-sd(lkl2$Likelihood_sum)
# lkl2$UB <- lkl2$Median_likelihood + lkl2$SD_likelihood
# lkl2$LB <- lkl2$Median_likelihood - lkl2$SD_likelihood
# lkl2$Dropped <- lkl2$Likelihood_sum < lkl2$LB | lkl2$Likelihood_sum > lkl2$UB


# param_df %>%
#   # filter(Iteration > burn) %>%
#   ggplot(aes(
#     x = Iteration,
#     y = Observed_likelihood,
#     group = interaction(Simulation, Chain),
#     colour = Model
#   )) +
#   geom_line() +
#   facet_grid(Scenario ~ Model)
# 
# param_df %>%
#   filter(Iteration > burn) %>%
#   ggplot(aes(
#     x = Iteration,
#     y = Complete_likelihood,
#     group = interaction(Simulation, Chain),
#     colour = Model
#   )) +
#   geom_line() +
#   facet_grid(Scenario ~ Model)

# param_df %>%
#   filter(Iteration > burn, Model == "MVN") %>%
#   ggplot(aes(
#     x = Iteration,
#     y = Complete_likelihood,
#     group = interaction(Simulation, Chain),
#     colour = factor(Chain)
#   )) +
#   geom_line() +
#   facet_wrap(~Chain) +
#   labs(title = "MVN in data generated from MVT",
#        subtitle = paste("Simulation", jj))

# ggsave("./MVNLikelihoodsScnMVTSim4.png")  

# facet_grid(Simulation + Scenario ~ Model) 


# param_df %>%
#   filter(Iteration > burn, Model == "MVT", Scenario == "MVT") %>%
#   dplyr::select(t_df_1, t_df_2, Iteration, Simulation, Chain, Scenario) %>%
#   pivot_longer(c(t_df_1, t_df_2), names_to = "Parameter") %>%
#   ggplot(aes(x = Iteration, y = value, colour = factor(Chain))) +
#   geom_line() +
#   facet_grid(Parameter ~ Simulation,  labeller = label_both) +
#   theme(axis.text.x = element_text(angle = 30))

# ggsave("../SampledDegreesOfFreedomMVTSim.png", width = 12, height = 8)

# param_df %>%
#   filter(Iteration > burn, Model == "MVT") %>%
#   ggplot(aes(x = Iteration, y = t_df_1, group = interaction(Simulation, Chain))) +
#   geom_line() +
#   facet_wrap(~Scenario)
#
# param_df %>%
#   filter(Iteration > burn, Model == "MVT") %>%
#   ggplot(aes(x = Iteration, y = t_df_2, group = interaction(Simulation, Chain))) +
#   geom_line() +
#   facet_wrap(~Scenario)
#
#
# param_df %>%
#   filter(Iteration > burn, Model == "MVT", Scenario == "MVT", Simulation == 1) %>%
#   ggplot(aes(x = Iteration, y = t_df_2, group = interaction(Simulation, Chain))) +
#   geom_line() +
#   facet_wrap(~Chain)

# param_df %>%
#   filter(Model == "MVT", Simulation == 1, Chain == 1, Iteration > burn) %>%
#   select(Iteration, Scenario, Simulation, Chain, t_df_1, t_df_2) %>%
#   group_by(Scenario, Simulation, Chain) %>%
#   pivot_longer(c(t_df_1, t_df_2), names_to = "Parameter") %>%
#   ggplot(aes(x = Iteration, y = value, colour = Parameter, group = Chain)) +
#   geom_line() +
#   facet_grid(Parameter ~ Scenario) +
#   labs(
#     title = "Simulations",
#     subtitle = "Sampled degrees of freedom in a single chain"
#   )
# 
# param_df %>%
#   filter(Model == "MVT", Iteration > burn, Simulation == 1) %>%
#   select(Iteration, Scenario, Simulation, Chain, t_df_1, t_df_2) %>%
#   group_by(Scenario, Simulation, Chain) %>%
#   pivot_longer(c(t_df_1, t_df_2), names_to = "Parameter") %>%
#   ggplot(aes(x = value, y = Chain, colour = Parameter, group = Chain)) +
#   geom_boxplot() +
#   facet_grid(Parameter ~ Scenario) +
#   labs(
#     title = "Simulations",
#     subtitle = "Sampled degrees of freedom in a single chain"
#   )
# 
# # ggsave("../sims_sampled_degrees_of_freedom_boxplot.png")
# 
# param_df %>%
#   filter(Iteration > burn) %>%
#   ggplot(aes(
#     x = Iteration,
#     y = Complete_likelihood,
#     group = Chain,
#     colour = Model
#   )) +
#   geom_line() +
#   facet_grid(Model + Simulation ~ Scenario, labeller = label_both)
# 
# ggsave("./Complete_liklihood.png")
# 
# # param_df %>%
# #   filter(Iteration > burn) %>%
# #   ggplot(aes(x = t_df_1, y = factor(Chain))) +
# #   geom_boxplot() +
# #   facet_wrap(~Scenario)
# #   # xlim(c(3.5, 4.0))
# 
# # ggsave("../Worying_df.png")
# #
# # used_lst[[1]][1:2, c(1, 78:81)]
# #
# # used_lst[[11]]
# #
# # colnames(used_lst[[1]])
# # colnames(used_lst[[11]])

# === Geweke within-chain ======================================================
# 
# mcmc_lst <- convergence_lst %>%
#   lapply(function(x) {
#     
#     # We use the likelihoods to decide if the chains are converged
#     cols_used <- c("Complete_likelihood", "Observed_likelihood")
#     
#     y <- as.mcmc(x[-c(1:eff_burn), colnames(x) %in% cols_used])
#     
#     # d <- ncol(x)
#     # if (x$Model[1] == "MVN") {
#     #   y <- as.mcmc(x[-c(1:eff_burn), -c(1, (d - 5):d)])
#     # } else {
#     #   y <- as.mcmc(x[-c(1:eff_burn), -c(1, (d - 3):d)])
#     # }
#     y
#   })
# 
# geweke_lst <- mcmc_lst %>%
#   # lapply(gewekePlot)
#   lapply(function(x) {
#     frac_1 <- 0.1
#     frac_2 <- 0.5
#     n_bins <- 20
#     p_value_threshold <- 0.05
#     threshold_line_colour <- "grey"
#     plt_title <- "Geweke diagnostic plot"
#     
#     # The preferred object type for interacting with coda functions
#     x <- coda::as.mcmc.list(x)
#     
#     # The vector of start iterations to calculate the Geweke statistic for
#     start_iter_vec <- floor(seq(
#       from = stats::start(x),
#       to = (stats::start(x) + stats::end(x)) / 2,
#       length = n_bins
#     ))
#     
#     # The matrix that will hold the Geweke stat
#     geweke_mat <- matrix(nrow = length(start_iter_vec), ncol = coda::nvar(x), dimnames = list(start_iter_vec, coda::varnames(x)))
#     
#     for (n in 1:length(start_iter_vec)) {
#       curr_geweke_diag <- coda::geweke.diag(stats::window(x, start = start_iter_vec[n]),
#                                             frac1 = frac_1,
#                                             frac2 = frac_2
#       )
#       
#       geweke_mat[n, ] <- curr_geweke_diag[[1]]$z
#     }
#     
#     # The 1.96 threshold for 0.05 significance on a standard normal distribution
#     c_limit <- stats::qnorm(1 - p_value_threshold / 2)
#     
#     # The variables to gather when moving from wide to long data (these are our
#     # parameters)
#     vars_to_gather <- coda::varnames(x)
#     
#     # The data.frame we will plot (transform to long data to use the ggplot2
#     # framework)
#     geweke_df <- data.frame(Start_iteration = start_iter_vec) %>%
#       cbind(geweke_mat) %>%
#       tidyr::gather_("Parameter", "Geweke_statistic", vars_to_gather)
#     
#     geweke_df
#   })
# 
# # geweke_df <- do.call(rbind, geweke_lst)
# 
# for (i in 1:length(geweke_lst)) {
#   curr_chain <- convergence_lst[[i]]$Chain %>% unique()
#   curr_sim <- convergence_lst[[i]]$Simulation %>% unique()
#   curr_scn <- convergence_lst[[i]]$Scenario %>% unique()
#   curr_model <- convergence_lst[[i]]$Model %>% unique()
#   if (length(curr_chain) > 1) {
#     stop("Too many chains.")
#   }
#   if (length(curr_sim) > 1) {
#     stop("Too many sims")
#   }
#   if (length(curr_scn) > 1) {
#     stop("Too many scenarios")
#   }
#   if (length(curr_model) > 1) {
#     stop("Too many models.")
#   }
#   geweke_lst[[i]]$Chain <- curr_chain # (1:n_chains)[i]
#   # geweke_lst[[i]]$Iteration <- used_lst[[i]]$Iteration
#   geweke_lst[[i]]$Simulation <- curr_sim
#   geweke_lst[[i]]$Scenario <- curr_scn
#   geweke_lst[[i]]$Model <- curr_model
#   
#   if (i == 1) {
#     geweke_df <- geweke_lst[[i]]
#   } else {
#     geweke_df <- rbind(geweke_df, geweke_lst[[i]])
#   }
# }
# 
# geweke_df$Start_iteration <- geweke_df$Start_iteration * thin + burn
# 
# 
# # Does the distribution of the Geweke statistics pass the Shapiro Wilks test of 
# # normality
# passed_geweke_table <- geweke_df %>% 
#   group_by(Chain, Simulation, Scenario, Model) %>% 
#   summarise(Normal = shapiro.test(Geweke_statistic)$p.value > 0.05) 
# 
# # Left join the ``Normal`` variable to the main df 
# geweke_df <- left_join(geweke_df, passed_geweke_table)
# 
# # write.csv(geweke_df, "./Simulations/GewekeDF.csv")
# 
# # The 1.96 threshold for 0.05 significance on a standard normal distribution
# c_limit <- stats::qnorm(1 - 0.05 / 2)
# 
# params_used <- c("Complete_likelihood") # unique(geweke_df$Parameter)[35:74]
# 
# # Plot the kept chains
# p_geweke_dropped <- geweke_df %>%
#   filter(Parameter %in% params_used, Normal) %>%
#   ggplot(aes(
#     x = Start_iteration,
#     y = Geweke_statistic,
#     colour = as.factor(Chain),
#     group = interaction(Chain, Parameter)
#   )) +
#   geom_line() +
#   geom_hline(yintercept = c_limit, linetype = "dashed", color = "red") +
#   geom_hline(yintercept = -c_limit, linetype = "dashed", color = "red") +
#   facet_grid(Simulation ~ Scenario + Model, labeller = label_both) +
#   labs(title = "Geweke statistics for the complete likelihood",
#        subtitle = "Non-converged chains dropped") +
#   theme(axis.text.x = element_text(angle = 30))
# 
# p_geweke_dropped
#   
# # Plot the kept chains
# p_geweke_included <- geweke_df %>%
#   filter(Parameter %in% params_used) %>%
#   ggplot(aes(
#     x = Start_iteration,
#     y = Geweke_statistic,
#     colour = as.factor(Chain),
#     group = interaction(Chain, Parameter)
#   )) +
#   geom_line() +
#   geom_hline(yintercept = c_limit, linetype = "dashed", color = "red") +
#   geom_hline(yintercept = -c_limit, linetype = "dashed", color = "red") +
#   facet_grid(Simulation ~ Scenario + Model, labeller = label_both) +
#   labs(title = "Geweke statistics for the complete likelihood",
#        subtitle = "Non-converged chains included")
# 
#  
# # 
# # # write.csv(output_df, "../SimOutputBeforeMangling.csv")
# # 
# # #
# # # my_df <- tibble(
# # #   iteration = convergence_lst[[1]]$Iteration,
# # #   # model_1 = convergence_lst[[1]]$Likelihood,
# # #   model_2 = convergence_lst[[2]]$Likelihood,
# # #   model_3 = convergence_lst[[3]]$Likelihood
# # # ) %>%
# # #   pivot_longer(-iteration, names_to = "model", values_to = "likelihood")
# # #
# # # my_df %>%
# # #   ggplot(aes(x = iteration, y = likelihood, colour = model)) +
# # #   geom_line() +
# # #   ggthemes::scale_color_colorblind()
# # #
# # # train <- my_data[my_data$fixed == 1, c(1,4,5)]
# # # test <- my_data[my_data$fixed == 0, c(1,4,5)]
# # #
# # # cross_val_tbl <- vfold_cv(train, v = 10)
# # #
# # # recipe_rf <- function(dataset) {
# # #   recipe(labels ~ ., data = dataset) %>%
# # #     # step_string2factor(all_nominal(), -all_outcomes()) %>%
# # #     # step_dummy(all_nominal(), -all_outcomes()) %>%
# # #     step_center(all_numeric()) %>%
# # #     step_scale(all_numeric()) %>%
# # #     prep(training = dataset)
# # # }
# # #
# # # recipe_rf(my_data)
# # #
# # # rf_fun <- function(split, id, try, tree) {
# # #
# # #   analysis_set <- split %>% analysis()
# # #   analysis_prepped <- analysis_set %>% recipe_rf()
# # #   analysis_baked <- analysis_prepped %>% bake(new_data = analysis_set)
# # #
# # #   model_rf <-
# # #     rand_forest(
# # #       mode = "classification",
# # #       mtry = try,
# # #       trees = tree
# # #     ) %>%
# # #     set_engine("ranger",
# # #                importance = "impurity"
# # #     ) %>%
# # #     fit(factor(labels) ~ ., data = analysis_baked)
# # #
# # #   assessment_set <- split %>% assessment()
# # #   assessment_prepped <- assessment_set %>% recipe_rf()
# # #   assessment_baked <- assessment_prepped %>% bake(new_data = assessment_set)
# # #
# # #   tibble(
# # #     "id" = id,
# # #     "truth" = factor(assessment_baked$labels),
# # #     "prediction" = model_rf %>%
# # #       predict(new_data = assessment_baked) %>%
# # #       unlist()
# # #   )
# # #
# # # }
# # #
# # # pred_rf <- map2_df(
# # #   .x = cross_val_tbl$splits,
# # #   .y = cross_val_tbl$id,
# # #   ~ rf_fun(split = .x, id = .y, try = 2, tree = 200)
# # # )
# # #
# # # pred_rf %>%
# # #   conf_mat(truth, prediction) %>%
# # #   summary() %>%
# # #   select(-.estimator) %>%
# # #   filter(.metric %in%
# # #            c("accuracy", "precision", "recall", "f_meas")) %>%
# # #   knitr::kable()
# # #
# # #
# # # my_data$labels <- factor(my_data$labels)
# # #
# # # recipe_rf(my_data)
# # #
# # # recipe(labels ~ X + Y, data = my_data)  %>%
# # #   step_string2factor(all_nominal(), -all_outcomes()) %>%
# # #   step_dummy(all_nominal(), -all_outcomes()) %>%
# # #   step_center(all_numeric()) %>%
# # #   step_scale(all_numeric())
# # #
# # # drop_seeds <- acceptance_full %>%
# # #   select(-c(Model, Chain, Scenario)) %>%
# # #   apply(1, function(x) {
# # #     any(x > 0.8 | x < 0.2)
# # #   }) %>%
# # #   which()

# === Model comparison =========================================================
# 
# # dropped_models$Seed <- dropped_models$Chain
# 
# mixture_output <- output_df[!is.na(output_df$BICmcmc), ]
# mixture_output$Chain <- mixture_output$Seed
# kept_mixture_output <- merge(mixture_output, new_dropped,
#   by = c("Model", "Scenario", "Chain", "Simulation")
# ) %>%
#   filter(!Dropped) %>%
#   select(-c(Chain, Dropped))
# 
# new_output_df <- rbind(output_df[is.na(output_df$BICmcmc), ], kept_mixture_output)
# 
# 
# # additional_drops_base <- which(new_output_df$Distance > 200)
# # additional_drops_mvt <- which(new_output_df$Distance > 60 & new_output_df$Model %in% c("MVT", "MVN") & new_output_df$Scenario == "MVT")
# 
# new_output_df %>%
#   # filter(Distance < 75) %>%
#   ggplot(aes(x = Distance, y = Model)) +
#   geom_boxplot() +
#   facet_wrap(~Scenario)+#, scales = "free_x") +
#   labs(
#     title = "Simulation results",
#     subtitle = "Squared Euclidean distance"
#   )
# 
# ggsave("./Simulations/Output/SimModelDistance.png")
# 
# # additional_f1_drops <- which(new_output_df$F1 < 0.2)
# 
# new_output_df %>%
#   # filter(Distance < 75) %>%
#   ggplot(aes(x = F1, y = Model)) +
#   geom_boxplot() +
#   facet_wrap(~Scenario) +
#   labs(
#     title = "Simulation results",
#     subtitle = "F1 score"
#   )
# 
# ggsave("./Simulations/Output/SimModelF1.png")
# 
# # new_output_df %>%
# #   ggplot(aes(x = F1, y = Model)) +
# #   geom_boxplot() +
# #   geom_jitter(color="grey", size=0.7, alpha=0.5) +
# #   facet_wrap(~Scenario) +
# #   labs(
# #     title = "Simulation results: F1"
# #   )
# # 
# # ggsave("./Current_results_f1_jittered.png")
# 
# new_output_df %>%
#   ggplot(aes(x = Time, y = Model)) +
#   geom_boxplot() +
#   facet_wrap(~Scenario) +
#   labs(
#     title = "Simulation results",
#     subtitle = paste0("Time (", R, " MCMC iterations)")
#   )
# 
# ggsave("./Simulations/Output/SimModelTime.png")
# 

data.table::fwrite(output_df, paste0(my_wd, "Simulations/Output/SimModelPerformance.csv"))
