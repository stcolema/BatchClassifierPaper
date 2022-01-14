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

# Let's plot the original dataset
plot_df <- X %>%
  as_tibble() %>%
  add_column(
    Label = truth,
    Batch = batch_vec
  )

p_raw_data <- plot_df %>%
  # mutate(Plot_label = as.numeric(factor(Label)) * fixed) %>%
  ggplot(aes(
    x = SPIKE,
    y = RBD,
    colour = factor(Label,
      labels = c("Seronegative", "Seropositive", "Unknown")
    )
  )) +
  geom_point(size = 0.7) +
  labs(
    title = "ELISA data",
    colour = "Class"
  ) +
  ggthemes::scale_color_colorblind()

ggsave(paste0(save_dir, "ELISA_observed.png"),
  plot = p_raw_data,
  width = 5,
  height = 3.75
)


p_controls <- plot_df[fixed == 1, ] %>%
  # mutate(Plot_label = as.numeric(factor(Label, levels = c("Historical controls", "COVID")))) %>%
  ggplot(aes(
    x = SPIKE,
    y = RBD,
    colour = factor(Label,
      labels = c("Seronegative", "Seropositive")
    )
  )) +
  geom_point(size = 0.7) +
  labs(
    title = "ELISA data",
    subtitle = "Controls",
    colour = "Class"
  ) +
  ggthemes::scale_color_colorblind()


ggsave(paste0(save_dir, "ELISA_controls.png"),
  plot = p_controls,
  width = 5,
  height = 3.75
)

plot_df %>%
  group_by(Batch, Label) %>%
  summarise(count = n()) %>%
  knitr::kable()

# Model inputs
K_max <- 2
B <- length(unique(batch_vec))
P <- ncol(X)

# The class proportions in the obsered labels (no longer used)
# class_weights <- table(fixed) / N

# Our initial labelling is 1 for known controls, 2 for known COVIDs
initial_labels <- rep(1, N)
initial_labels[negative_controls] <- 1
initial_labels[positive_controls] <- 2

# We then sample with equal probability the class for the remaining items
initial_labels[-c(negative_controls, positive_controls)] <- sample(1:2,
  size = N - N_positive - N_negative,
  replace = TRUE,
  prob = c(1, 1)
)

# Convert the batch vector to integers
batch_vec <- batch_vec %>%
  factor() %>%
  as.numeric()

# Hyperparameter for the group weights prior
concentration <- rep(1, K_max)

# === Prior hyperparameters ====================================================

# This sets the scale of the batch mean effect
m_scale <- c(1, 0.1, 0.01)
n_m_scales <- length(m_scale)

# Choice of rho and theta to have a common expectation
rho <- c(2, 20, 10) + 1
theta <- c(1, 10, 5)

# Quick check that the rhos and thetas align
n_theta <- length(theta)
n_rho <- length(rho)
rhos_match_thetas <- n_rho == n_theta
if (!rhos_match_thetas) {
  stop("Number of rhos to try and number of thetas to try not matching.")
}

hyper_params <- data.frame(
  m_scale = rep(m_scale, n_theta),
  rho = rep(rho, each = n_m_scales),
  theta = rep(theta, each = n_m_scales)
)

# How many combinations
n_hyperparam_combinations <- nrow(hyper_params)

# === MCMC =====================================================================

# MCMC parameters
R <- 50000
thin <- 100
burn <- 20000
eff_burn <- burn / thin

# Number of chains to run
n_chains <- 5L
mvt_samples <- mvn_samples <- vector("list", n_chains * n_hyperparam_combinations)

# What chain are we currently saving
chain_count <- 0

for (j in seq(1, n_hyperparam_combinations)) {

  # Hyperparameters for the prior distribution of the batch effects
  .m_scale <- hyper_params$m_scale[j]
  .rho <- hyper_params$rho[j]
  .theta <- hyper_params$theta[j]

  # Generate a Markov chain for each model and record acceptance rates separately
  for (i in 1:n_chains) {
    set.seed(i)

    chain_count <- chain_count + 1

    # Modelling
    t0 <- Sys.time()
    mvt_samples[[chain_count]] <- .mvt <- batchSemiSupervisedMixtureModel(
      X,
      R,
      thin,
      initial_labels,
      fixed,
      batch_vec,
      type = "MVT",
      alpha = 1,
      mu_proposal_window = 0.12**2,
      cov_proposal_window = 200,
      m_proposal_window = 0.12**2,
      S_proposal_window = 75,
      t_df_proposal_window = 25,
      m_scale = .m_scale,
      rho = .rho,
      theta = .theta
    )

    t1 <- Sys.time()

    # Modelling MVN
    set.seed(i)
    t2 <- Sys.time()
    mvn_samples[[chain_count]] <- .mvn <- batchSemiSupervisedMixtureModel(
      X,
      R,
      thin,
      initial_labels,
      fixed,
      batch_vec,
      "MVN",
      alpha = 1,
      mu_proposal_window = 0.12**2,
      cov_proposal_window = 350,
      m_proposal_window = 0.12**2,
      S_proposal_window = 125,
      m_scale = .m_scale,
      rho = .rho,
      theta = .theta
    )

    t3 <- Sys.time()

    # Model likelihood
    likelihood_df_local <- data.frame(
      "mvt_complete" = .mvt$complete_likelihood,
      "mvt_observed" = .mvt$observed_likelihood,
      "mvn_complete" = .mvn$complete_likelihood,
      "mvn_observed" = .mvn$observed_likelihood,
      "mvt_bic" = .mvt$BIC,
      "mvn_bic" = .mvn$BIC,
      "Chain" = i,
      "Iteration" = (1:(R / thin)) * thin,
      "m_scale" = .m_scale,
      "rho" = .rho,
      "theta" = .theta
    ) %>%
      pivot_longer(c(
        mvt_complete,
        mvt_observed,
        mvn_complete,
        mvn_observed,
        mvt_bic,
        mvn_bic
      ),
      names_sep = "_",
      names_to = c("model", "param")
      )

    # Acceptance rates
    mvt_cluster_acceptance_df <- data.frame(
      Mu = .mvt$mu_acceptance_rate,
      Cov = .mvt$cov_acceptance_rate,
      DF = .mvt$t_df_acceptance_rate,
      K = seq(1, length(.mvt$mu_acceptance_rate)),
      Model = "MVT",
      Chain = i,
      m_scale = .m_scale,
      rho = .rho,
      theta = .theta
    )

    mvt_batch_acceptance_df <- data.frame(
      m = .mvt$m_acceptance_rate,
      S = .mvt$S_acceptance_rate,
      B = seq(1, length(.mvt$m_acceptance_rate)),
      Model = "MVT",
      Chain = i,
      m_scale = .m_scale,
      rho = .rho,
      theta = .theta
    )

    mvn_cluster_acceptance_df <- data.frame(
      Mu = .mvn$mu_acceptance_rate,
      Cov = .mvn$cov_acceptance_rate,
      DF = NA,
      K = seq(1, length(.mvn$mu_acceptance_rate)),
      Model = "MVN",
      Chain = i,
      m_scale = .m_scale,
      rho = .rho,
      theta = .theta
    )

    mvn_batch_acceptance_df <- data.frame(
      m = .mvn$m_acceptance_rate,
      S = .mvn$S_acceptance_rate,
      B = seq(1, length(.mvn$m_acceptance_rate)),
      Model = "MVN",
      Chain = i,
      m_scale = .m_scale,
      rho = .rho,
      theta = .theta
    )


    if (i == 1 & j == 1) {
      batch_acceptance_df <- rbind(mvn_batch_acceptance_df, mvt_batch_acceptance_df)
      cluster_acceptance_df <- rbind(mvn_cluster_acceptance_df, mvt_cluster_acceptance_df)
      likelihood_df <- likelihood_df_local
    } else {
      cluster_acceptance_df <- rbind(
        cluster_acceptance_df,
        mvn_cluster_acceptance_df,
        mvt_cluster_acceptance_df
      )
      batch_acceptance_df <- rbind(
        batch_acceptance_df,
        mvn_batch_acceptance_df,
        mvt_batch_acceptance_df
      )

      likelihood_df <- rbind(likelihood_df, likelihood_df_local)
    }

    # For saving numbers with a decimal place
    m_scale_sci_notation <- formatC(.m_scale, format = "e", digits = 0)

    # Save the samples for the current chain
    saveRDS(.mvt, file = paste0(
      my_wd,
      "/Analysis/Sweden/Outputs/sweden_mvt_chain_",
      i,
      "_m_scale_",
      m_scale_sci_notation,
      "_rho_",
      .rho,
      "_theta_",
      .theta,
      ".rds"
    ))

    saveRDS(.mvn, file = paste0(
      my_wd,
      "/Analysis/Sweden/Outputs/sweden_mvn_chain_",
      i,
      "_m_scale_",
      m_scale_sci_notation,
      "_rho_",
      .rho,
      "_theta_",
      .theta,
      ".rds"
    ))
  }
}

# === Write outputs ============================================================

# Save the various data frames that we use to check model behaviour
write.csv(cluster_acceptance_df,
  paste0(my_wd, "/Analysis/Sweden/ModelChecks/groupParamAcceptance.csv"),
  row.names = FALSE
)

write.csv(batch_acceptance_df,
  paste0(my_wd, "/Analysis/Sweden/ModelChecks/batchParamAcceptance.csv"),
  row.names = FALSE
)

write.csv(likelihood_df,
  paste0(my_wd, "/Analysis/Sweden/ModelChecks/likelihoods.csv"),
  row.names = FALSE
)
