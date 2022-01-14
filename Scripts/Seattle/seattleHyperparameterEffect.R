#
# Plots for the ELISA data from
# Dingens, A.S., Crawford, K.H.D., Adler, A. et al. (2020).
# We plot the
#   * the observed data and the controls
#   * complete log-likelihood to check within/across chain convergence
#   * the inferred and observed data
#   * allocation probability proportion over time
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
library(stringr)
library(tibble)

setMyTheme()
set.seed(1)

# Data lives in the repo:
# "https://github.com/chr1swallace/seroprevalence-paper/blob/master/adjusted-data.RData"
# Please download before proceeding if it is not already in the Data directory.
# Read in the ELISA data
rbd_screen_results <- read.csv("./Data/Seattle/Original/RBD_screen_results.csv")

# ELISA controls (not seronegative controls)
dropped_samples <- which(rbd_screen_results$category == "no sera")
rbd_screen_results <- rbd_screen_results[-dropped_samples, ]

save_dir <- "./Analysis/Seattle/Outputs/Plots/"

# Log transform the data and move to matrix format
X <- log(matrix(rbd_screen_results$Screen_OD450, ncol = 1))
colnames(X) <- "RBD"

# Number of samples
N <- nrow(X)

# Batch labels
batch_vec <- rbd_screen_results$screen_batch

# Types present in data
types <- factor(rbd_screen_results$category)


# Find the controls
negative_controls <- which(types == "pre-2020 sera pool")
positive_controls <- which(types == "CR3022 antibody")
controls <- c(negative_controls, positive_controls)

N_positive <- length(positive_controls)
N_negative <- length(negative_controls)
N_controls <- length(controls)

initial_labels <- rep(0, N)

# The non-control data is used when estimating seroprevalence
non_control_data <- X[-controls, , drop = F]

# Find which observations have a known label
fixed <- rep(0, N)
fixed[c(positive_controls, negative_controls)] <- 1

initial_labels[fixed] <- 1

# Our batch variable
batches_present <- rbd_screen_results$screen_batch %>%
  unique()

batch_vec <- rbd_screen_results$screen_batch %>%
  as.factor() %>%
  as.numeric()


observed_labels <- rep("Unknown", N)
observed_labels[positive_controls] <- "Seropositive"
observed_labels[negative_controls] <- "Seronegative"

# Let's plot the original dataset
plot_df <- X %>%
  as_tibble() %>%
  add_column(
    Label = observed_labels,
    Batch = batch_vec,
    Fixed = fixed
  )

control_df <- plot_df %>%
  filter(Fixed == 1)

control_df$Controls <- "Controls"
plot_df$Controls <- "Full dataset"

# p_data <- rbind(plot_df, control_df) %>%
#   ggplot(aes(x = Batch, y = RBD, colour = Label)) +
#   geom_jitter(size = 0.7) +
#   facet_grid(~Controls) +
#   labs(
#     # title = "ELISA data",
#     colour = "Group"
#   ) +
#   ggthemes::scale_color_colorblind()

p_data <- plot_df %>%
  ggplot(aes(x = Batch, y = RBD, colour = Label)) +
  geom_jitter(size = 0.7) +
  labs(
    # title = "ELISA data",
    colour = "Group"
  ) +
  ggthemes::scale_color_colorblind()

ggsave(paste0(save_dir, "ELISA_observed_and_controls.png"),
  plot = p_data,
  width = 8,
  height = 5
)

#
# p_controls <- plot_df[fixed == 1, ] %>%
#   # mutate(Plot_label = as.numeric(factor(Label, levels = c("Historical controls", "COVID")))) %>%
#   ggplot(aes(
#     x = SPIKE,
#     y = RBD,
#     colour = factor(Label,
#                     labels = c("Seronegative", "Seropositive")
#     )
#   )) +
#   geom_point(size = 0.7) +
#   labs(
#     title = "ELISA data",
#     subtitle = "Controls",
#     colour = "Class"
#   ) +
#   ggthemes::scale_color_colorblind()
#
#
# ggsave(paste0(save_dir, "ELISA_controls.png"),
#        plot = p_controls,
#        width = 5,
#        height = 3.75
# )

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
R <- 15000
thin <- 25
burn <- 5000
eff_burn <- burn / thin

# === Extract data =============================================================

chain_used <- 3

my_list <- list(1)

# Tibble to hold results
chain_tib <- hyper_params %>%
  as_tibble() %>%
  mutate(
    Chain = chain_used,
    R = R,
    thin = thin,
    burn = burn,
    MVT_samples = my_list,
    MVN_samples = my_list,
    MVN_prob = my_list,
    MVN_pred = my_list,
    MVT_prob = my_list,
    MVT_pred = my_list
  )

chain_dir <- "./Analysis/Seattle/Outputs/Chains/"

for (i in seq.int(n_hyperparam_combinations)) {
  .m_scale <- chain_tib$m_scale[i]
  .rho <- chain_tib$rho[i]
  .theta <- chain_tib$theta[i]

  m_scale_sci_notation <- formatC(.m_scale, format = "e", digits = 0)

  f_mvt <- paste0(
    chain_dir,
    "seattle_",
    "mvt",
    "_chain_",
    chain_used,
    "_m_scale_",
    m_scale_sci_notation,
    "_rho_",
    .rho,
    "_theta_",
    .theta,
    ".rds"
  )

  mvt_samples <- readRDS(f_mvt)

  f_mvn <- paste0(
    chain_dir,
    "seattle_",
    "mvn",
    "_chain_",
    chain_used,
    "_m_scale_",
    m_scale_sci_notation,
    "_rho_",
    .rho,
    "_theta_",
    .theta,
    ".rds"
  )

  mvn_samples <- readRDS(f_mvn)


  # Allocations
  mvn_prob <- calcAllocProb(mvn_samples$alloc, eff_burn)
  mvn_pred <- predictClass(mvn_prob)

  mvt_prob <- calcAllocProb(mvt_samples$alloc, eff_burn)
  mvt_pred <- predictClass(mvt_prob)

  # Inferred datasets
  mvn_inferred_data <- rowMeans(mvn_samples$batch_corrected_data[, , -c(1:eff_burn), drop = F], dims = 2) %>%
    as_tibble(.name_repair = "minimal") %>%
    set_colnames(colnames(X)) %>%
    add_column(
      "Label" = factor(mvn_pred),
      "Prob" = apply(mvn_prob, 1, function(x) {
        x[which.max(x)]
      }),
      "Batch" = factor(batch_vec),
      "Type" = "Inferred",
      "Fixed" = factor(fixed),
      "Model" = "MVN",
      "Chain" = chain_used,
      m_scale = .m_scale,
      rho = .rho,
      theta = .theta
    )
  
  mvt_inferred_data <- rowMeans(mvt_samples$batch_corrected_data[, , -c(1:eff_burn), drop = F], dims = 2) %>%
    as_tibble(.name_repair = "minimal") %>%
    set_colnames(colnames(X)) %>%
    add_column(
      "Label" = factor(mvt_pred),
      "Prob" = apply(mvt_prob, 1, function(x) {
        x[which.max(x)]
      }),
      "Batch" = factor(batch_vec),
      "Type" = "Inferred",
      "Fixed" = factor(fixed),
      "Model" = "MVT",
      "Chain" = chain_used,
      m_scale = .m_scale,
      rho = .rho,
      theta = .theta
    )
  
  .inf_df <- rbind(mvn_inferred_data, mvt_inferred_data)
  
  if(i == 1) {
    inferred_data <- .inf_df
  } else {
    inferred_data <- rbind(inferred_data, .inf_df)
  }
  
  
  # Find the allocations across time. here batch and time are one and the same.

  # Associate the allocation probabilities with the week
  .week_df <- data.frame(
    Batch = batch_vec[-controls],
    m_scale = .m_scale,
    rho = .rho,
    theta = .theta,
    chain = chain_used,
    MVN_prob = mvn_prob[-controls, 2],
    MVT_prob = mvt_prob[-controls, 2],
    MVN_pred = mvn_pred[-controls] - 1,
    MVT_pred = mvt_pred[-controls] - 1
  )

  if (i == 1) {
    week_prob_df <- .week_df
  } else {
    week_prob_df <- rbind(week_prob_df, .week_df)
  }

  # Sampled parameters
  mvn_sampled_data <- samplesToDF(mvn_samples, "MVN", R, thin,
    keep_allocation = FALSE
  )
  mvn_sampled_data$Model <- "MVN"
  mvn_sampled_data$m_scale <- .m_scale
  mvn_sampled_data$rho <- .rho
  mvn_sampled_data$theta <- .theta
  mvn_sampled_data$chain <- chain_used
  mvn_sampled_data$t_df_1 <- mvn_sampled_data$t_df_2 <- NA

  mvt_sampled_data <- samplesToDF(mvt_samples, "MVT", R, thin,
    keep_allocation = FALSE
  )

  mvt_sampled_data$Model <- "MVT"
  mvt_sampled_data$m_scale <- .m_scale
  mvt_sampled_data$rho <- .rho
  mvt_sampled_data$theta <- .theta
  mvt_sampled_data$chain <- chain_used

  # Bind the sampled parameters together and move the degrees of freedom into
  # the df
  .samples_df <- rbind(mvn_sampled_data, mvt_sampled_data)[, c(1:31, 40:41, 32:39)]

  if (i == 1) {
    sampled_params <- .samples_df
  } else {
    sampled_params <- rbind(sampled_params, .samples_df)
  }

  # Save everything to the tibble
  chain_tib$MVT_samples[[i]] <- mvt_samples
  chain_tib$MVN_samples[[i]] <- mvn_samples

  chain_tib$MVN_prob[[i]] <- mvn_prob
  chain_tib$MVN_pred[[i]] <- mvn_pred

  chain_tib$MVT_prob[[i]] <- mvt_prob
  chain_tib$MVT_pred[[i]] <- mvt_pred
}

# === Plotting =================================================================

## == Seroprevalence ===========================================================

non_control_data <- rbd_screen_results[-controls, ]


dingens_pred <- non_control_data[, c(7, 12)] %>%
  mutate(
    Batch = as.numeric(factor(screen_batch)),
    Label = if_else(Hit_5SD == "True", 1, 0)
  ) %>%
  group_by(Batch) %>%
  summarise(Dingen_est = sum(Label))

seroprevalence_df <- week_prob_df %>%
  pivot_longer(c(MVT_prob, MVT_pred, MVN_prob, MVN_pred),
               names_sep = "_", #
               names_to = c("Model", ".value")
  ) %>%
  group_by(Batch, Model, m_scale, rho, theta) %>%
  summarise(Num_seropositive = sum(prob),
            Count = n())


p_alloc_time_incl_dingens <- seroprevalence_df %>%
  pivot_wider(
    names_from = c(Model, m_scale, rho, theta),
    values_from = Num_seropositive
  ) %>%
  add_column(Dingens = dingens_pred$Dingen_est) %>%
  pivot_longer(-c(Batch, Count)) %>% 
  group_by(Batch, name) %>% 
  summarise(Seroprevalence = value / Count) %>%
  ggplot(aes(x = Batch, y = 100 * Seroprevalence, colour = name)) +
  geom_point() +
  labs(y = "Seroprevalence (%)")

# p_alloc_time_incl_dingens <- seroprevalence_df %>%
#   select(-Count) %>% 
#   pivot_wider(
#     names_from = c(Model, m_scale, rho, theta),
#     values_from = Num_seropositive
#   ) %>%
#   add_column(Dingens = dingens_pred$Dingen_est) %>%
#   pivot_longer(-Batch) %>%
#   ggplot(aes(x = Batch, y = value, colour = name)) +
#   geom_point()

p_alloc_time <- week_prob_df %>%
  pivot_longer(c(MVT_prob, MVT_pred, MVN_prob, MVN_pred),
    names_sep = "_", #
    names_to = c("Model", ".value")
  ) %>%
  group_by(Batch, Model, m_scale, rho, theta) %>%
  summarise(Seroprevalence = 100 * sum(prob) / n()) %>%
  ggplot(aes(x = Batch, y = Seroprevalence, colour = Model)) +
  geom_point() +
  geom_line(aes(colour = interaction(Model, m_scale, rho, theta))) +
  # geom_smooth(se = F, method = "gam") +
  labs(
    title = "Seattle data",
    y = "Seroprevalence (%)",
    # y = "Number of serorpositive allocations",
    colour = "Method"
  )

ggsave(paste0(save_dir, "Positive_allocations_across_time_different_hypers.png"),
  height = 5,
  width = 5,
  plot = p_alloc_time
)

ggsave(paste0(save_dir, "Seroprevalence_different_hypers_incl_dingens.png"),
  height = 5,
  width = 5,
  plot = p_alloc_time_incl_dingens
)

## == Sampled parameters =======================================================

sampled_params %>%
  select(pi_1, pi_2, Iteration, Model, m_scale, rho, theta, chain) %>%
  pivot_longer(-c(Iteration, Model, m_scale, rho, theta, chain)) %>%
  ggplot(aes(
    x = value,
    # fill = name,
    fill = interaction(name, Model, m_scale, rho, theta)
  )) +
  geom_density(alpha = 0.6)

p_weights <- sampled_params %>%
  filter(Iteration > burn) %>%
  select(pi_1, pi_2, Iteration, Model, m_scale, rho, theta, chain) %>%
  pivot_longer(-c(Iteration, Model, m_scale, rho, theta, chain)) %>%
  ggplot(aes(
    x = value,
    # fill = name,
    fill = interaction(name, Model)
  )) +
  geom_density(alpha = 0.6) +
  facet_grid(m_scale + rho + theta ~ Model)

p_group_means <- sampled_params %>%
  filter(Iteration > burn) %>%
  select(Mu_11, Mu_21, Iteration, Model, m_scale, rho, theta, chain) %>%
  pivot_longer(-c(Iteration, Model, m_scale, rho, theta, chain)) %>%
  ggplot(aes(
    x = value,
    fill = name
  )) +
  # fill = interaction(name, Model, m_scale, rho, theta))) +
  geom_density(alpha = 0.6) +
  facet_grid(m_scale + rho + theta ~ Model)

p_batch_means <- sampled_params %>%
  filter(Iteration > burn) %>%
  select(
    m_11,
    m_21,
    m_31,
    m_41,
    Iteration,
    Model,
    m_scale,
    rho,
    theta,
    chain
  ) %>%
  pivot_longer(-c(Iteration, Model, m_scale, rho, theta, chain)) %>%
  ggplot(aes(
    x = value,
    fill = name
  )) +
  # fill = interaction(name, Model, m_scale, rho, theta))) +
  geom_density(alpha = 0.6) +
  facet_grid(m_scale + rho + theta ~ Model)


p_batch_scales <- sampled_params %>%
  filter(Iteration > burn) %>%
  select(
    S_11,
    S_21,
    S_31,
    S_41,
    Iteration,
    Model,
    m_scale,
    rho,
    theta,
    chain
  ) %>%
  pivot_longer(-c(Iteration, Model, m_scale, rho, theta, chain)) %>%
  ggplot(aes(
    x = value,
    fill = name
  )) +
  # fill = interaction(name, Model, m_scale, rho, theta))) +
  geom_density(alpha = 0.6) +
  facet_grid(m_scale + rho + theta ~ Model)

p_combined_means <- sampled_params %>%
  filter(Iteration > burn) %>%
  select(
    Mu_11.m_11,
    Mu_21.m_11,
    Mu_11.m_21,
    Mu_21.m_21,
    Mu_11.m_31,
    Mu_21.m_31,
    Mu_11.m_41,
    Mu_21.m_41,
    Iteration,
    Model,
    m_scale,
    rho,
    theta,
    chain
  ) %>%
  pivot_longer(-c(Iteration, Model, m_scale, rho, theta, chain)) %>%
  ggplot(aes(
    x = value,
    fill = name
  )) +
  # fill = interaction(name, Model, m_scale, rho, theta))) +
  geom_density(alpha = 0.6) +
  facet_grid(m_scale + rho + theta ~ Model)

p_combined_cov <- sampled_params %>%
  filter(Iteration > burn) %>%
  select(
    Sigma_111.S_11,
    Sigma_211.S_11,
    Sigma_111.S_21,
    Sigma_211.S_21,
    Sigma_111.S_31,
    Sigma_211.S_31,
    Sigma_111.S_41,
    Sigma_211.S_41,
    Iteration,
    Model,
    m_scale,
    rho,
    theta,
    chain
  ) %>%
  pivot_longer(-c(Iteration, Model, m_scale, rho, theta, chain)) %>%
  ggplot(aes(
    x = value,
    fill = name
  )) +
  # fill = interaction(name, Model, m_scale, rho, theta))) +
  geom_density(alpha = 0.6) +
  facet_grid(m_scale + rho + theta ~ Model)

# Save plots
ggsave(paste0(save_dir, "SampledParameters/group_weights.png"),
  height = 6,
  width = 7,
  plot = p_weights
)

ggsave(paste0(save_dir, "SampledParameters/group_means.png"),
  height = 6,
  width = 7,
  plot = p_group_means
)

ggsave(paste0(save_dir, "SampledParameters/batch_means.png"),
  height = 6,
  width = 7,
  plot = p_batch_means
)

ggsave(paste0(save_dir, "SampledParameters/batch_scales.png"),
  height = 6,
  width = 7,
  plot = p_batch_scales
)

ggsave(paste0(save_dir, "SampledParameters/combined_means.png"),
  height = 6,
  width = 10,
  plot = p_combined_means
)

ggsave(paste0(save_dir, "SampledParameters/combined_cov.png"),
  height = 6,
  width = 10,
  plot = p_combined_cov
)

# === Inferred dataset =========================================================


observed_labels <- initial_labels
observed_labels[fixed == 0] <- 3

m_scale_used <- 0.1
rho_used <- 11
theta_used <- 5

mvt_pred_used <- inferred_data$Label[
  inferred_data$m_scale == m_scale_used &
    inferred_data$rho == rho_used &
    inferred_data$theta == theta_used &
    inferred_data$Model == "MVT"
]

observed_data <- X %>%
  as_tibble() %>%
  add_column(
    "Label" = mvt_pred_used,
    "Prob" = 1,
    "Batch" = factor(batch_vec),
    "Type" = "Observed",
    "Fixed" = factor(fixed),
    "Model" = "Observed",
    Chain = NA,
    m_scale = NA, 
    rho = NA,
    theta = NA
  )

p_inferred_all <- inferred_data %>% 
  ggplot(aes(x = Batch, y = RBD, colour = Label, alpha = Prob)) +
  geom_jitter() + 
  facet_grid(m_scale + rho + theta ~ Model) +
  ggthemes::scale_color_colorblind() +
  labs(
    title = "Seattle: inferred datasets"
  )

p_inferred_and_observed <- inferred_data %>% 
  filter(m_scale == m_scale_used, rho == rho_used, theta == theta_used) %>% 
  rbind(observed_data) %>% 
  ggplot(aes(x = Batch, y = RBD, colour = Label, alpha = Prob)) +
  geom_jitter() + 
  facet_grid(~ Model) +
  ggthemes::scale_color_colorblind() +
  labs(
    title = "Seattle: inferred and observed data"
  )


ggsave(paste0(save_dir, "inferred_datasets.png"),
       plot = p_inferred_all,
       height = 6,
       width = 10
)

ggsave(paste0(save_dir, "inferred_data_and_observed.png"),
       plot = p_inferred_and_observed,
       height = 6,
       width = 10
)

# === Parameter inference ======================================================

# Allocations
mvt_prob <- calcAllocProb(mvt_samples$alloc, eff_burn)
mvt_pred <- predictClass(mvt_prob)

# Inferred datasets
mvt_inferred_data <- rowMeans(mvt_samples$batch_corrected_data[, , -c(1:eff_burn)], dims = 2) %>%
  as_tibble(.name_repair = "minimal") %>%
  set_colnames(colnames(X)) %>%
  add_column(
    "Label" = factor(mvt_pred),
    "Prob" = apply(mvt_prob, 1, function(x) {
      x[which.max(x)]
    }),
    "Batch" = factor(batch_vec),
    "Type" = "Inferred",
    "Fixed" = factor(fixed),
    "Model" = "MVT"
  )


mvt_inferred_data %>%
  ggplot(aes(x = SPIKE, y = RBD)) +
  geom_point(aes(shape = Fixed, colour = Label, alpha = Prob)) +
  scale_alpha_continuous(range = c(0.4, 1.0)) +
  # geom_density_2d(aes(colour = Predicted)) +
  scale_color_viridis_d() +
  labs(
    title = "Elisa data - batch adjusted",
    subtitle = "Predicted labels",
    caption = "log transformed"
  )

# Check the correlation for the SPIKE and RBD measurements in the original and
# the inferred datasets
cor(X)
cor(rowMeans(mvt_samples$batch_corrected_data[, , -c(1:eff_burn)], dims = 2))

# # The MVT inferred dataset (uncertain points highlighted)
# p_mvt_inferred <- mvt_inferred_data %>%
#   ggplot(aes(x = SPIKE, y = RBD)) +
#   geom_point(aes(
#     shape = Fixed,
#     colour = factor(Label, labels = c("Seronegative", "Seropositive")),
#     alpha = 1 - Prob
#   ),
#   size = 1
#   ) +
#   scale_alpha_continuous(range = c(0.4, 1.0)) +
#   # facet_wrap(~Model, ncol = 1) +
#   # geom_density_2d(aes(colour = Predicted)) +
#   ggthemes::scale_color_colorblind() +
#   labs(
#     title = "ELISA data",
#     subtitle = "Inferred dataset and labels for MVT model",
#     caption = "log transformed",
#     colour = "Class",
#     alpha = "1 - Probability"
#   )
#
# ggsave("./Data/MVT_inferred_dataset.png",
#   plot = p_mvt_inferred,
#   width = 5,
#   height = 3.75
# )
#
# p_inferred <- rbind(mvt_inferred_data, mvn_inferred_data) %>%
#   mutate(Label = factor(Label, labels = c("Uninfected", "COVID"))) %>%
#   ggplot(aes(x = SPIKE, y = RBD)) +
#   geom_point(aes(shape = Fixed, colour = Label, alpha = Prob), size = 0.7) +
#   scale_alpha_continuous(range = c(0.4, 1.0)) +
#   facet_wrap(~Model, ncol = 1) +
#   # geom_density_2d(aes(colour = Predicted)) +
#   ggthemes::scale_color_colorblind() +
#   labs(
#     title = "Elisa data",
#     subtitle = "Inferred datasets and labels",
#     caption = "log transformed"
#   )
#
observed_labels <- initial_labels
observed_labels[fixed == 0] <- 3

observed_data <- X %>%
  as_tibble() %>%
  add_column(
    "Label" = observed_labels,
    "Prob" = 1,
    "Batch" = factor(batch_vec),
    "Type" = "Observed",
    "Fixed" = factor(fixed),
    "Model" = "Observed"
  )

p_res <- rbind(observed_data, mvt_inferred_data) %>%
  mutate(Label = factor(Label, levels = c(1, 2, 3), labels = c("Seronegative", "Seropositive", "Unknown"))) %>%
  ggplot(aes(x = SPIKE, y = RBD)) +
  geom_point(aes(shape = Fixed, colour = Label, alpha = Prob), size = 0.4) +
  scale_alpha_continuous(range = c(0.4, 1.0)) +
  facet_wrap(~Model, ncol = 2) +
  # geom_density_2d(aes(colour = Predicted)) +
  ggthemes::scale_color_colorblind() +
  labs(
    # title = "Elisa data",
    # subtitle = "Observed and inferred datasets and labels",
    caption = "log transformed",
    colour = "Group",
    alpha = "Probability"
  )

p_res <- rbind(observed_data, mvt_inferred_data) %>%
  mutate(Label = factor(Label, levels = c(1, 2, 3), labels = c("Seronegative", "Seropositive", "Unknown"))) %>%
  ggplot(aes(x = SPIKE, y = RBD)) +
  geom_point(aes(shape = Fixed, colour = Label, alpha = 1 - Prob), size = 0.4) +
  scale_alpha_continuous(range = c(0.4, 1.0)) +
  facet_wrap(~Model, ncol = 1) +
  # geom_density_2d(aes(colour = Predicted)) +
  ggthemes::scale_color_colorblind() +
  labs(
    title = "Elisa data",
    subtitle = "Observed and inferred datasets and labels",
    caption = "log transformed",
    colour = "Class",
    alpha = "1 - Probability"
  )


rbind(observed_data, mvt_inferred_data) %>%
  mutate(Label = factor(Label, levels = c(1, 2, 3), labels = c("Seronegative", "Seropositive", "Unknown"))) %>%
  ggplot(aes(x = SPIKE, y = Model)) +
  geom_point(aes(shape = Fixed, colour = Label, alpha = Prob), size = 0.7) +
  scale_alpha_continuous(range = c(0.4, 1.0)) +
  # facet_grid(~Batch) +
  # geom_density_2d(aes(colour = Predicted)) +
  ggthemes::scale_color_colorblind() +
  labs(
    title = "Elisa data",
    subtitle = "Observed and inferred datasets and labels",
    caption = "log transformed",
    colour = "Class",
    alpha = "Probability"
  )

rbind(observed_data, mvt_inferred_data) %>%
  pivot_longer(c(SPIKE, RBD), names_to = "Antibody", values_to = "Measurement") %>%
  mutate(Label = factor(Label, levels = c(1, 2, 3), labels = c("Seronegative", "Seropositive", "Unknown"))) %>%
  ggplot(aes(x = Measurement, y = Model)) +
  geom_point(aes(shape = Fixed, colour = Label, alpha = Prob), size = 0.7) +
  # geom_density(aes(colour = Label)) +
  scale_alpha_continuous(range = c(0.4, 1.0)) +
  facet_grid(~Antibody) +
  # geom_density_2d(aes(colour = Predicted)) +
  ggthemes::scale_color_colorblind() +
  labs(
    title = "Elisa data",
    subtitle = "Observed and inferred datasets and labels",
    caption = "log transformed",
    colour = "Class",
    alpha = "Probability"
  )

observed_with_mvt_label <- observed_data
observed_with_mvt_label$Label <- mvt_inferred_data$Label

p_density <- rbind(observed_with_mvt_label, mvt_inferred_data) %>%
  pivot_longer(c(SPIKE, RBD), names_to = "Antibody", values_to = "Measurement") %>%
  mutate(Label = factor(Label, levels = c(1, 2, 3), labels = c("Seronegative", "Seropositive", "Unknown"))) %>%
  ggplot(aes(x = Measurement)) +
  # geom_point(aes(shape = Fixed, colour = Label, alpha = Prob), size = 0.7) +
  geom_density(aes(fill = Label), alpha = 0.4) +
  scale_alpha_continuous(range = c(0.4, 1.0)) +
  facet_grid(Model ~ Antibody) +
  # geom_density_2d(aes(colour = Predicted)) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  labs(
    title = "Sweden data",
    subtitle = "Observed and inferred densities for antibody responses",
    colour = "Group"
  ) +
  geom_hline(data = data.frame())

ggsave(paste0(save_dir, "ELISA_observed_and_inferred_densities.png"),
  plot = p_density,
  width = 8,
  height = 7
)

ggsave(paste0(save_dir, "ELISA_observed_and_inferred_data.png"),
  plot = p_res,
  width = 8,
  height = 7
)


ggsave(paste0(save_dir, "ELISA_observed_and_inferred_data.png"),
  plot = p_res,
  width = 6,
  height = 4
)


library(ggExtra)

p1 <- observed_data %>%
  mutate(Label = factor(Label, levels = c(1, 2, 3), labels = c("Seronegative", "Seropositive", "Unknown"))) %>%
  ggplot(aes(x = SPIKE, y = RBD)) +
  geom_point(aes(shape = Fixed, colour = Label, alpha = 1 - Prob), size = 0.4) +
  scale_alpha_continuous(range = c(0.4, 1.0)) +
  ggthemes::scale_color_colorblind() +
  labs(
    colour = "Group"
  )

p2 <- mvt_inferred_data %>%
  mutate(Label = factor(Label, levels = c(1, 2, 3), labels = c("Seronegative", "Seropositive", "Unknown"))) %>%
  ggplot(aes(x = SPIKE, y = RBD)) +
  geom_point(aes(shape = Fixed, colour = Label, alpha = 1 - Prob), size = 0.4) +
  scale_alpha_continuous(range = c(0.4, 1.0)) +
  ggthemes::scale_color_colorblind() +
  labs(
    colour = "Group"
  )

# install.packages("patchwork")
library(patchwork)

p_obs <- ggMarginal(p1, type = "histogram")
p_inf <- ggMarginal(p2, type = "histogram")

# === Allocations across time ==================================================

neg_control_df <- m[negative_controls, c(6, 7)]

neg_summary_stats <- neg_control_df %>%
  summarise(
    Mean_spike = mean(SPIKE),
    Mean_rbd = mean(RBD),
    SD_spike = sd(SPIKE),
    SD_rbd = sd(RBD)
  )

cutoffs <- data.frame(
  Threshold = c("3SD", "5SD"),
  Spike = c(
    neg_summary_stats$Mean_spike + 3 * neg_summary_stats$SD_spike,
    neg_summary_stats$Mean_spike + 5 * neg_summary_stats$SD_spike
  ),
  RBD = c(
    neg_summary_stats$Mean_rbd + 3 * neg_summary_stats$SD_rbd,
    neg_summary_stats$Mean_rbd + 5 * neg_summary_stats$SD_rbd
  )
)

trad_inference <- observed_data %>%
  mutate(
    SPIKE_3SD_pos = exp(SPIKE) > cutoffs$Spike[1],
    SPIKE_5SD_pos = exp(SPIKE) > cutoffs$Spike[2],
    RBD_3SD_pos = exp(RBD) > cutoffs$RBD[1],
    RBD_5SD_pos = exp(RBD) > cutoffs$RBD[2]
  )

trad_inference %>%
  ggplot(aes(x = SPIKE, y = RBD)) +
  geom_point() +
  geom_hline(yintercept = c(cutoffs$RBD[1], cutoffs$RBD[2]), lty = 2) +
  geom_vline(xintercept = c(cutoffs$Spike[1], cutoffs$Spike[2]), lty = 2)

p_res_w_thresholds <- p_res +
  geom_hline(yintercept = c(cutoffs$RBD[1], cutoffs$RBD[2]), lty = 2) +
  geom_vline(xintercept = c(cutoffs$Spike[1], cutoffs$Spike[2]), lty = 2)

ggsave(paste0(save_dir, "ELISA_observed_and_inferred_data_sd_thres.png"),
  plot = p_res_w_thresholds,
  width = 8,
  height = 7
)

non_control_data <- m[-controls, ]

# The non-control data has the week of the sample collection included in the
# sample ID
weeks_str <- non_control_data$Sample.ID

# Extract the week
weeks <- weeks_str %>%
  str_match_all("Wk(\\d{1,2})") %>%
  do.call(rbind, .) %>%
  `[`(, 2) %>%
  as.numeric()

weeks_present <- unique(weeks)

# Associate the allocation probabilities with the week
week_prob_df <- data.frame(
  Week = weeks,
  Seronegative_prob = mvt_prob[-controls, 1],
  Seropositive_prob = mvt_prob[-controls, 2],
  Allocation = mvt_pred[-controls],
  Spike_3SD = trad_inference$SPIKE_3SD_pos[-controls],
  Spike_5SD = trad_inference$SPIKE_5SD_pos[-controls],
  RBD_3SD = trad_inference$RBD_3SD_pos[-controls],
  RBD_5SD = trad_inference$RBD_5SD_pos[-controls]
)


# Associate the allocation probabilities with the week
week_prob_df <- data.frame(
  Week = weeks,
  Seronegative_prob = mvt_prob[-controls, 1],
  Seropositive_prob = mvt_prob[-controls, 2],
  Allocation = mvt_pred[-controls],
  Spike_3SD = trad_inference$SPIKE_3SD_pos[-controls],
  Spike_5SD = trad_inference$SPIKE_5SD_pos[-controls],
  RBD_3SD = trad_inference$RBD_3SD_pos[-controls],
  RBD_5SD = trad_inference$RBD_5SD_pos[-controls]
)

week_prob_df$MVT_prob <- mvt_prob[-controls, 2]

(colSums(week_prob_df[, c(5:9)]) / (N - N_controls)) %>% round(3)

p_alloc_time_all <- week_prob_df %>%
  pivot_longer(c(Spike_3SD, Spike_5SD, RBD_3SD, RBD_5SD, MVT_prob)) %>%
  group_by(Week, name) %>%
  summarise(Num_seropositive = sum(value)) %>%
  ggplot(aes(x = Week, y = Num_seropositive, colour = name)) +
  geom_point() +
  # geom_smooth(se = F, method = "gam") +
  labs(
    title = "Seropositive allocations across time",
    y = "Number of serorpositive allocations",
    colour = "Method"
  ) +
  ggthemes::scale_color_colorblind()

ggsave(paste0(save_dir, "Positive_allocations_across_time_different_options.png"),
  height = 5,
  width = 5,
  plot = p_alloc_time_all
)


# Plot the number of seropositive allocations as a function of time
p_alloc_time <- week_prob_df %>%
  group_by(Week) %>%
  summarise(
    Seronegative = sum(Seronegative_prob), # sum(Allocation == 1),
    Seropositive = sum(Seropositive_prob) # sum(Allocation == 2)
  ) %>%
  ggplot(aes(x = Week, y = Seropositive)) +
  geom_point() +
  labs(
    title = "Seropositive allocations across time",
    y = "Number of serorpositive allocations"
  )

ggsave(paste0(save_dir, "Positive_allocations_across_time.png"),
  height = 5,
  width = 5,
  plot = p_alloc_time
)


# ggsave(paste0("./Data/sweden_full_chain_", chain_used, "_mvn_v_mvt_inferred.png"),
#   plot = p_all_datasets,
#   height = 7,
#   width = 6
# )
#
#
# rbind(observed_data, mvt_inferred_data, mvn_inferred_data) %>%
#   mutate(Label = factor(Label, levels = c(1, 2, 3), labels = c("Seronegative", "Seropositive", "Unknown"))) %>%
#   ggplot(aes(x = SPIKE, y = RBD)) +
#   geom_point(aes(shape = Fixed, colour = Label, alpha = Prob), size = 0.4) +
#   scale_alpha_continuous(range = c(0.4, 1.0)) +
#   facet_wrap(~Model, ncol = 1) +
#   # geom_density_2d(aes(colour = Predicted)) +
#   ggthemes::scale_color_colorblind() +
#   labs(
#     title = "Elisa data",
#     subtitle = "Observed and inferred datasets and labels",
#     caption = "log transformed",
#     colour = "Class",
#     alpha = "Probability"
#   )
#
# ggsave(paste0("./Data/sweden_full_chain_", chain_used, "_mvn_v_mvt_inferred_norm_prob.png"), height = 7, width = 6)
#
#
# # Degrees of freedom for the MVT
# for (i in 1:n_chains) {
#   .x <- data.frame(mvt_samples[[i]]$t_df[-c(1:eff_burn), ]) %>%
#     set_colnames(c("DF_1", "DF_2"))
#   .x$Chain <- i
#   .x$Model <- "MVT"
#
#   if (i == 1) {
#     df_df <- .x
#   } else {
#     df_df <- rbind(df_df, .x)
#   }
# }
#
# p_df <- df_df %>%
#   pivot_longer(c(DF_1, DF_2), names_to = "Parameter", values_to = "Value") %>%
#   ggplot(aes(x = Value)) +
#   geom_histogram() +
#   facet_grid(Chain ~ Parameter, scales = "free") +
#   labs(title = "Posterior distribution for degrees of freedom")
#
# ggsave("./Sweden_sampled_degrees_of_freedom.png",
#   plot = p_df,
#   height = 6,
#   width = 4
# )
#
#
# samples_used <- which(mvt_prob[, 1] < 0.6 & mvt_prob[, 1] > 0.4 |
#   mvn_prob[, 1] < 0.6 & mvn_prob[, 1] > 0.4)
#
# row_order <- findOrder(mvt_prob[samples_used, ])
#
#
# for (i in 1:n_chains) {
#   .prob <- calcAllocProb(mvn_samples[[i]]$alloc[samples_used, , ], eff_burn)
#
#   .x <- .prob[row_order, ] %>%
#     as_tibble(rownames = "Y") %>%
#     set_colnames(c("Y", 1, 2)) %>%
#     pivot_longer(-Y, names_to = "X")
#
#   .x$Chain <- i
#   .x$Model <- "MVN"
#
#   if (i == 1) {
#     mvn_alloc_df <- .x
#   } else {
#     mvn_alloc_df <- rbind(mvn_alloc_df, .x)
#   }
# }
#
# for (i in 1:n_chains) {
#   .prob <- calcAllocProb(mvt_samples[[i]]$alloc[samples_used, , ], eff_burn)
#
#   .x <- .prob[row_order, ] %>%
#     as_tibble(rownames = "Y") %>%
#     set_colnames(c("Y", 1, 2)) %>%
#     pivot_longer(-Y, names_to = "X")
#
#   .x$Chain <- i
#   .x$Model <- "MVT"
#
#   if (i == 1) {
#     mvt_alloc_df <- .x
#   } else {
#     mvt_alloc_df <- rbind(mvt_alloc_df, .x)
#   }
# }
#
# mvt_alloc_df$Y <- as.numeric(mvt_alloc_df$Y)
# mvn_alloc_df$Y <- as.numeric(mvn_alloc_df$Y)
#
# mvn_alloc_df %>%
#   ggplot(aes(x = as.factor(X), y = as.numeric(Y), fill = value)) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "#146EB4") +
#   facet_grid(~Chain)
#
# mvt_alloc_df %>%
#   ggplot(aes(x = as.factor(X), y = as.numeric(Y), fill = value)) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "#146EB4") +
#   facet_grid(~Chain)
#
# p_alloc <- rbind(mvn_alloc_df, mvt_alloc_df) %>%
#   ggplot(aes(x = as.factor(X), y = as.numeric(Y), fill = value)) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "#146EB4") +
#   facet_grid(Model ~ Chain, labeller = label_both) +
#   theme(
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank(),
#     panel.grid = element_blank(),
#     axis.title.y = element_blank(),
#     # axis.title.x = element_text(size = 10.5),
#     plot.title = element_text(face = "bold"),
#     # plot.subtitle = element_text(size = 14),
#     # strip.text.x = element_text(size = 10.5),
#     # legend.text = element_text(size = 10.5)
#   ) +
#   labs(
#     title = "Allocation probabilities",
#     subtitle = "Samples with allocation probability < 0.6 in chain 3 of either model",
#     x = "Class"
#   )
#
#
# samples_used <- 1:N # which(mvt_prob[,1] < 0.6 | mvn_prob[,1] < 0.6)
#
# row_order <- findOrder(mvt_prob[samples_used, ])
#
#
# for (i in 1:n_chains) {
#   .prob <- calcAllocProb(mvn_samples[[i]]$alloc[samples_used, , ], eff_burn)
#
#   .x <- .prob[row_order, ] %>%
#     as_tibble(rownames = "Y") %>%
#     set_colnames(c("Y", 1, 2)) %>%
#     pivot_longer(-Y, names_to = "X")
#
#   .x$Chain <- i
#   .x$Model <- "MVN"
#
#   if (i == 1) {
#     mvn_alloc_df <- .x
#   } else {
#     mvn_alloc_df <- rbind(mvn_alloc_df, .x)
#   }
# }
#
# for (i in 1:n_chains) {
#   .prob <- calcAllocProb(mvt_samples[[i]]$alloc[samples_used, , ], eff_burn)
#
#   .x <- .prob[row_order, ] %>%
#     as_tibble(rownames = "Y") %>%
#     set_colnames(c("Y", 1, 2)) %>%
#     pivot_longer(-Y, names_to = "X")
#
#   .x$Chain <- i
#   .x$Model <- "MVT"
#
#   if (i == 1) {
#     mvt_alloc_df <- .x
#   } else {
#     mvt_alloc_df <- rbind(mvt_alloc_df, .x)
#   }
# }
#
# mvt_alloc_df$Y <- as.numeric(mvt_alloc_df$Y)
# mvn_alloc_df$Y <- as.numeric(mvn_alloc_df$Y)
#
# mvn_alloc_df %>%
#   ggplot(aes(x = as.factor(X), y = as.numeric(Y), fill = value)) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "#146EB4") +
#   facet_grid(~Chain)
#
# mvt_alloc_df %>%
#   ggplot(aes(x = as.factor(X), y = as.numeric(Y), fill = value)) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "#146EB4") +
#   facet_grid(~Chain)
#
# rbind(mvn_alloc_df, mvt_alloc_df) %>%
#   ggplot(aes(x = as.factor(X), y = as.numeric(Y), fill = value)) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "#146EB4") +
#   facet_grid(Model ~ Chain, labeller = label_both) +
#   theme(
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank(),
#     panel.grid = element_blank(),
#     axis.title.y = element_blank(),
#     # axis.title.x = element_text(size = 10.5),
#     plot.title = element_text(face = "bold"),
#     # plot.subtitle = element_text(size = 14),
#     # strip.text.x = element_text(size = 10.5),
#     # legend.text = element_text(size = 10.5)
#   ) +
#   labs(
#     title = "Allocation probabilities",
#     # subtitle = "Samples with allocation probability < 0.6 in chain 3 of either model",
#     x = "Class"
#   )
#
# ggsave("./Alloc_probs_across_chains_models_sweden.png")
# ggsave("./Subset_alloc_probs_across_chains_models_sweden.png", plot = p_alloc)
#
#
# # observed_data[samples_used, ],
# rbind(
#   mvt_inferred_data[samples_used, ],
#   mvn_inferred_data[samples_used, ]
# ) %>%
#   mutate(Label = factor(Label,
#     levels = c(1, 2, 3),
#     labels = c("Seronegative", "Seropositive", "Unknown")
#   )) %>%
#   ggplot(aes(x = SPIKE, y = RBD)) +
#   geom_point(aes(shape = Fixed, colour = Label, alpha = 1 - Prob), size = 0.4) +
#   scale_alpha_continuous(range = c(0.4, 1.0)) +
#   facet_wrap(~Model, ncol = 1) +
#   # geom_density_2d(aes(colour = Predicted)) +
#   ggthemes::scale_color_colorblind() +
#   labs(
#     title = "Elisa data",
#     subtitle = "Observed and inferred datasets and labels",
#     caption = "log transformed"
#   )
#
#
# # === Inferred group density ===================================================
#
# mvn_cov <- mvn_samples$covariance[, , -c(1:eff_burn)] %>%
#   rowMeans(dims = 2L)
#
# mvn_mu <- mvn_samples$means[, , -c(1:eff_burn)] %>%
#   rowMeans(dims = 2L)
#
# mvt_cov <- mvt_samples$covariance[, , -c(1:eff_burn)] %>%
#   rowMeans(dims = 2L)
#
# mvt_mu <- mvt_samples$means[, , -c(1:eff_burn)] %>%
#   rowMeans(dims = 2L)
#
# mvt_S <- mvt_samples$batch_scale[, , -c(1:eff_burn)] %>%
#   rowMeans(dims = 2L)
#
# mvt_m <- mvt_samples$batch_shift[, , -c(1:eff_burn)] %>%
#   rowMeans(dims = 2L)
#
# mvt_df <- mvt_samples$t_df[-c(1:eff_burn), ] %>% colMeans()
#
# library(mvtnorm)
#
# n_sampled <- 1e5
# mvt_g1 <- rmvt(n_sampled, sigma = mvt_cov[1:2, 1:2], df = mvt_df[1], delta = mvt_mu[, 1])
# mvt_g2 <- rmvt(n_sampled, sigma = mvt_cov[1:2, 3:4], df = mvt_df[2], delta = mvt_mu[, 2])
#
# mvn_g1 <- rmvnorm(n_sampled, sigma = mvn_cov[1:2, 1:2], mean = mvn_mu[, 1])
# mvn_g2 <- rmvnorm(n_sampled, sigma = mvn_cov[1:2, 3:4], mean = mvn_mu[, 2])
#
# mvt_g1 <- data.frame(mvt_g1) %>% add_column("Class" = "Sero-negative")
# mvt_g2 <- data.frame(mvt_g2) %>% add_column("Class" = "Sero-positive")
# mvt_densities <- rbind(mvt_g1, mvt_g2) %>% add_column("Model" = "MVT")
#
# mvn_g1 <- data.frame(mvn_g1) %>% add_column("Class" = "Sero-negative")
# mvn_g2 <- data.frame(mvn_g2) %>% add_column("Class" = "Sero-positive")
# mvn_densities <- rbind(mvn_g1, mvn_g2) %>% add_column("Model" = "MVN")
#
# density_df <- rbind(mvn_densities, mvt_densities)
# colnames(density_df) <- c("SPIKE", "RBD", "Class", "Model")
#
#
# p_density <- rbind(mvn_densities, mvt_densities) %>%
#   ggplot(aes(x = X1, y = X2, colour = Class)) +
#   geom_density_2d() +
#   facet_grid(~Model) +
#   labs(x = "SPIKE", y = "RBD")
#
# p_density
#
# inferred_data <- rbind(mvt_inferred_data, mvn_inferred_data)
#
# inferred_data$Class <- factor(inferred_data$Label, levels = c(1, 2), labels = c("Sero-negative", "Sero-positive"))
#
# p_data <- inferred_data %>%
#   ggplot(aes(x = SPIKE, y = RBD, colour = Class)) +
#   geom_point(alpha = 0.1) +
#   scale_alpha_continuous(range = c(0.1, 1)) +
#   facet_wrap(~Model) +
#   ggthemes::scale_color_colorblind() +
#   labs(title = "Density from rmv(t/norm) + geom_density_2d")
#
#
# p_density_flavour_1 <- p_data +
#   geom_density_2d(data = density_df, mapping = aes(color = Class)) +
#   facet_wrap(~Model)
#
# ggsave("./Sweden_inferred_densities_v1.png",
#   plot = p_density_flavour_1,
#   height = 6,
#   width = 10
# )
#
# #
# p_density +
#   geom_point(data = inferred_data, mapping = aes(x = SPIKE, y = RBD, colour = Class)) +
#   facet_grid(~Model)
#
#
# confidenceEllipse <- function(mu = c(0, 0), Sigma = matrix(c(1, 0, 0, 1), 2, 2), confidenceLevel = 0.95) {
#   radius <- sqrt(2 * stats::qf(confidenceLevel, 2, Inf))
#   chol_decomp <- chol(Sigma)
#   angles <- (0:100) * 2 * pi / 100
#   unit.circle <- cbind(cos(angles), sin(angles))
#   ellipse <- t(mu + radius * t(unit.circle %*% chol_decomp))
#   colnames(ellipse) <- c("X1", "X2")
#   as.data.frame(ellipse)
# }
#
#
# mvn_ellipse_seronegative <- confidenceEllipse(
#   mu = mvn_mu[, 1],
#   Sigma = mvn_cov[1:2, 1:2],
#   confidenceLevel = 0.95
# )
#
# mvn_ellipse_seronegative$Class <- "Sero-negative"
#
# mvn_ellipse_seropositive <- confidenceEllipse(
#   mu = mvn_mu[, 2],
#   Sigma = mvn_cov[1:2, 3:4],
#   confidenceLevel = 0.95
# )
#
# mvn_ellipse_seropositive$Class <- "Sero-positive"
#
# mvn_ellipse <- rbind(mvn_ellipse_seronegative, mvn_ellipse_seropositive)
# mvn_ellipse$Model <- "MVN"
#
# mvt_ellipse_seronegative <- confidenceEllipse(
#   mu = mvt_mu[, 1],
#   Sigma = mvt_cov[1:2, 1:2],
#   confidenceLevel = 0.95
# )
#
# mvt_ellipse_seronegative$Class <- "Sero-negative"
#
# mvt_ellipse_seropositive <- confidenceEllipse(
#   mu = mvt_mu[, 2],
#   Sigma = mvt_cov[1:2, 3:4],
#   confidenceLevel = 0.95
# )
#
# mvt_ellipse_seropositive$Class <- "Sero-positive"
#
# mvt_ellipse <- rbind(mvt_ellipse_seronegative, mvt_ellipse_seropositive)
# mvt_ellipse$Model <- "MVT"
#
# p_density_flavour_2 <- p_data +
#   geom_path(mvn_ellipse, mapping = aes(x = X1, y = X2, colour = Class)) +
#   geom_path(mvt_ellipse, mapping = aes(x = X1, y = X2, colour = Class)) +
#   ggthemes::scale_color_colorblind() +
#   labs(title = "Density using Paul's confidence ellipse function")
#
# ggsave("./Sweden_inferred_densities_v2.png",
#   plot = p_density_flavour_2,
#   height = 6,
#   width = 5
# )
