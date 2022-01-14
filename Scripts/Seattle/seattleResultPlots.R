#
# Plots for the ELISA data from
# Dingens, A.S., Crawford, K.H.D., Adler, A. et al. (2020).
# We plot the
#   * the observed and the inferred datasets
#   * seroprevalence over time (here batch)
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
library(patchwork)

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
R <- 25000
thin <- 50
burn <- 15000
eff_burn <- burn / thin

# Number of chains to run
n_chains <- 10L

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

  # 95% credible interval
  mvn_ci_seropositive <- mvn_samples$alloc[, 2, ] %>% 
    apply(1, quantile, probs = c(0.025, 0.975)) %>% 
    t()
  
  mvt_ci_seropositive <- mvt_samples$alloc[, 2, ] %>% 
    apply(1, quantile, probs = c(0.025, 0.975)) %>% 
    t()
  
  mvn_seropositive <- data.frame("Estimate" = mvn_prob[-controls, 2], 
    "CI_0.025" = mvn_ci_seropositive[-controls, 1],
    "CI_0.975" = mvn_ci_seropositive[-controls, 2],
    Batch = batch_vec[-controls],
    "Model" = "MVN",
    "Chain" = chain_used,
    m_scale = .m_scale,
    rho = .rho,
    theta = .theta
  )
  
  mvt_seropositive <- data.frame("Estimate" = mvt_prob[-controls, 2], 
    "CI_0.025" = mvt_ci_seropositive[-controls, 1],
    "CI_0.975" = mvt_ci_seropositive[-controls, 2],
    Batch = batch_vec[-controls],
    "Model" = "MVT",
    "Chain" = chain_used,
    m_scale = .m_scale,
    rho = .rho,
    theta = .theta
  )
  
  seropositive_df <- rbind(mvn_seropositive, mvt_seropositive)
  
  if (i == 1) {
    seroprevalence_prob_df <- seropositive_df
  } else {
    seroprevalence_prob_df <- rbind(seroprevalence_prob_df, seropositive_df)
  }
  
  
  
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

  if (i == 1) {
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
    Label = if_else(Hit_5SD == "True", 1, 0),
    Model = "Dingens_est"
  ) %>%
  group_by(Batch, Model) %>%
  summarise(Seroprevalence = 100 * sum(Label) / n()) %>% 
  mutate("97.5%" = Seroprevalence,
             "02.5%" = Seroprevalence,
         Name = "Dingens")

my_sero_df <- seroprevalence_prob_df  %>%
  group_by(Batch, Model, m_scale, rho, theta) %>%
  summarise(
    # Seroprevalence = 100 * sum(Estimate > 0.5) / n(),
    Seroprevalence = 100 * sum(Estimate) / n(),
    "97.5%" = 100 * sum(CI_0.975) / n(),
    "02.5%" = 100 * sum(CI_0.025) / n(),
    Name = paste0(Model, "_", m_scale, "_", rho, "_", theta)) %>% 
  filter(Model != "MVN")

p_sero_w_ci <- rbind(dingens_pred, my_sero_df) %>% 
  ggplot(aes(x = Batch, 
    group = Name)) +
  # geom_line(aes(y = Seroprevalence, 
  #               colour = factor(Model, labels = c("Dingens et al.", "MVT"))),
  #           lty = 2) +
  geom_point(aes(y = Seroprevalence, 
                 colour = factor(Model, labels = c("Dingens et al.", "MVT")))) +
  geom_errorbar(data = my_sero_df,
                mapping = aes(ymin=`02.5%`, ymax=`97.5%`), width=.1, colour = "#00BFC4") +
  # geom_line(aes(y = `97.5%`), lty = 2) +
  # geom_line(aes(y = `02.5%`), lty = 2) +
  labs(
    y = "Seroprevalence (%)",
    colour = "Model:"
  )
# +
  # theme(legend.position = "bottom")

dingens_pred <- non_control_data[, c(7, 12)] %>%
  mutate(
    Batch = as.numeric(factor(screen_batch)),
    Label = if_else(Hit_5SD == "True", 1, 0)
  ) %>%
  group_by(Batch) %>%
  summarise(Dingen_est = 100 * sum(Label) / n())

seroprevalence_df <- week_prob_df %>%
  pivot_longer(c(MVT_prob, MVT_pred, MVN_prob, MVN_pred),
    names_sep = "_", #
    names_to = c("Model", ".value")
  ) %>%
  group_by(Batch, Model, m_scale, rho, theta) %>%
  summarise(Seroprevalence = 100 * sum(prob) / n())

p_sero <- seroprevalence_df %>%
  pivot_wider(
    names_from = c(Model, m_scale, rho, theta),
    values_from = Seroprevalence
  ) %>%
  add_column(Dingens = dingens_pred$Dingen_est) %>%
  pivot_longer(-c(Batch), values_to = "Seroprevalence") %>%
  mutate(Model = ifelse(grepl("MVN_", name), "MVN", ifelse(grepl("MVT_", name), "MVT", "Dingens et al."))) %>%
  filter(Model != "MVN") %>%
  ggplot(aes(x = Batch, y = Seroprevalence, colour = Model)) +
  geom_point() +
  geom_line() +
  labs(
    y = "Seroprevalence (%)",
    colour = "Model:"
  ) +
  theme(legend.position = "bottom")

ggsave(paste0(save_dir, "seroprevalence.png"),
  height = 5,
  width = 5,
  plot = p_sero
)

# === Example chain ============================================================

# Chain used based on complete likelihood trace plots
chain_used <- 3

# MCMC samples
mvt_samples <- readRDS(paste0(
  "./Analysis/Seattle/Outputs/Chains/seattle_mvt_chain_",
  chain_used,
  "_m_scale_1e-01_rho_11_theta_5.rds"
))

# === Parameter inference ======================================================

# Allocations
mvt_prob <- calcAllocProb(mvt_samples$alloc, eff_burn)
mvt_pred <- predictClass(mvt_prob)

# Inferred datasets
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
    "Model" = "MVT"
  )


p_inferred <- mvt_inferred_data %>%
  ggplot(aes(x = Batch, y = RBD)) +
  geom_jitter(aes(
    colour = factor(Label,
      labels = c("Seronegative", "Seropositive")
    ),
    alpha = Prob,
    shape = factor(Fixed, label = c("False", "True")),
    size = factor(Fixed, label = c("False", "True"))
  ),
  # size = 1
  ) +
  scale_alpha_continuous(range = c(0.4, 1.0)) +
  # geom_density_2d(aes(colour = Predicted)) +
  ggthemes::scale_color_colorblind() +
  labs(
    # title = "Elisa data - batch adjusted",
    # subtitle = "Predicted labels",
    # caption = "log transformed"
    x = "Batch",
    y = "Log-corrected RBD OD",
    colour = "Group",
    shape = "Control",
    size = "Control",
    alpha = "Probability of\nallocation"
  ) +
  scale_size_discrete(range = c(1.0, 1.6))

p_observed <- plot_df %>%
  mutate(Decision_boundary = rbd_screen_results$Cutoff_5SD) %>% 
  ggplot(aes(
    x = Batch, y = RBD, colour = factor(Label,
      labels = c("Seronegative", "Seropositive", "Unknown")
    ),
    shape = factor(Fixed, label = c("False", "True")),
    size = factor(Fixed, label = c("False", "True"))
  )) +
  geom_jitter(size = 1) +
  labs(
    # title = "ELISA data",
    x = "Batch",
    y = "Log RBD OD",
    colour = "Group",
    shape = "Control"
  ) +
  ggthemes::scale_color_colorblind() +
  scale_size_discrete(range = c(1.0, 1.6))

p_patchwork <- (p_observed / p_inferred) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect")
# theme(legend.position = "bottom")

ggsave(paste0(save_dir, "observedAndMVTInferredData.png"),
  plot = p_patchwork,
  height = 8,
  width = 7
)

patch1 <- (p_observed / p_inferred) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect")



layout = "
AB
CC
"

p_seattle <- wrap_plots(p_sero, p_observed, p_inferred, design =layout) + #widths = c(2, 1), ncol = 2) + # (p_observed / p_inferred) | p_sero +
  plot_annotation(tag_levels = "A")
# +
  # plot_layout(guides = "collect")

ggsave(paste0(save_dir, "seattleResults.png"),
  plot = p_seattle,
  height = 7,
  width = 8
)

p_seattle_w_ci <- wrap_plots(p_observed, 
                             p_inferred, 
                             p_sero_w_ci, 
                             design =layout) + #widths = c(2, 1), ncol = 2) + # (p_observed / p_inferred) | p_sero +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

ggsave(paste0(save_dir, "seattleResultsWCI.png"),
       plot = p_seattle_w_ci,
       height = 7,
       width = 8
)

# === Compare decision boundaries ==============================================

mvt_inferred_data_2 <- mvt_inferred_data %>% 
  mutate(Decision_boundary = log(rbd_screen_results$Cutoff_5SD),
         Label = factor(Label, 
                        levels = c(1, 2, 3),
                        labels = c("Seronegative", "Seropositive", "Unknown")))

p_new <- plot_df %>%
  select(-Controls) %>% 
  mutate(Decision_boundary = log(rbd_screen_results$Cutoff_5SD), 
         Prob = 1, 
         Model = NA, 
         Type = "Observed") %>% 
  rbind(mvt_inferred_data_2) %>%
  ggplot(aes(x = Type, y = RBD)) +
  geom_jitter(aes(
    shape = factor(Fixed, label = c("False", "True")),
    colour = factor(Label,
                    labels = c("Seronegative", "Seropositive", "Unknown")
    ),
    alpha = Prob,
    size = Fixed
  )
  ) +
  scale_alpha_continuous(range = c(0.4, 1.0)) +
  # geom_density_2d(aes(colour = Predicted)) +
  ggthemes::scale_color_colorblind() +
  labs(
    # title = "Elisa data - batch adjusted",
    # subtitle = "Predicted labels",
    # caption = "log transformed"
    x = "Type",
    y = "Log RBD OD",
    colour = "Group",
    shape = "Control",
    alpha = "Probability of\nallocation",
    yintercept = "5SD decision boundary"
  ) +
  geom_hline(aes(yintercept = Decision_boundary), colour = "red", lty = 2) +
  facet_wrap(~Batch, labeller = label_both) +
  scale_size_discrete(range = c(1.0, 1.6))

p_new

dingens <- plot_df %>% 
  select(-Controls) %>% 
  mutate(Decision_boundary = log(rbd_screen_results$Cutoff_5SD),
         Label = factor(ifelse(rbd_screen_results$Hit_5SD == "True", 2, 1),
                        levels = c(1, 2, 3),
                        labels = c("Seronegative", "Seropositive", "Unknown")),                        
         Prob = 1,
         Model = "5SD",
         Type = "Dingens et al."
  )

plot_df %>%
  select(-Controls) %>% 
  mutate(Decision_boundary = log(rbd_screen_results$Cutoff_5SD), 
         Prob = 1, 
         Model = NA, 
         Type = "Observed") %>% 
  rbind(mvt_inferred_data_2, dingens) %>%
  ggplot(aes(x = Type, y = RBD)) +
  geom_jitter(aes(
    shape = factor(Fixed, label = c("False", "True")),
    colour = factor(Label,
                    labels = c("Seronegative", "Seropositive", "Unknown")
    ),
    alpha = Prob,
    size = Fixed
  )
  ) +
  scale_alpha_continuous(range = c(0.4, 1.0)) +
  # geom_density_2d(aes(colour = Predicted)) +
  ggthemes::scale_color_colorblind() +
  labs(
    # title = "Elisa data - batch adjusted",
    # subtitle = "Predicted labels",
    # caption = "log transformed"
    x = "Type",
    y = "Log RBD OD",
    colour = "Group",
    shape = "Control",
    alpha = "Probability of\nallocation",
    yintercept = "5SD decision boundary"
  ) +
  geom_hline(aes(yintercept = Decision_boundary), colour = "red", lty = 2) +
  facet_wrap(~Batch, labeller = label_both) +
  scale_size_discrete(range = c(1.0, 1.6))
