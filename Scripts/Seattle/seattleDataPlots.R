#
# Plots for the observed and inferred data from Dopico et al., 2021.
# We plot the 
#   * the observed data and the controls
#   * the inferred data
# This forms figure 4 in the paper
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


# === MCMC =====================================================================

# MCMC parameters
R <- 15000
thin <- 25
burn <- 5000
eff_burn <- burn / thin

# Chain used based on complete likelihood trace plots
chain_used <- 3

# MCMC samples
mvt_samples <- readRDS(paste0("./Analysis/Seattle/Outputs/Chains/seattle_mvt_chain_",
                              chain_used, 
                              "_m_scale_1e-01_rho_11_theta_5.rds")
)
# mvn_samples <- readRDS("./Analysis/Outputs/sweden_mvn_chain_4.rds")



# # === Parameter inference ======================================================

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
  ggplot(aes(x = Batch, y = RBD))  +
  geom_jitter(aes(shape = factor(Fixed, label = c("False", "True")),
                 colour = factor(Label,
                                 labels = c("Seronegative", "Seropositive")
                 ),
                 alpha = Prob),
             size = 1) +
  scale_alpha_continuous(range = c(0.4, 1.0)) +
  # geom_density_2d(aes(colour = Predicted)) +
  ggthemes::scale_color_colorblind() +  
  labs(
    # title = "Elisa data - batch adjusted",
    # subtitle = "Predicted labels",
    # caption = "log transformed"
    x = "Batch",
    y = "Log-transformed batch-corrected RBD OD",
    colour = "Group",
    shape = "Control",
    alpha = "Probability of\nallocation"
  )

p_observed <- plot_df %>% 
  ggplot(aes(x = Batch, y = RBD, colour = factor(Label,
                                                 labels = c("Seronegative", "Seropositive", "Unknown")
  ), 
  shape = factor(Fixed, label = c("False", "True")),
  )
  ) +
  geom_jitter(size = 1) +
  # facet_grid(~Controls) +
  labs(
    # title = "ELISA data",
    x = "Batch",
    y = "Log-transformed observed RBD OD",
    colour = "Group",
    shape = "Control"
  ) +
  ggthemes::scale_color_colorblind()

p_patchwork <- p_observed / p_inferred +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect") 
# theme(legend.position = "bottom")

ggsave(paste0(save_dir, "observedAndMVTInferredData.png"),
       plot = p_patchwork,
       height = 8,
       width = 7
)
