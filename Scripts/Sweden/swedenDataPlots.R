#
# Plots for the ELISA data from Dopico et al., 2021.
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
library(patchwork)

setMyTheme()
set.seed(1)

# Flag indicating if the paper branch of the batch mixture model package is used
# or the most recent
batchPackageIsCurrent <- FALSE

# Data lives in the repo:
# "https://github.com/chr1swallace/seroprevalence-paper/blob/master/adjusted-data.RData"
# Please download before proceeding if it is not already in the Data directory.
# Read in the ELISA data
data_file <- "./Data/Sweden/adjusted-data.RData"
load(data_file)

save_dir <- "./Analysis/Sweden/Outputs/Plots/"

# These samples are poorly behaved, drop
drop_sample_in_12 <- which((m$Sample.ID) %in% 1:2)
m <- m[-drop_sample_in_12, ]

# Drop the "patient 4" samples
patient_4 <- which(m$type == "Patient 4")
patient_4_data <- m[patient_4, ]
m <- m[-patient_4, ]


# Find the controls
negative_controls <- which(m$type == "Historical controls")
positive_controls <- which((m$type == "COVID") | (m$type == "Patient 4"))
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
    Batch = batch_vec,
    Fixed = fixed
  )

control_df <- plot_df %>%
  filter(Fixed == 1)

control_df$Controls <- "Controls"
plot_df$Controls <- "Full dataset"

p_data <- rbind(plot_df, control_df) %>%
  ggplot(aes(x = SPIKE, y = RBD, colour = factor(Label,
    labels = c("Seronegative", "Seropositive", "Unknown")
  ))) +
  geom_point(size = 0.7) +
  facet_grid(~Controls) +
  labs(
    # title = "ELISA data",
    x = "Log-transformed observed SPIKE measurement",
    y = "Log-transformed observed RBD measurement",
    colour = "Class"
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
R <- 50000
thin <- 100
burn <- 20000
eff_burn <- burn / thin

# Chain used based on complete likelihood trace plots
chain_used <- 3

# MCMC samples
mvt_samples <- readRDS(paste0(
  "./Analysis/Sweden/Outputs/Chains/sweden_mvt_chain_",
  chain_used,
  "_m_scale_1e-01_rho_11_theta_5.rds"
))
# mvn_samples <- readRDS("./Analysis/Outputs/sweden_mvn_chain_4.rds")



# === Parameter inference ======================================================

# Allocations
if (batchPackageIsCurrent) {
  mvt_samples$R <- R
  mvt_samples$thin <- thin

  mvt_prob <- calcAllocProb(mvt_samples, burn)
} else {
  mvt_prob <- calcAllocProb(mvt_samples$alloc, eff_burn)
}
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


p_inferred <- mvt_inferred_data %>%
  ggplot(aes(
    x = SPIKE,
    y = RBD,
    shape = factor(Fixed, label = c("False", "True")),
  )) +
  geom_point(
    aes(
      colour = factor(Label,
        levels = c(1, 2, 0),
        labels = c("Seronegative", "Seropositive", "Unknown")

        # labels = c("Seronegative", "Seropositive")
      ),
      alpha = Prob,
      size = factor(Fixed, label = c("False", "True"))
    ) # , size = 1
  ) +
  scale_alpha_continuous(range = c(0.4, 1.0)) +
  # geom_density_2d(aes(colour = Predicted)) +
  ggthemes::scale_color_colorblind() +
  labs(
    # title = "Elisa data - batch adjusted",
    # subtitle = "Predicted labels",
    # caption = "log transformed"
    x = "Log-transformed batch-corrected SPIKE OD",
    y = "Log-transformed batch-corrected RBD OD",
    colour = "Class",
    shape = "Control",
    alpha = "Probability of\nallocation",
    size = "Control"
  ) +
  scale_size_discrete(range = c(1.0, 1.6))

p_observed <- plot_df %>%
  ggplot(aes(
    x = SPIKE, y = RBD, colour = factor(Label,
      labels = c("Seronegative", "Seropositive", "Unknown")
    ),
    shape = factor(Fixed, label = c("False", "True")),
  )) +
  geom_point(aes(size = factor(Fixed, label = c("False", "True")))) +
  # facet_grid(~Controls) +
  labs(
    # title = "ELISA data",
    x = "Log-transformed observed SPIKE OD",
    y = "Log-transformed observed RBD OD",
    colour = "Class",
    shape = "Control",
    size = "Control"
  ) +
  ggthemes::scale_color_colorblind() +
  scale_size_discrete(range = c(1.0, 1.6))

p_patchwork <- p_observed / p_inferred +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect")
# theme(legend.position = "bottom")

ggsave(paste0(save_dir, "observedAndMVTInferredData.png"),
  plot = p_patchwork,
  height = 8,
  width = 7
)
