# visualisation.R
# Using the output of ``modelling.R`` create the plots to interpret the
# performance of the methods on the simulations and check convergence in the
# Bayesian methods.
#
# === Setup ====================================================================

# install and call various packages
cran_packages <- c(
  "tidyr",
  "ggplot2",
  "dplyr",
  "magrittr",
  "devtools",
  "tibble",
  "data.table"
)

install.packages(setdiff(cran_packages, rownames(installed.packages())))

library(tidyr)
library(ggplot2)
library(dplyr)
library(magrittr)
library(tibble)
library(coda)
library(data.table)
library(patchwork)
library(mdiHelpR)
library(BatchMixtureModel)

# Set the gglot2 theme
setMyTheme()

# Set the random seed
set.seed(1)

# Save the generated files to
save_path <- "./Simulations/ELISA_like/Output/"

# The MCMC number of iterations, etc.
R <- 50000
thin <- 100
burn <- 20000
eff_burn <- burn / thin

# The number of chains/sims/scenarios
n_sim <- 10
n_chains <- 10

scenario_string <- "ELISA-like"

# === Model comparison =========================================================

# Read in the data regarding performance of each model
output_df <- fread("./Simulations/ELISA_like/Output/SimModelPerformance.csv")

# Niver labels for the models
output_df$Model <- factor(output_df$Model,
                          levels = unique(output_df$Model),
                          labels = c("MVN", "MVT", "RF", "SVM", "LR", "LR (batch corrected)")
)

# The fruit of the various convergence steps, the data frame indicating which
# chains have converged.
keep_chain_info <- fread("./Simulations/ELISA_like/Output/Convergence/GewekeDF.csv", select = 4:7) %>%
  distinct()

keep_chain_info <- fread("./Simulations/ELISA_like/Output/Convergence/GewekeAndLikelihood.csv")

dir.create("./Simulations/ELISA_like/Output/ModelPerformance/")

# Separate out the mixture and ML results
mixture_output <- output_df[!is.na(output_df$BICmcmc), ]
ml_output <- output_df[is.na(output_df$BICmcmc), ]

# Add information about which chains have behaved
keep_chain_info$Seed <- keep_chain_info$Chain
mixture_output_new <- mixture_output %>%
  left_join(keep_chain_info, by = c("Model", "Simulation", "Seed")) %>%
  select(-Chain)

# We keep all the ML results
ml_output$Kept <- T

# The used output data frame
new_output_df <- rbind(mixture_output_new, ml_output)

# Plot of the distance between the probability allocation matrix and the truth
p_distance <- new_output_df %>%
  filter(Kept) %>%
  ggplot(aes(x = Distance, y = Model)) +
  geom_boxplot(fill = "gold") +
  labs(
    x = "Squared Euclidean distance"
    # title = "ELISA-like: Simulation results",
    # subtitle = "Squared Euclidean distance"
  )

ggsave("./Simulations/ELISA_like/Output/ModelPerformance/SimModelDistance.png",
       plot = p_distance,
       height = 4.5,
       width = 7
)

# Plot of the F1 score
p_f1 <- new_output_df %>%
  filter(Kept) %>%
  ggplot(aes(x = F1, y = Model)) +
  geom_boxplot(fill = "gold") +
  labs(
    x = "F1 score"
    # title = "ELISA-like: Simulation results",
    # subtitle = "F1 score"
  )

ggsave("./Simulations/ELISA_like/Output/ModelPerformance/SimModelF1.png",
       plot = p_f1,
       height = 4.5,
       width = 7
)

# Plot of the F1 score
p_f1_violin <- new_output_df %>%
  filter(Kept) %>%
  ggplot(aes(x = F1, y = Model)) +
  geom_violin(fill = "gold") +
  labs(
    x = "F1 score"
    # title = "ELISA-like: Simulation results",
    # subtitle = "F1 score"
  )

ggsave("./Simulations/ELISA_like/Output/ModelPerformance/SimModelF1Violin.png",
       plot = p_f1_violin,
       height = 4.5,
       width = 7
)

# Plot of the estimated seroprevalence
p_sero <- new_output_df %>%
  filter(Kept) %>%
  ggplot(aes(x = Seroprevalence, y = Model)) +
  geom_boxplot(fill = "gold") +
  labs(
    x = "Seroprevalence (%)"
    # title = "ELISA-like: Simulation results",
    # subtitle = "F1 score"
  )

ggsave("./Simulations/ELISA_like/Output/ModelPerformance/SimModelSeroprevalence.png",
       plot = p_sero,
       height = 4.5,
       width = 7
)

p_sero_diff <- new_output_df %>%
  filter(Kept) %>%
  ggplot(aes(x = Seroprevalence_diff, y = Model)) +
  geom_boxplot(fill = "gold") +
  labs(
    x = "Difference in estimated seroprevalence to truth (%)"
    # title = "ELISA-like: Simulation results",
    # subtitle = "F1 score"
  )

ggsave("./Simulations/ELISA_like/Output/ModelPerformance/SimModelSeroprevalenceDiff.png",
       plot = p_sero_diff,
       height = 4.5,
       width = 7
)

# new_output_df %>%
#   filter(Kept) %>%
#   ggplot(aes(x = F1, y = Model)) +
#   geom_boxplot() +
#   geom_jitter(color="grey", size=0.7, alpha=0.5) +
#   facet_wrap(~Scenario) +
#   labs(
#     title = "Simulation results: F1"
#   )
#
# ggsave("./Current_results_f1_jittered.png")

p_time <- new_output_df %>%
  filter(Kept) %>%
  ggplot(aes(x = Time, y = Model)) +
  geom_boxplot(fill = "gold") +
  labs(
    title = "ELISA-like: Simulation results",
    subtitle = paste0("Time (", R, " MCMC iterations)")
  )

ggsave("./Simulations/ELISA_like/Output/ModelPerformance/SimModelTime.png",
       plot = p_time,
       height = 4.5,
       width = 7
)

# Save the output file.
fwrite(new_output_df, "./Simulations/ELISA_like/Output/SimModelPerformanceConvergenceIncluded.csv")


p_f1 / p_distance

# === Data example sim 3 ========================================================

elisa_sim <- read.csv("./Simulations/ELISA_like/Data/elisaLikeSeed3.csv")

elisa_sim$Labels <- factor(elisa_sim$Labels)
elisa_sim$Batch <- factor(elisa_sim$Batch)
elisa_sim$Fixed <- factor(elisa_sim$Fixed)
controls <- which(elisa_sim$Fixed == 1)

control_df <- elisa_sim[controls, ]

control_df$Controls <- "Controls"
elisa_sim$Controls <- "Full dataset"

p_data <- elisa_sim %>% 
  mutate(Plot_labels = ifelse(Fixed == 1, Labels, 3)) %>% 
  ggplot(aes(x = X_observed, 
             y = Y_observed,
             colour = factor(Plot_labels, labels = c("Seronegative", "Seropositive", "Unknown")),
             shape = factor(Fixed, labels = c("False", "True")))) +
  geom_point(size = 0.8) +
  # facet_wrap(~Controls) + 
  ggthemes::scale_color_colorblind() +
  labs(
    x = "Observed SPIKE OD",
    y = "Observed RBD OD",
    shape = "Control",
    colour = "Group"
  )



layout = "
AACC
BBCC
"

p_elisaLike <- wrap_plots(p_f1, p_distance, p_data, design =layout) + #widths = c(2, 1), ncol = 2) + # (p_observed / p_inferred) | p_sero +
  plot_annotation(tag_levels = "A")
# +



p_elisaLike

# === Inferred dataset =========================================================

mvt_samples <- readRDS("./Simulations/ELISA_like/Output/Chains/Sim3/MVT_chain_7_m_scale_1e-01_rho_11_theta_5.rds")

# Allocations
mvt_prob <- calcAllocProb(mvt_samples$alloc, eff_burn)
mvt_pred <- predictClass(mvt_prob)

# Inferred datasets
mvt_inferred_data <- rowMeans(mvt_samples$batch_corrected_data[, , -c(1:eff_burn)], dims = 2) %>%
  as_tibble(.name_repair = "minimal") %>%
  set_colnames(c("X_inferred", "Y_inferred")) %>%
  add_column(
    "Label" = factor(mvt_pred),
    "Prob" = apply(mvt_prob, 1, function(x) {
      x[which.max(x)]
    }),
    "Batch" = factor(elisa_sim$Batch),
    "Type" = "Inferred",
    "Fixed" = factor(elisa_sim$Fixed),
    "Model" = "MVT"
  )

p_inferred <- mvt_inferred_data %>%
  ggplot(aes(x = X_inferred, y = Y_inferred))  +
  geom_point(aes(shape = factor(Fixed, label = c("False", "True")),
                 colour = factor(Label,
                                 labels = c("Seronegative", "Seropositive")
                 ),
                 alpha = Prob),
             size = 0.8) +
  scale_alpha_continuous(range = c(0.4, 1.0)) +
  # geom_density_2d(aes(colour = Predicted)) +
  ggthemes::scale_color_colorblind() +  
  labs(
    # title = "Elisa data - batch adjusted",
    # subtitle = "Predicted labels",
    # caption = "log transformed"
    x = "Log-transformed batch-corrected SPIKE OD",
    y = "Log-transformed batch-corrected RBD OD",
    colour = "Group",
    shape = "Control",
    alpha = "Probability of\nallocation"
  )

layout = "
AACC
BBDD
"

p_elisaLike <- wrap_plots(p_f1, p_distance, p_data, p_inferred, design =layout) + #widths = c(2, 1), ncol = 2) + # (p_observed / p_inferred) | p_sero +
  plot_annotation(tag_levels = "A")

ggsave(paste0("./Simulations/ELISA_like/Output/elisaLikeResults.png"),
       plot = p_elisaLike,
       height = 8,
       width = 7
)
