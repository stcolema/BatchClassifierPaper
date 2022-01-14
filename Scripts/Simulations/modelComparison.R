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
library(mdiHelpR)
library(BatchMixtureModel)

# Set the gglot2 theme
setMyTheme()

# Set the random seed
set.seed(1)

# Save the generated files to
save_path <- "./Simulations/Output/"

# The MCMC number of iterations, etc.
R <- 15000
thin <- 25
burn <- 7500
eff_burn <- burn / thin

# The generated data and the description tibble
data_path <- "./Simulations/Generated_data/"
scenario_descriptions <- readRDS(paste0(data_path, "scenario_descriptions.rds"))

# The number of chains/sims/scenarios
n_scn <- nrow(scenario_descriptions)
n_sim <- 10
n_chains <- 10

# Nicer forms of the scenario names for title, labels, etc.
scenario_labels <- c(
  "Base case",
  "No batch effects",
  "Varying batch\nsize",
  "Varying batch\neffects",
  "Varying class\nrepresentation",
  "Multivariate t\ngenerated"
) %>% set_names(scenario_descriptions$Scenario)

scenario_strings <- c(
  "Base case",
  "No batch effects",
  "Varying batch size",
  "Varying batch effects",
  "Varying group representation",
  "Multivariate t generated"
) %>% set_names(scenario_descriptions$Scenario)

# === Model comparison =========================================================

# Read in the data regarding performance of each model
output_df <- fread("./Simulations/Output/SimModelPerformance.csv")

# Niver labels for the models
output_df$Model <- factor(output_df$Model,
  levels = unique(output_df$Model),
  labels = c("MVN", "MVT", "RF", "SVM", "LR", "LR - BC")
)

# The fruit of the various convergence steps, the data frame indicating which
# chains have converged.
keep_chain_info <- fread("./Simulations/Output/Convergence/SimsChainBehaviourDF.csv",
                         select = c("Chain", "Simulation", "Scenario", "Model", "Kept"))

dir.create("./Simulations/Output/ModelPerformance/")

# Separate out the mixture and ML results
mixture_output <- output_df[!is.na(output_df$BICmcmc), ]
ml_output <- output_df[is.na(output_df$BICmcmc), ]

# Add information about which chains have behaved
keep_chain_info$Seed <- keep_chain_info$Chain
mixture_output_new <- mixture_output %>%
  left_join(keep_chain_info, by = c("Model", "Scenario", "Simulation", "Seed")) %>%
  select(-Chain)

# We keep all the ML results
ml_output$Kept <- T

# The used output data frame
new_output_df <- rbind(mixture_output_new, ml_output)

# Plot of the difference between the estimated and true seroprevalence
p_sero_diff <- new_output_df %>%
  filter(Kept) %>%
  ggplot(aes(x = Seroprevalence_diff, y = Model)) +
  geom_boxplot(fill = "gold") +
  facet_wrap(~Scenario,
             labeller = labeller(Scenario = scenario_labels),
             scales = "free_x"
  ) +
  labs(
    x = "Difference between estimated and true seroprevalence (%)"
    # title = "Simulation results",
    # subtitle = "Squared Euclidean distance"
  )

ggsave("./Simulations/Output/ModelPerformance/SimModelSeroDiff.png",
  plot = p_sero_diff,
  height = 4.5,
  width = 7
)


# Plot of the distance between the probability allocation matrix and the truth
p_distance <- new_output_df %>%
  filter(Kept) %>%
  ggplot(aes(x = Distance, y = Model)) +
  geom_boxplot(fill = "gold") +
  facet_wrap(~Scenario,
    labeller = labeller(Scenario = scenario_labels),
    scales = "free_x"
  ) +
  labs(
    title = "Simulation results",
    subtitle = "Squared Euclidean distance"
  )

ggsave("./Simulations/Output/ModelPerformance/SimModelDistance.png",
  plot = p_distance,
  height = 4.5,
  width = 7
)

# Plot of the F1 score
p_f1 <- new_output_df %>%
  filter(Kept) %>%
  ggplot(aes(x = F1, y = Model)) +
  geom_boxplot(fill = "gold") +
  facet_wrap(~Scenario,
    labeller = labeller(Scenario = scenario_labels),
    scales = "free_x"
  ) +
  labs(
    title = "Simulation results",
    subtitle = "F1 score"
  )

ggsave("./Simulations/Output/ModelPerformance/SimModelF1.png",
  plot = p_f1,
  height = 4.5,
  width = 7
)

# Plot of the F1 score
p_f1_common_axis <- new_output_df %>%
  filter(Kept) %>%
  ggplot(aes(x = F1, y = Model)) +
  geom_boxplot(fill = "gold") +
  facet_wrap(~Scenario, labeller = labeller(Scenario = scenario_labels)) +
  labs(
    title = "Simulation results",
    subtitle = "F1 score"
  )

ggsave("./Simulations/Output/ModelPerformance/SimModelF1CommonAxis.png",
  plot = p_f1_common_axis,
  height = 4.5,
  width = 7
)

# Plot of the F1 score
p_f1_violin <- new_output_df %>%
  filter(Kept) %>%
  ggplot(aes(x = F1, y = Model)) +
  geom_violin(fill = "gold") +
  facet_wrap(~Scenario,
    labeller = labeller(Scenario = scenario_labels),
    scales = "free_x"
  ) +
  labs(
    title = "Simulation results",
    subtitle = "F1 score"
  )

ggsave("./Simulations/Output/ModelPerformance/SimModelF1Violin.png",
  plot = p_f1_violin,
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
  facet_wrap(~Scenario, labeller = labeller(Scenario = scenario_labels)) +
  labs(
    title = "Simulation results",
    subtitle = paste0("Time (", R, " MCMC iterations)")
  )

ggsave("./Simulations/Output/ModelPerformance/SimModelTime.png",
  plot = p_time,
  height = 4.5,
  width = 7
)

new_output_df %>%
  filter(Kept) %>%
  group_by(Model) %>% 
  summarise(Mean = mean(Time ), SD = sd(Time))
  
# Save the output file.
fwrite(new_output_df, "./Simulations/Output/SimModelPerformanceConvergenceIncluded.csv")

# === Paper plots ==============================================================

library(patchwork)
library(scales)

# Plot of the distance between the probability allocation matrix and the truth
p_distance <- new_output_df %>%
  filter(Kept) %>%
  ggplot(aes(x = Distance, y = Model)) +
  geom_boxplot(fill = "gold") +
  facet_wrap(~Scenario,
             labeller = labeller(Scenario = scenario_labels),
             scales = "free_x"
  ) + 
  labs(x = "Squared Euclidean distance") + 
  scale_x_continuous(breaks = pretty_breaks(3))


# Plot of the F1 score
p_f1 <- new_output_df %>%
  filter(Kept) %>%
  ggplot(aes(x = F1, y = Model)) +
  geom_boxplot(fill = "gold") +
  facet_wrap(~Scenario,
             labeller = labeller(Scenario = scenario_labels),
             scales = "free_x"
  ) + 
  labs(x = "F1 score") + 
  scale_x_continuous(breaks = pretty_breaks(3))

p_sero <- new_output_df %>%
  filter(Kept) %>%
  ggplot(aes(x = Seroprevalence_diff, y = Model)) +
  geom_boxplot(fill = "gold") +
  facet_wrap(~Scenario,
             labeller = labeller(Scenario = scenario_labels),
             scales = "free_x"
  ) + 
  labs(x = "Difference to true seroprevalence (%)") + 
  scale_x_continuous(breaks = pretty_breaks(3))

p_time <- new_output_df %>%
  filter(Kept) %>%
  ggplot(aes(x = Time, y = Model)) +
  geom_boxplot(fill = "gold") +
  facet_wrap(~Scenario,
             labeller = labeller(Scenario = scenario_labels),
             scales = "free_x"
  ) + 
  labs(x = "Time (s)") + 
  scale_x_continuous(breaks = pretty_breaks(3))


p_f1_dist_1 <- new_output_df  %>%
  filter(Kept) %>%
  pivot_longer(c(F1, Distance), names_to = "Score", values_to = "Value") %>% 
  ggplot(aes(x = Value, y = Model)) +
  geom_boxplot(fill = "gold") +
  facet_grid(Scenario ~ Score, 
    labeller = labeller(Scenario = scenario_labels), 
    scales = "free_x"
  )

ggsave("./Simulations/Output/ModelPerformance/SimModelF1Dist.png",
       plot = p_f1_dist_1,
       height = 7,
       width = 5
)

patchwork <- p_f1 / p_distance
patchwork + plot_annotation(tag_levels = 'A')

ggsave("./Simulations/Output/ModelPerformance/SimModelF1DistPatch.png",
       height = 7,
       width = 5
)

patchwork <- (p_f1 / p_distance / p_sero)
p_comp <- patchwork + plot_annotation(tag_levels = 'A')

ggsave("./Simulations/Output/ModelPerformance/SimModelComparison.png",
  plot = p_comp,
  height = 10,
  width = 6
)
