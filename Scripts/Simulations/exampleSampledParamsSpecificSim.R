# parameters.R
# Using the output of ``modelling.R`` and ``findRemainingMisbehavingChains.R``
# create the plots of the sampled distributions from the simulations.
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

# Path to project directory
my_wd <- "./"

# Save the generated files to
save_path <- paste0(my_wd, "Simulations/Output/SampledParameters/")

# The MCMC number of iterations, etc.
R <- 15000
thin <- 25
burn <- 7500
eff_burn <- burn / thin

# The generated data and the description tibble
data_path <- paste0(my_wd, "Simulations/Generated_data/")
scenario_descriptions <- readRDS(paste0(data_path, "scenario_descriptions.rds"))

# Simulation outputs
param_df <- fread(paste0(my_wd, "Simulations/Output/SampledParameters.csv"))

# The fruit of the various convergence steps, the data frame indicating which
# chains have converged.
keep_chain_info <- fread(paste0(my_wd, "Simulations/Output/Convergence/SimsChainBehaviourDF.csv"))

# Performance scores
output_df <- fread(paste0(my_wd, "Simulations/Output/SimModelPerformance.csv"))

# The number of chains/sims/scenarios
n_scn <- nrow(scenario_descriptions)
n_sim <- 10
n_chains <- 10
scenarios <- scenario_descriptions$Scenario

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
  "Varying class representation",
  "Multivariate t generated"
) %>% set_names(scenario_descriptions$Scenario)

# === Parameter plots =========================================================

# Reduce down to the relevant parts of the data.frames
param_df_rel <- param_df[, c(1:7, 16:35, 79:84)] %>%
  left_join(keep_chain_info) %>%
  filter(Kept)

models_used <- param_df_rel$Model %>%
  unique()


ii <- 1
# for (ii in 1:n_scn) {
curr_scn <- scenario_descriptions$Scenario[ii]


# dir.create(save_dir)

true_mu_k <- scenario_descriptions$mu_k[[ii]]
true_m_b <- scenario_descriptions$m_b[[ii]]
true_S_b <- scenario_descriptions$S_b[[ii]]
true_t_df <- scenario_descriptions$dof[[ii]]

if (any(is.na(true_t_df))) {
  true_t_df <- NULL
}

.df <- param_df_rel %>%
  filter(Scenario == curr_scn, Simulation == 1) # %>%
# left_join(geweke_rel, by = c("Chain", "Simulation", "Scenario", "Model"))

mu_df <- .df %>%
  select(
    Iteration,
    Kept,
    Chain,
    Simulation,
    Scenario,
    Model,
    contains("Mu_")
  ) %>%
  filter(Iteration > burn) %>%
  pivot_longer(contains("Mu_"))

group_details <- mu_df$name %>%
  sapply(function(x) {
    c(substr(x, 4, 4), substr(x, 5, 5))
  }) %>%
  as.matrix() %>%
  t() %>%
  set_colnames(c("Class", "Dim"))

mu_df_full <- mu_df %>% cbind(group_details)

p_mu <- mu_df_full %>%
  ggplot(aes(
    x = value,
    fill = factor(Class), # labels = c("Group 1", "Group 2")),
    group = interaction(name, Chain, Simulation, Scenario, Model)
  )) +
  geom_density(alpha = 0.4) +
  facet_grid(Model ~ Dim, labeller = label_both) +
  labs(
    # title = scenario_strings[ii],
    # subtitle = "Sampled class mean",
    x = "Sampled class mean",
    y = "Density",
    fill = "Class"
  ) +
  geom_vline(xintercept = true_mu_k, colour = "red", lty = 2) +
  ggthemes::scale_fill_colorblind()


m_df <- .df %>%
  select(
    Iteration,
    Kept,
    Chain,
    Simulation,
    Scenario,
    Model,
    contains("m_")
  ) %>%
  filter(Iteration > burn) %>%
  pivot_longer(contains("m_"))

batch_details <- m_df$name %>%
  sapply(function(x) {
    c(substr(x, 3, 3), substr(x, 4, 4))
  }) %>%
  as.matrix() %>%
  t() %>%
  set_colnames(c("Batch", "Dim"))

m_df_full <- m_df %>% cbind(batch_details)


p_m <- m_df_full %>%
  ggplot(aes(
    x = value,
    fill = factor(Batch),
    group = interaction(name, Chain, Simulation, Scenario, Model)
  )) +
  geom_density(alpha = 0.4) +
  facet_grid(Model ~ Dim, labeller = label_both) +
  labs(
    # title = scenario_strings[ii],
    # subtitle = "Sampled class mean",
    x = "Sampled batch mean-effect",
    y = "Density",
    fill = "Batch"
  ) +
  geom_vline(xintercept = true_m_b, colour = "red", lty = 2) +
  ggthemes::scale_fill_colorblind()

S_df <- .df %>%
  select(
    Iteration,
    Kept,
    Chain,
    Simulation,
    Scenario,
    Model,
    contains("S_")
  ) %>%
  filter(Iteration > burn) %>%
  pivot_longer(contains("S_"))

S_df_full <- S_df %>% cbind(batch_details)

p_S <- S_df_full %>%
  ggplot(aes(
    x = value,
    fill = factor(Batch),
    group = interaction(name, Chain, Simulation, Scenario, Model)
  )) +
  geom_density(alpha = 0.4) +
  facet_grid(Model ~ Dim, labeller = label_both) +
  labs(
    # title = scenario_strings[ii],
    # subtitle = "Sampled class mean",
    x = "Sampled batch scaling",
    y = "Density",
    fill = "Batch"
  ) +
  geom_vline(xintercept = true_S_b, colour = "red", lty = 2) +
  ggthemes::scale_fill_colorblind()

t_df_df <- .df %>%
  select(
    Iteration,
    Kept,
    Chain,
    Simulation,
    Scenario,
    Model,
    contains("t_df_")
  ) %>%
  filter(Iteration > burn, Model == "MVT") %>%
  pivot_longer(contains("t_df_"))

df_details <- t_df_df$name %>%
  sapply(function(x) {
    c(substr(x, 6, 6))
  }) %>%
  matrix(ncol = 1) %>%
  set_colnames(c("Class"))

t_df_df_full <- t_df_df %>% cbind(df_details)

p_df <- t_df_df_full %>%
  ggplot(aes(
    x = value,
    fill = factor(Class),
    group = interaction(name, Chain, Simulation, Scenario, Model)
  )) +
  geom_density(alpha = 0.4) +
  facet_grid(~Class, labeller = label_both) +
  labs(
    # title = scenario_strings[ii],
    # subtitle = "Sampled Class degrees of freedom",
    x = "Sampled class degrees of freedom",
    y = "Density",
    fill = "Class"
  ) +
  geom_vline(xintercept = true_t_df, colour = "red", lty = 2) +
  ggthemes::scale_fill_colorblind()

p_patch <- p_mu + p_m + p_S + p_df +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


save_dir <- paste0(my_wd, "Simulations/Output/SampledParameters/")
ggsave(paste0(save_dir, "base_case_sim_1_sampled_params.png"),
  plot = p_patch,
  height = 6,
  width = 8
)
