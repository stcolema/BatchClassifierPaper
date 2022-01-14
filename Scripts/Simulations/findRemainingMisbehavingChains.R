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
  "tibble"
)

install.packages(setdiff(cran_packages, rownames(installed.packages())))

library(tidyr)
library(ggplot2)
library(dplyr)
library(magrittr)
library(tibble)
library(coda)
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

columns_used <- c(
  "Iteration",
  "Observed_likelihood",
  "Complete_likelihood",
  "Scenario",
  "Simulation",
  "Chain",
  "Model"
)

# Simulation outputs
param_df <- data.table::fread("./Simulations/Output/SampledParameters.csv",
  select = columns_used
)
geweke_df <- read.csv("./Simulations/Output/Convergence/GewekeDF.csv")
output_df <- read.csv("./Simulations/Output/SimModelPerformance.csv")

# Keep only the relevant parts of the geweke matrix
geweke_rel <- geweke_df # geweke_df[, c(4:8)] %>%
#     distinct()

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
  "Varying class representation",
  "Multivariate t generated"
) %>% set_names(scenario_descriptions$Scenario)

# === Likelihood plots =========================================================

# Reduce down to the relevant parts of the data.frames
param_df_rel <- param_df  %>%
  # select(
  #   Iteration,
  #   Observed_likelihood,
  #   Complete_likelihood,
  #   Scenario,
  #   Simulation,
  #   Chain,
  #   Model
  # ) %>%
  left_join(geweke_rel)

models_used <- param_df_rel$Model %>%
  unique()


# The following (Model, Simulation, Scenario) combinations contain at least one
# mis-behaving chain based on the likelihood plots even after reducing the
# number of chains using the Geweke statistic and a test for Normality.
misbehaving_cases <- tribble(
  ~Scenario, ~Sim, ~Model,
  "MVT", 7L, "MVN"
  
  # "MVT", 1L, "MVN",
  # "MVT", 1L, "MVT",
  # "MVT",  10L,  "MVN",
  # "varying_batch_effect", 1L, "MVN",
  # "varying_batch_representation", 3L, "MVT",
  # "varying_batch_size", 1L, "MVN"
)

n_misbehaved <- nrow(misbehaving_cases)

for (ii in seq(1, n_misbehaved)) {
  curr_scn <- misbehaving_cases$Scenario[ii]
  curr_sim <- misbehaving_cases$Sim[ii]
  curr_model <- misbehaving_cases$Model[ii]

  rel_df <- param_df_rel %>%
    filter(Simulation == curr_sim 
           & Model == curr_model 
           & Scenario == curr_scn 
           & Iteration > burn 
           & Kept)

  if (ii == 1) {
    find_chains_df <- rel_df
  } else {
    find_chains_df <- rbind(find_chains_df, rel_df)
  }
}

find_chains_df %>%
  ggplot(aes(x = Iteration, y = Complete_likelihood, colour = factor(Chain))) +
  geom_line() +
  facet_wrap(~ Simulation + Model + Scenario, labeller = label_both) +
  scale_color_viridis_d()

# These are the misbehaving chains in each case
misbehaving_cases$Chain <- c(10)

# Quick check that dropping these chains is correct
for (ii in seq(1, n_misbehaved)) {
  curr_scn <- misbehaving_cases$Scenario[ii]
  curr_sim <- misbehaving_cases$Sim[ii]
  curr_model <- misbehaving_cases$Model[ii]
  curr_chain <- misbehaving_cases$Chain[ii]
  rel_df <- param_df_rel %>%
    filter(Simulation == curr_sim & Model == curr_model & Scenario == curr_scn & Iteration > burn & Kept)

  if (ii == 1) {

    # The samples for the well-behaved chaines
    good_chains_df <- rel_df %>%
        filter(Chain != curr_chain)

    # All the chains
    all_chains_df <- rel_df 
  } else {
    good_chains_df <- rbind(good_chains_df, filter(rel_df, Chain != curr_chain))
    all_chains_df <- rbind(all_chains_df, rel_df)
  }
}

good_chains_df %>%
  ggplot(aes(x = Iteration, y = Complete_likelihood, colour = factor(Chain))) +
  geom_line() +
  facet_wrap(~ Simulation + Model + Scenario, labeller = label_both) +
  scale_color_viridis_d()

all_chains_df %>%
  ggplot(aes(x = Iteration, y = Complete_likelihood, colour = factor(Chain))) +
  geom_line() +
  facet_wrap(~ Simulation + Model + Scenario, labeller = label_both) +
  scale_color_viridis_d()

# Now change these chains to not kept in the output df
for (ii in seq(1, n_misbehaved)) {
  curr_scn <- misbehaving_cases$Scenario[ii]
  curr_sim <- misbehaving_cases$Sim[ii]
  curr_model <- misbehaving_cases$Model[ii]
  curr_chain <- misbehaving_cases$Chain[ii]

  output_df$Kept[output_df$Model == curr_model &
    output_df$Simulation == curr_sim &
    output_df$Scenario == curr_scn &
    output_df$Seed == curr_chain] <- F
}

# === Model comparison =========================================================

dtt$model <- factor(dtt$model, levels = c("mb", "ma", "mc"), labels = c("MBB", "MAA", "MCC"))

output_df$Model %>% unique()
output_df$Model <- factor(output_df$Model,
  levels = unique(output_df$Model),
  labels = c("MVN", "MVT", "RF", "SVM", "LR", "LR (batch corrected)")
)

output_df %>%
  filter(Kept) %>%
  ggplot(aes(x = Distance, y = Model)) +
  geom_boxplot(fill = "gold") +
  facet_wrap(~Scenario, labeller = labeller(Scenario = scenario_labels), scales = "free_x") +
  labs(
    title = "Simulation results",
    subtitle = "Squared Euclidean distance"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))

ggsave("./Simulations/Output/ModelPerformance/SimModelDistanceReduced.png",
  height = 4.5, width = 7
)


# additional_f1_drops <- which(new_output_df$F1 < 0.2)

output_df %>%
  filter(Kept) %>%
  ggplot(aes(x = F1, y = Model)) +
  geom_boxplot(fill = "gold") +
  facet_wrap(~Scenario, labeller = labeller(Scenario = scenario_labels), scales = "free_x") +
  labs(
    title = "Simulation results",
    subtitle = "F1 score"
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))

ggsave("./Simulations/Output/ModelPerformance/SimModelF1Reduced.png",
  height = 4.5, width = 7
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
#
# output_df %>%
#   filter(Kept) %>%
#   ggplot(aes(x = Time, y = Model)) +
#   geom_boxplot(fill = "gold") +
#   facet_wrap(~Scenario, labeller = labeller(Scenario = scenario_labels)) +
#   labs(
#     title = "Simulation results",
#     subtitle = paste0("Time (", R, " MCMC iterations)")
#   )
#
# ggsave("./Simulations/Output/ModelPerformance/SimModelTimeReduced.png")

# write.csv(output_df, "./Simulations/Output/SimModelPerformance.csv")

# === Parameter plots ==========================================================

# geweke_df <- read.csv("./Simulations/Output/Convergence/GewekeDF.csv", row.names = 1)

geweke_rel <- geweke_df[, c(4:8)]
param_df_rel <- param_df[, c(1:7, 16:35, 79:84)]

for (ii in 1:n_scn) {
  curr_scn <- scenario_descriptions$Scenario[ii]

  save_dir <- paste0("./Simulations/Output/", curr_scn)
  dir.create(save_dir)

  true_mu_k <- scenario_descriptions$mu_k[[ii]]
  true_m_b <- scenario_descriptions$m_b[[ii]]
  true_S_b <- scenario_descriptions$S_b[[ii]]
  true_t_df <- scenario_descriptions$dof[[ii]]

  .df <- param_df_rel %>%
    filter(Scenario == curr_scn) %>%
    left_join(geweke_rel, by = c("Chain", "Simulation", "Scenario", "Model"))

  p_mu <- .df %>%
    select(
      Iteration,
      Kept,
      Chain,
      Simulation,
      Scenario,
      Model,
      contains("Mu_")
    ) %>%
    filter(Iteration > burn, Kept) %>%
    pivot_longer(contains("Mu_")) %>%
    ggplot(aes(
      x = Iteration,
      y = value,
      group = interaction(Chain, Simulation, Scenario, Model)
    )) +
    geom_line() +
    facet_wrap(~name) +
    labs(
      title = scenario_strings[ii],
      subtitle = "Sampled class mean"
    ) +
    geom_hline(yintercept = true_mu_k, colour = "red", lty = 2)

  ggsave(paste0(save_dir, "/sampled_class_means.png"), plot = p_mu)

  p_m <- .df %>%
    select(
      Iteration,
      Kept,
      Chain,
      Simulation,
      Scenario,
      Model,
      contains("m_")
    ) %>%
    filter(Iteration > burn, Kept) %>%
    pivot_longer(contains("m_")) %>%
    ggplot(aes(
      x = Iteration,
      y = value,
      group = interaction(Chain, Simulation, Scenario, Model)
    )) +
    geom_line() +
    facet_wrap(~name, labeller = label_bquote(cols = .(name))) +
    geom_hline(yintercept = 0, colour = "red") +
    labs(
      title = scenario_strings[ii],
      subtitle = "Sampled batch mean effects"
    ) +
    geom_hline(yintercept = true_m_b, colour = "red", lty = 2)

  ggsave(paste0(save_dir, "/sampled_batch_means.png"), plot = p_m)

  p_S <- .df %>%
    select(
      Iteration,
      Kept,
      Chain,
      Simulation,
      Scenario,
      Model,
      contains("S_")
    ) %>%
    filter(Iteration > burn, Kept) %>%
    pivot_longer(contains("S_")) %>%
    ggplot(aes(
      x = Iteration,
      y = value,
      group = interaction(Chain, Simulation, Scenario, Model)
    )) +
    geom_point() +
    facet_wrap(~name, labeller = label_bquote(cols = .(name))) +
    # geom_hline(yintercept = 1.0, colour = "red") +
    labs(
      title = scenario_strings[ii],
      subtitle = "Sampled batch scale effects"
    ) +
    geom_hline(yintercept = true_S_b, colour = "red", lty = 2)

  ggsave(paste0(save_dir, "/sampled_batch_scale.png"), plot = p_S)

  p_df <- .df %>%
    select(
      Iteration,
      Kept,
      Chain,
      Simulation,
      Scenario,
      Model,
      contains("t_df_")
    ) %>%
    filter(Iteration > burn, Kept, Model == "MVT") %>%
    pivot_longer(contains("t_df_")) %>%
    ggplot(aes(
      x = Iteration,
      y = value,
      group = interaction(Chain, Simulation, Scenario, Model)
    )) +
    geom_point() +
    facet_wrap(~name, labeller = label_bquote(cols = .(name)), scales = "free_y", ncol = 1) +
    # geom_hline(yintercept = 0, colour = "red") +
    labs(
      title = scenario_strings[ii],
      subtitle = "Sampled class degrees of freedom"
    )

  if (!is.na(true_t_df)) {
    p_df <- p_df +
      geom_hline(yintercept = true_t_df, colour = "red", lty = 2)
  }

  ggsave(paste0(save_dir, "/sampled_class_dof.png"), plot = p_df)
}

# Time table
output_df %>%
  group_by(Model) %>%
  summarise(m_t = mean(Time), sd_t = sd(Time))
