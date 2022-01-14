# completeLikelihood.R
# Using the output of ``gewekeConvergence.R`` create the plots of the model fit,
# measured using the complete likelihood, and find chains which are not
# converged and were not found by the Geweke diagnostic test.
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

library(patchwork)

# Set the gglot2 theme
setMyTheme()

# Set the random seed
set.seed(1)

# Path from working directory to project
my_wd <- "./"

# Save the generated files to
save_path <- paste0(my_wd, "Simulations/Output/")

# The MCMC number of iterations, etc.
R <- 15000
thin <- 25
burn <- 7500
eff_burn <- burn / thin

# The generated data and the description tibble
data_path <- paste0(my_wd, "Simulations/Generated_data/")
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
  "Varying class representation",
  "Multivariate t generated"
) %>% set_names(scenario_descriptions$Scenario)

# === Likelihood plots =========================================================

used_columns <- c(
  "Iteration",
  "Observed_likelihood",
  "Complete_likelihood",
  "Scenario",
  "Simulation",
  "Chain",
  "Model"
)

param_df <- fread(
  paste0(my_wd, "Simulations/Output/SampledParameters.csv"),
  select = used_columns
)

keep_chain_info <- fread(paste0(my_wd, "Simulations/Output/Convergence/GewekeDF.csv"))

param_df_rel <- param_df %>%
  left_join(keep_chain_info)

models_used <- param_df_rel$Model %>%
  unique()

dir.create(paste0(my_wd, "Simulations/Output/LikelihoodPlots/"))
for (ii in 1:n_scn) {
  scn_df <- param_df_rel %>%
    filter((Scenario == scenario_descriptions$Scenario[ii]) & (Iteration > burn))

  for (.mod in models_used) {
    p_observed_full <- scn_df %>%
      filter(Model == .mod) %>%
      ggplot(aes(x = Iteration, y = Observed_likelihood, group = Chain)) +
      facet_wrap(~Simulation) +
      geom_line() +
      labs(
        title = scenario_strings[ii],
        subtitle = .mod,
        y = "Observed likelihood"
      ) +
      theme(axis.text.x = element_text(angle = 30)) +
      theme(panel.spacing.x = unit(1.2, "lines"))

    file_name <- paste0(scenario_descriptions$Scenario[ii], "_", .mod, "_observed_likelihood.png")
    ggsave(paste0(my_wd, "Simulations/Output/LikelihoodPlots/", file_name), p_observed_full)

    p_complete_full <- scn_df %>%
      filter(Model == .mod) %>%
      ggplot(aes(x = Iteration, y = Complete_likelihood, group = Chain)) +
      facet_wrap(~Simulation) +
      geom_line() +
      labs(
        title = scenario_strings[ii],
        subtitle = .mod,
        y = "Complete likelihood"
      ) +
      theme(axis.text.x = element_text(angle = 30)) +
      theme(panel.spacing.x = unit(1.2, "lines"))

    file_name <- paste0(scenario_descriptions$Scenario[ii], "_", .mod, "_complete_likelihood.png")
    ggsave(paste0(my_wd, "Simulations/Output/LikelihoodPlots/", file_name), p_complete_full)

    p_observed_conv <- scn_df %>%
      filter(Model == .mod & Kept) %>%
      ggplot(aes(x = Iteration, y = Observed_likelihood, group = Chain)) +
      facet_wrap(~Simulation) +
      geom_line() +
      labs(
        title = scenario_strings[ii],
        subtitle = .mod,
        y = "Observed likelihood"
      ) +
      theme(axis.text.x = element_text(angle = 30)) +
      theme(panel.spacing.x = unit(1.2, "lines"))

    file_name <- paste0(scenario_descriptions$Scenario[ii], "_", .mod, "_observed_likelihood_converged.png")
    ggsave(paste0(my_wd, "Simulations/Output/LikelihoodPlots/", file_name), p_observed_conv)

    p_complete_conv <- scn_df %>%
      filter(Model == .mod & Kept) %>%
      ggplot(aes(x = Iteration, y = Complete_likelihood, group = Chain)) +
      facet_wrap(~Simulation) +
      geom_line() +
      labs(
        title = scenario_strings[ii],
        subtitle = .mod,
        y = "Complete likelihood"
      ) +
      theme(axis.text.x = element_text(angle = 30)) +
      theme(panel.spacing.x = unit(1.2, "lines"))

    file_name <- paste0(scenario_descriptions$Scenario[ii], "_", .mod, "_complete_likelihood_converged.png")
    ggsave(paste0(my_wd, "Simulations/Output/LikelihoodPlots/", file_name), p_complete_conv)
  }
}

# === Misbehaving chains =======================================================

# The following (Model, Simulation, Scenario) combinations contain at least one
# mis-behaving chain based on the likelihood plots even after reducing the
# number of chains using the Geweke statistic and a test for Normality.
misbehaving_cases <- tribble(
  ~Scenario, ~Simulation, ~Model,
  "MVT", 7L, "MVN",
  # "MVT", 6L, "MVN",
  "MVT", 10L, "MVN"

  # "MVT", 7L, "MVN",
  # "MVT", 10L, "MVN",
  # "MVT", 10L, "MVN",
  # "MVT", 10L, "MVN",
  # "MVT",  10L,  "MVN",
  # "varying_batch_effect", 1L, "MVN"
)

n_misbehaved <- nrow(misbehaving_cases)

for (ii in seq(1, n_misbehaved)) {
  curr_scn <- misbehaving_cases$Scenario[ii]
  curr_sim <- misbehaving_cases$Simulation[ii]
  curr_model <- misbehaving_cases$Model[ii]

  rel_df <- param_df_rel %>%
    filter(Simulation == curr_sim & Model == curr_model & Scenario == curr_scn & Iteration > burn & Kept)

  if (ii == 1) {
    find_chains_df <- rel_df
  } else {
    find_chains_df <- rbind(find_chains_df, rel_df)
  }
}

find_chains_df %>%
  # filter(Simulation == 10) %>%
  select(Scenario, Simulation, Chain, Kept) %>%
  distinct()

find_chains_df %>%
  filter(Simulation == 7, Chain %in% 1:9) %>%
  ggplot(aes(x = Iteration, y = Complete_likelihood, colour = factor(Chain))) +
  geom_line() +
  facet_wrap(~ Simulation + Model + Scenario, labeller = label_both, scales = "free") +
  scale_color_viridis_d()

find_chains_df %>%
  filter(Simulation == 10, Chain %in% c(2, 3, 6, 9)) %>%
  ggplot(aes(x = Iteration, y = Complete_likelihood, colour = factor(Chain))) +
  geom_line() +
  facet_wrap(~ Simulation + Model + Scenario, labeller = label_both, scales = "free") +
  scale_color_viridis_d()

dropped <- c(10, 1, 4, 8, 10)

misbehaving_cases <- misbehaving_cases[c(rep(1, 1), rep(2, 4)), ]
misbehaving_cases$Chain <- dropped

# These are the misbehaving chains in each case
# misbehaving_cases$Chain <- c(1, 1, 1, 6, 1)

# misbehaving_cases <- misbehaving_cases[c(3, 1, 2, 2, 2),]
# misbehaving_cases$Chain <- c(1, 1, 5, 8, 8)

# Quick check that dropping these chains is correct
for (ii in seq(1, n_misbehaved)) {
  curr_scn <- misbehaving_cases$Scenario[ii]
  curr_sim <- misbehaving_cases$Simulation[ii]
  curr_model <- misbehaving_cases$Model[ii]
  curr_chain <- misbehaving_cases %>%
    filter(
      Scenario == curr_scn,
      Simulation == curr_sim,
      Model == curr_model
    ) %>%
    select(Chain) %>%
    unlist()

  # The current simulation in the current scenario for the current model
  rel_df <- param_df_rel %>%
    filter(Simulation == curr_sim &
      Model == curr_model &
      Scenario == curr_scn &
      Iteration > burn &
      Kept)

  if (ii == 1) {

    # The samples for the well-behaved chains
    good_chains_df <- rel_df %>%
      filter(Chain != curr_chain)

    # All the chains
    all_chains_df <- rel_df
  } else {
    good_chains_df <- rbind(good_chains_df, filter(rel_df, !Chain %in% curr_chain))
    all_chains_df <- rbind(all_chains_df, rel_df)
  }

  # Update the information about which chains are good to keep
  keep_chain_info$Kept[
    keep_chain_info$Model == curr_model &
      keep_chain_info$Simulation == curr_sim &
      keep_chain_info$Scenario == curr_scn &
      keep_chain_info$Chain %in% curr_chain
  ] <- FALSE
}

# good_chains_df %>%
#   filter(Simulation == 10) %>%
#   select(-c(Iteration, Observed_likelihood, Complete_likelihood)) %>%
# distinct()

p_good_chains <- good_chains_df %>%
  ggplot(aes(x = Iteration, y = Complete_likelihood, colour = factor(Chain))) +
  geom_line() +
  facet_wrap(~ Simulation + Model + Scenario, labeller = label_both) +
  scale_color_viridis_d() +
  labs(
    title = "Simulations with misbehaving chains",
    subtitle = "Well-behaved chains"
  )

p_all_chains <- all_chains_df %>%
  ggplot(aes(x = Iteration, y = Complete_likelihood, colour = factor(Chain))) +
  geom_line() +
  facet_wrap(~ Simulation + Model + Scenario, labeller = label_both) +
  scale_color_viridis_d() +
  labs(
    title = "Simulations with misbehaving chains",
    subtitle = "Manually identified chains included"
  )

# Save the trace plots of the complete log-likelihood for the simulations
# where poorly behaved chains were manually identified.
ggsave(paste0(my_wd, "Simulations/Output/Convergence/well_behaved_chains.png"),
  plot = p_good_chains,
  width = 7,
  height = 5
)

ggsave(paste0(my_wd, "Simulations/Output/Convergence/manually_identified_chains.png"),
  plot = p_all_chains,
  width = 7,
  height = 5
)


# Save the data.frame indicating which chains are well-behaved
fwrite(
  keep_chain_info,
  paste0(my_wd, "Simulations/Output/Convergence/SimsChainBehaviourDF.csv")
)

# === Supplementary material plot ==============================================

p_1 <- param_df_rel %>%
  filter(Scenario == "MVT" &
    Model == "MVN" &
    Iteration > burn &
    Simulation == 7) %>%
  ggplot(aes(x = Iteration, y = Complete_likelihood, group = Chain)) +
  geom_line() +
  labs(y = "Complete log-likelihood") +
  scale_x_continuous(breaks = c(7500, 11250, 15000))
# +
#   facet_wrap(~Simulation, ncol = 1)

p_2 <- param_df_rel %>%
  filter(Scenario == "MVT" &
    Model == "MVN" &
    Iteration > burn &
    Simulation == 7 &
    Kept) %>%
  ggplot(aes(x = Iteration, y = Complete_likelihood, group = Chain)) +
  geom_line() +
  labs(y = "Complete log-likelihood") +
  scale_x_continuous(breaks = c(7500, 11250, 15000))
# +
# facet_wrap(~Simulation, ncol = 1)

p_3 <- param_df_rel %>%
  select(-Kept) %>%
  left_join(keep_chain_info) %>%
  filter(Scenario == "MVT" &
    Model == "MVN" &
    Iteration > burn &
    Simulation == 7 &
    Kept) %>%
  ggplot(aes(x = Iteration, y = Complete_likelihood, group = Chain)) +
  geom_line() +
  labs(y = "Complete log-likelihood") +
  scale_x_continuous(breaks = c(7500, 11250, 15000)) +
  ylim(-2714, -2543)
# +
#   facet_wrap(~Simulation, ncol = 1)

patchwork <- (p_1 + p_2 + p_3)
p_comp <- patchwork + plot_annotation(tag_levels = "A")

ggsave(
  paste0(
    my_wd,
    "Simulations/Output/Convergence/chain_reduction_process_MVT_MVN_7.png"
  ),
plot = p_comp,
width = 7.5,
height = 5
)
