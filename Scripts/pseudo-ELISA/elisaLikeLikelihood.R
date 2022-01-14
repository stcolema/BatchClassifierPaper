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
library(patchwork)

# on_hpc <- TRUE
# if(on_hpc) {
#   setwd("/home/sdc56/rds/hpc-work/BatchMixtureModel/")
# }

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

# === Likelihood plots =========================================================

param_df <- fread("./Simulations/ELISA_like/Output/SampledParameters.csv")
keep_chain_info <- fread("./Simulations/ELISA_like/Output/Convergence/GewekeDF.csv", select = 4:7) %>%
  distinct()

param_df_rel <- param_df %>%
  select(
    Iteration,
    Observed_likelihood,
    Complete_likelihood,
    Simulation,
    Chain,
    Model
  ) %>%
  left_join(keep_chain_info)


models_used <- param_df_rel$Model %>%
  unique()

dir.create("./Simulations/ELISA_like/Output/LikelihoodPlots/")
# === Manual reduction =========================================================
# Some chains settle in a local minimum, remove these
mvt_chains_dropped <- data.frame(
  Simulation = c(rep(1,8), rep(2, 9), rep(3, 1), rep(4, 8), rep(5, 8), rep(7, 1), rep(9, 8)),
  Chain = c(c(1, 2, 4, 5, 6, 7, 8, 10),
            c(1, 2, 3, 5, 6, 7, 8, 9, 10),
            c(3), 
            c(2:8, 10), 
            c(1, c(4:10)), 
            c(7), 
            c(1:3, 5, 7:10)
            ),
  Kept = FALSE,
  Model = "MVT"
)

# param_df_rel %>%
#   filter(Model == "MVT", Iteration > burn, Kept, Simulation == 9, Chain %in% 1:6:10) %>%
#   ggplot(aes(x = Iteration, y = Complete_likelihood, group = Chain, colour = factor(Chain))) +
#   facet_wrap(~Simulation) +
#   geom_line() +
#   labs(
#     title = "ELISA-like",
#     y = "Complete likelihood"
#   ) +
#   theme(axis.text.x = element_text(angle = 30)) +
#   theme(panel.spacing.x = unit(1.2, "lines")) +
#   ggthemes::scale_color_colorblind()

mvn_chains_dropped <- data.frame(
  Simulation = c(rep(1,7), rep(2, 9), rep(5, 8), rep(9, 6)),
  Chain = c(1, 2, 4, 6:9, 1:3, 5:10, 1, 3:5, 7:10, 1, 2, 5, 7, 8, 10),
  Kept = FALSE,
  Model = "MVN"
)

chains_dropped <- rbind(mvn_chains_dropped, mvt_chains_dropped)

new_df <- left_join(param_df_rel, chains_dropped, by = c("Chain", "Simulation", "Model")) %>% 
  mutate(Kept.y = ifelse(is.na(Kept.y), T, F ),
         Kept = ifelse(Kept.x & Kept.y, T, F)
  ) %>% select(-c(Kept.x, Kept.y))

# === Plotting =================================================================
p_lst <- vector("list", 3)

for (.mod in models_used) {
  p_observed_full <- param_df_rel %>%
    filter(Model == .mod, Iteration > burn) %>%
    ggplot(aes(x = Iteration, y = Observed_likelihood, group = Chain)) +
    facet_wrap(~Simulation) +
    geom_line() +
    labs(
      title = "ELISA-like",
      subtitle = .mod,
      y = "Observed likelihood"
    ) +
    theme(axis.text.x = element_text(angle = 30)) +
    theme(panel.spacing.x = unit(1.2, "lines"))
  
  file_name <- paste0(.mod, "_observed_likelihood.png")
  ggsave(paste0("./Simulations/ELISA_like/Output/LikelihoodPlots/", file_name), 
    p_observed_full)
  
  p_complete_full <- param_df_rel %>%
    filter(Model == .mod, Iteration > burn) %>%
    ggplot(aes(x = Iteration, y = Complete_likelihood, group = Chain)) +
    facet_wrap(~Simulation) +
    geom_line() +
    labs(
      title = "ELISA-like",
      subtitle = .mod,
      y = "Complete likelihood"
    ) +
    theme(axis.text.x = element_text(angle = 30)) +
    theme(panel.spacing.x = unit(1.2, "lines"))
  
  if(.mod == "MVT") {
    p_lst[[1]] <- p_complete_full
  }
  
  file_name <- paste0(.mod, "_complete_likelihood.png")
  ggsave(paste0("./Simulations/ELISA_like/Output/LikelihoodPlots/", file_name), p_complete_full)
  
  p_observed_conv <- param_df_rel %>%
    filter(Model == .mod, Iteration > burn, Kept) %>%
    ggplot(aes(x = Iteration, y = Observed_likelihood, group = Chain)) +
    facet_wrap(~Simulation) +
    geom_line() +
    labs(
      title = "ELISA-like",
      subtitle = .mod,
      y = "Observed likelihood"
    ) +
    theme(axis.text.x = element_text(angle = 30)) +
    theme(panel.spacing.x = unit(1.2, "lines"))
  
  file_name <- paste0(.mod, "_observed_likelihood_converged.png")
  ggsave(paste0("./Simulations/ELISA_like/Output/LikelihoodPlots/", file_name),
    p_observed_conv)
  
  p_complete_conv <- param_df_rel %>%
    filter(Model == .mod, Iteration > burn, Kept) %>%
    ggplot(aes(x = Iteration, y = Complete_likelihood, group = Chain)) +
    facet_wrap(~Simulation) +
    geom_line() +
    labs(
      title = "ELISA-like",
      subtitle = .mod,
      y = "Complete likelihood"
    ) +
    theme(axis.text.x = element_text(angle = 30)) +
    theme(panel.spacing.x = unit(1.2, "lines"))
  
  if(.mod == "MVT") {
    p_lst[[2]] <- p_complete_conv
  }
  
  file_name <- paste0(.mod, "_complete_likelihood_converged.png")
  ggsave(paste0("./Simulations/ELISA_like/Output/LikelihoodPlots/", file_name), p_complete_conv)
  
  p_complete_true_conv <- new_df %>%
    filter(Model == .mod, Iteration > burn, Kept) %>%
    ggplot(aes(x = Iteration, y = Complete_likelihood, group = Chain)) +
    facet_wrap(~Simulation) +
    geom_line() +
    labs(
      title = "ELISA-like",
      subtitle = .mod,
      y = "Complete likelihood"
    ) +
    theme(axis.text.x = element_text(angle = 30)) +
    theme(panel.spacing.x = unit(1.2, "lines"))
  
  if(.mod == "MVT") {
    p_lst[[3]] <- p_complete_true_conv
  }
  
  file_name <- paste0(.mod, "_complete_likelihood_true_converged.png")
  ggsave(paste0("./Simulations/ELISA_like/Output/LikelihoodPlots/", file_name), p_complete_true_conv)
  
}

convergence_df <- new_df[, 4:7] %>% 
  distinct()

fwrite(convergence_df, file = "./Simulations/ELISA_like/Output/Convergence/GewekeAndLikelihood.csv")

p_1 <- param_df_rel %>%
  filter(Model == "MVT", Iteration > burn) %>%
  ggplot(aes(x = Iteration, y = Complete_likelihood, group = Chain)) +
  facet_wrap(~Simulation, ncol = 1) +
  geom_line() +
  labs(
    # title = "ELISA-like",
    # subtitle = .mod,
    y = "Complete log-likelihood"
  ) +
  theme(axis.text.x = element_text(angle = 30)) +
  theme(panel.spacing.x = unit(1.2, "lines")) +
  scale_y_continuous(breaks = scales::pretty_breaks(3))


p_2 <- param_df_rel %>%
  filter(Model == "MVT", Iteration > burn, Kept) %>%
  ggplot(aes(x = Iteration, y = Complete_likelihood, group = Chain)) +
  facet_wrap(~Simulation, ncol = 1) +
  geom_line() +
  labs(
    # title = "ELISA-like",
    # subtitle = .mod,
    y = "Complete log-likelihood"
  ) +
  theme(axis.text.x = element_text(angle = 30)) +
  theme(panel.spacing.x = unit(1.2, "lines")) +
  scale_y_continuous(breaks = scales::pretty_breaks(3))


p_3 <- new_df %>%
  filter(Model == "MVT", Iteration > burn, Kept) %>%
  ggplot(aes(x = Iteration, y = Complete_likelihood, group = Chain)) +
  facet_wrap(~Simulation, ncol = 1) +
  geom_line() +
  labs(
    y = "Complete log-likelihood"
  ) +
  theme(axis.text.x = element_text(angle = 30)) +
  theme(panel.spacing.x = unit(1.2, "lines")) +
  scale_y_continuous(breaks = scales::pretty_breaks(3))


p_patchwork <- (p_1 | p_2 | p_3) +
  plot_annotation(tag_levels = 'A') 

ggsave("./Simulations/ELISA_like/Output/LikelihoodPlots/MVT_model_likelihoods_comparison.png", 
       p_patchwork, height = 8, width = 6)
