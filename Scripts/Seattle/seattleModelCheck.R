
library(dplyr)
library(ggplot2)
library(tidyr)

mdiHelpR::setMyTheme()

dataset <- "Seattle"

# MCMC parameters
R <- 25000
thin <- 50
burn <- 15000
eff_burn <- burn / thin

# Number of chains to run
n_chains <- 10L

home_dir <- paste0("./Analysis/", dataset, "/ModelChecks/")
save_dir <- paste0("./Analysis/", dataset, "/Outputs/Plots/")

batch_acceptance <- read.csv(paste0(home_dir, "batchParamAcceptance.csv"))
group_acceptance <- read.csv(paste0(home_dir, "groupParamAcceptance.csv"))
likelihood_df <- read.csv(paste0(home_dir, "likelihoods.csv"))

p_batch <- batch_acceptance %>%
  pivot_longer(c(S, m), names_to = "Parameter", values_to = "Acceptance_rate") %>%
  ggplot(aes(x = interaction(Parameter, B), y = Acceptance_rate, colour = Parameter)) +
  geom_point() +
  geom_boxplot() +
  facet_grid(m_scale + rho + theta ~ Model)


p_group <- group_acceptance %>%
  pivot_longer(c(Cov, Mu, DF), names_to = "Parameter", values_to = "Acceptance_rate") %>%
  ggplot(aes(x = interaction(Parameter, K), y = Acceptance_rate, colour = Parameter)) +
  geom_point() +
  geom_boxplot() +
  facet_grid(m_scale + rho + theta ~ Model, scales = "free_x")


likelihood_df %>%
  filter(param == "complete") %>%
  ggplot(aes(x = Iteration, y = value, group = factor(Chain))) +
  geom_line() +
  facet_grid(m_scale + rho + theta ~ model)

p_likelihood <- likelihood_df %>%
  filter(param == "complete", Iteration > burn) %>%
  ggplot(aes(x = Iteration, y = value, colour = factor(Chain))) +
  geom_line() +
  facet_grid(m_scale + rho + theta ~ model, labeller = label_both)
# +
  # ggthemes::scale_color_colorblind()


p_good_chains <- likelihood_df %>%
  filter(param == "complete", Iteration > burn, Chain %in% c(3, 5, 7, 10)) %>%
  ggplot(aes(x = Iteration, y = value, colour = factor(Chain))) +
  geom_line() +
  facet_grid(m_scale + rho + theta ~ model, labeller = label_both) +
  ggthemes::scale_color_colorblind()

p_good_chains

ggsave(paste0(save_dir, "batchAcceptanceRate.png"),
  plot = p_batch,
  height = 7,
  width = 7
)

ggsave(paste0(save_dir, "groupAcceptanceRate.png"),
  plot = p_group,
  height = 7,
  width = 7
)

ggsave(paste0(save_dir, "completeLikelihood.png"),
  plot = p_likelihood,
  height = 7,
  width = 7
)

ggsave(paste0(save_dir, "completeLikelihoodGoodChains.png"),
  plot = p_good_chains,
  height = 7,
  width = 7
)
