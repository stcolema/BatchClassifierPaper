
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

mdiHelpR::setMyTheme()

dataset <- "Sweden"
burn <- 20000
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


p_likelihood <- likelihood_df %>%
  filter(param == "complete", Iteration > burn) %>%
  ggplot(aes(x = Iteration, y = value, colour = factor(Chain))) +
  geom_line() +
  facet_grid(m_scale + rho + theta ~ model, labeller = label_both) +
  ggthemes::scale_color_colorblind()


p_good_chains <- likelihood_df %>%
  filter(param == "complete", Iteration > burn, Chain %in% c(3, 5)) %>%
  ggplot(aes(x = Iteration, y = value, colour = factor(Chain))) +
  geom_line() +
  facet_grid(m_scale + rho + theta ~ model, labeller = label_both) +
  ggthemes::scale_color_colorblind()


p_mvt_complete_likelihood_all <-
  likelihood_df %>%
  mutate(lambda = m_scale, alpha = rho, beta = theta) %>% 
  filter(param == "complete", Iteration > burn, model == "mvt") %>%
  ggplot(aes(x = Iteration, y = value, colour = factor(Chain))) +
  geom_line() +
  facet_grid(alpha + beta ~ lambda, 
             # labeller = label_bquote(rows = alpha ~ "=" ~ .(alpha)))
             labeller = label_bquote(
               cols = lambda ~ "=" ~ .(lambda),
               rows = alpha ~ "=" ~ .(alpha) ~ "," ~ beta ~ "=" ~ .(beta))
  ) +
  labs(colour = "Chain", y = "Complete log-likelihood") +
  ggthemes::scale_color_colorblind()


p_mvt_complete_likelihood_good <-
  likelihood_df %>%
  mutate(lambda = m_scale, alpha = rho, beta = theta) %>% 
  filter(param == "complete", Iteration > burn, Chain %in% c(1, 3, 4, 5), model == "mvt") %>%
  ggplot(aes(x = Iteration, y = value, colour = factor(Chain))) +
  geom_line() +
  facet_grid(alpha + beta ~ lambda, 
    # labeller = label_bquote(rows = alpha ~ "=" ~ .(alpha)))
    labeller = label_bquote(
      cols = lambda ~ "=" ~ .(lambda),
      rows = alpha ~ "=" ~ .(alpha) ~ "," ~ beta ~ "=" ~ .(beta))
    ) +
  labs(colour = "Chain", y = "Complete log-likelihood") +
  ggthemes::scale_color_colorblind()

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

ggsave(paste0(save_dir, "MVTCompleteLikelihoodGoodChains.png"),
       plot = p_mvt_complete_likelihood_good,
       height = 7,
       width = 7
)

ggsave(paste0(save_dir, "MVTCompleteLikelihoodAll.png"),
       plot = p_mvt_complete_likelihood_all,
       height = 7,
       width = 7
)


p_patch <- p_mvt_complete_likelihood_all / p_mvt_complete_likelihood_good +
  plot_annotation(tag_levels = "A")

ggsave(paste0(save_dir, "MVTCompleteLikelihoodAllandGood.png"),
       plot = p_patch,
       height = 10,
       width = 7
)
