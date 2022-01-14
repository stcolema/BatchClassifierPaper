# acceptance.R
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

# The working directory rellative to the project
my_wd <- "./"

# Save the generated files to
save_path <- paste0(my_wd, "Simulations/Output/" )

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


# === Acceptance ===============================================================

acceptance_full <- read.csv(
  paste0(my_wd, "Simulations/Output/SimAcceptanceRates.csv")
)

# Pivot into a longer form, drop the phi parameters (we don't use the MSN)
long_acceptance_full <- acceptance_full %>%
  select(-c(phi_1, phi_2)) %>%
  pivot_longer(-c(Model, Chain, Scenario, Simulation),
               values_to = "Acceptance_rate", 
               names_to = "Parameter"
               )

# Should be the number of scenarios x number of saved samples x number of clusters
# (i.e. )
MVN_df <- which(is.na(long_acceptance_full$Acceptance_rate))

# Box plot
p_acceptance <- long_acceptance_full[-MVN_df, ] %>%
  ggplot(aes(x = Parameter, y = Acceptance_rate)) +
  geom_boxplot(fill = "gold") +
  facet_grid(Scenario ~ Model, 
    scales = "free_x",
    labeller = labeller(
     Model = label_both,
     Scenario = scenario_labels
    )
  ) + 
  labs(title = "Parameter acceptance rates for Bayesian models",
       subtitle = "All simulations and chains pooled",
       y = "Acceptance rate"
  ) +
  theme(axis.text.x = element_text(angle = 30)) +
  ylim(0, 1)

ggsave(paste0(my_wd, "Simulations/Output/SimAcceptanceRates.png"),
  plot = p_acceptance,
  height = 8, 
  width = 7
)
