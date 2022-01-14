# gewekeConvergence.R
# Using the output of ``modelling.R`` check within-chain convergence using the
# Geweke statistic for the complete log-likelihood.
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

# The path to the project directory
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

# === Geweke within-chain ======================================================

# load("./Simulations/Output/ConvergenceList.Rdata")
convergence_lst <- readRDS(file = paste0(my_wd, "Simulations/Output/ConvergenceList.rds"))

dir.create(paste0(my_wd, "Simulations/Output/Convergence/"))

mcmc_lst <- convergence_lst %>%
  lapply(function(x) {

    # We use the likelihoods to decide if the chains are converged
    cols_used <- c("Complete_likelihood", "Observed_likelihood")

    y <- as.mcmc(x[-c(1:eff_burn), colnames(x) %in% cols_used])

    y
  })

geweke_lst <- mcmc_lst %>%
  # lapply(gewekePlot)
  lapply(function(x) {
    frac_1 <- 0.1
    frac_2 <- 0.5
    n_bins <- 20
    p_value_threshold <- 0.05
    threshold_line_colour <- "grey"
    plt_title <- "Geweke diagnostic plot"

    # The preferred object type for interacting with coda functions
    x <- coda::as.mcmc.list(x)

    # The vector of start iterations to calculate the Geweke statistic for
    start_iter_vec <- floor(seq(
      from = stats::start(x),
      to = (stats::start(x) + stats::end(x)) / 2,
      length = n_bins
    ))

    # The matrix that will hold the Geweke stat
    geweke_mat <- matrix(nrow = length(start_iter_vec), ncol = coda::nvar(x), dimnames = list(start_iter_vec, coda::varnames(x)))

    for (n in 1:length(start_iter_vec)) {
      curr_geweke_diag <- coda::geweke.diag(stats::window(x, start = start_iter_vec[n]),
        frac1 = frac_1,
        frac2 = frac_2
      )

      geweke_mat[n, ] <- curr_geweke_diag[[1]]$z
    }

    # The 1.96 threshold for 0.05 significance on a standard normal distribution
    c_limit <- stats::qnorm(1 - p_value_threshold / 2)

    # The variables to gather when moving from wide to long data (these are our
    # parameters)
    vars_to_gather <- coda::varnames(x)

    # The data.frame we will plot (transform to long data to use the ggplot2
    # framework)
    geweke_df <- data.frame(Start_iteration = start_iter_vec) %>%
      cbind(geweke_mat) %>%
      tidyr::gather_("Parameter", "Geweke_statistic", vars_to_gather)

    geweke_df
  })

# geweke_df <- do.call(rbind, geweke_lst)

for (i in 1:length(geweke_lst)) {
  curr_chain <- convergence_lst[[i]]$Chain %>% unique()
  curr_sim <- convergence_lst[[i]]$Simulation %>% unique()
  curr_scn <- convergence_lst[[i]]$Scenario %>% unique()
  curr_model <- convergence_lst[[i]]$Model %>% unique()
  if (length(curr_chain) > 1) {
    stop("Too many chains.")
  }
  if (length(curr_sim) > 1) {
    stop("Too many sims")
  }
  if (length(curr_scn) > 1) {
    stop("Too many scenarios")
  }
  if (length(curr_model) > 1) {
    stop("Too many models.")
  }
  geweke_lst[[i]]$Chain <- curr_chain # (1:n_chains)[i]
  # geweke_lst[[i]]$Iteration <- convergence_lst[[i]]$Iteration
  geweke_lst[[i]]$Simulation <- curr_sim
  geweke_lst[[i]]$Scenario <- curr_scn
  geweke_lst[[i]]$Model <- curr_model

  if (i == 1) {
    geweke_df <- geweke_lst[[i]]
  } else {
    geweke_df <- rbind(geweke_df, geweke_lst[[i]])
  }
}

geweke_df$Start_iteration <- geweke_df$Start_iteration * thin + burn


# Does the distribution of the Geweke statistics pass the Shapiro Wilks test of
# normality
passed_geweke_table_distn <- geweke_df %>%
  group_by(Chain, Simulation, Scenario, Model) %>%
  summarise(Kept = shapiro.test(Geweke_statistic)$p.value > 0.05)

passed_geweke_table <- geweke_df %>% 
  filter(Start_iteration == 7525, Parameter == "Complete_likelihood") %>% 
  group_by(Chain, Simulation, Model, Scenario) %>% 
  mutate(p.value = 2 * pnorm(-abs(-Geweke_statistic)), Kept = p.value > 0.05)

# Left join the ``Kept`` variable to the main df
geweke_df_updated <- geweke_df %>% 
  filter(Start_iteration == 7525, Parameter == "Complete_likelihood") %>% 
  left_join(passed_geweke_table) %>% 
  select(-Start_iteration, -Parameter) 

data.table::fwrite(geweke_df_updated, 
  paste0(my_wd, "./Simulations/Output/Convergence/GewekeDF.csv")
)

plot_df <- left_join(geweke_df, passed_geweke_table_distn)


# The 1.96 threshold for 0.05 significance on a standard normal distribution
c_limit <- stats::qnorm(1 - 0.05 / 2)

params_used <- c("Complete_likelihood") # unique(geweke_df$Parameter)[35:74]

# Plot the kept chains
p_geweke_dropped <- plot_df %>%
  filter(Parameter %in% params_used, Kept) %>%
  ggplot(aes(
    x = Start_iteration,
    y = Geweke_statistic,
    # colour = as.factor(Chain),
    group = interaction(Chain, Scenario, Simulation, Model)
  )) +
  geom_line() +
  geom_hline(yintercept = c_limit, linetype = "dashed", color = "red") +
  geom_hline(yintercept = -c_limit, linetype = "dashed", color = "red") +
  # facet_grid(Simulation ~ Scenario + Model, labeller = label_both) +
  labs(
    title = "Geweke statistics for the complete likelihood",
    subtitle = "Non-converged chains dropped",
    x = "Start iteration",
    y = "Geweke statistic"
  )

# Plot the kept chains
p_geweke_included <- plot_df %>%
  filter(Parameter %in% params_used) %>%
  ggplot(aes(
    x = Start_iteration,
    y = Geweke_statistic,
    # colour = as.factor(Chain),
    # group = interaction(Chain, Parameter),
    group = interaction(Chain, Scenario, Simulation, Model)
  )) +
  geom_line() +
  geom_hline(yintercept = c_limit, linetype = "dashed", color = "red") +
  geom_hline(yintercept = -c_limit, linetype = "dashed", color = "red") +
  # facet_grid(Simulation ~ Scenario + Model, labeller = label_both) +
  labs(
    title = "Geweke statistics for the complete likelihood",
    subtitle = "Non-converged chains included",
    x = "Start iteration",
    y = "Geweke statistic"
  )


ggsave(paste0(my_wd, "Simulations/Output/Convergence/SimGewekeDropped.png"),
  plot = p_geweke_dropped,
  width = 7,
  height = 7
)
ggsave("./Simulations/Output/Convergence/SimGewekeAll.png",
  plot = p_geweke_included,
  width = 7,
  height = 7
)
