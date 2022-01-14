#!/usr/bin/Rscript
#
#
#
#
#

# === Setup ====================================================================

# install and call various packages
cran_packages <- c(
  "foreach",
  "tidyr",
  "ggplot2",
  "dplyr",
  "magrittr",
  "devtools",
  "tibble"
)

install.packages(setdiff(cran_packages, rownames(installed.packages())))

github_packages <- c(
  "mdiHelpR",
  "BatchMixtureModel"
)

github_repos <- c(
  "stcolema/",
  "stcolema/"
)

for (i in 1:length(github_packages)) {
  pkg <- github_packages[i]
  if (!pkg %in% rownames(installed.packages())) {
    devtools::install_github(paste0(github_repos[i], pkg))
  }
}

library(foreach)
library(tidyr)
library(ggplot2)
library(dplyr)
library(magrittr)
library(tibble)
library(mdiHelpR)
library(BatchMixtureModel)

# Set the gglot2 theme
setMyTheme()

# Set the random seed
set.seed(1)

# Save the generated files to
save_path <- "./Simulations/Generated_data/"

dir.create(save_path, showWarnings = F)

# === Simulation scenarios =====================================================

# The base case parameters
n <- 500
k <- 2
b <- 5
p <- 2
mu_k <- c(-2, 2)
sigma_k <- rep(1.25, k) # rep(2, k)
pi_k <- c(0.75, 0.25) # rep(1 / k, k)
m_b <- rep(0.5, b) * (-1)^(1:b)
s_b <- c(1.2, 1.5, 1.2, 1.5, 1.2) # rep(1.4, b)
# class_weights <- rep(1 / k, k)
batch_weights <- rep(1 / b, b)
frac_known <- 0.2
dof <- NA
type <- "MVN"



scn_df <- tibble(
  Scenario = c("base_case", "no_batch_effects", "varying_batch_size", "varying_batch_effect", "varying_batch_representation", "MVT"),
  N = c(n, n, n, n, n, n),
  K = c(k, k, k, k, k, k),
  B = c(b, b, b, b, b, b),
  P = c(p, p, p, p, p, p),
  mu_k = list(mu_k, mu_k, mu_k, mu_k, mu_k, mu_k),
  sigma_k = list(sigma_k, sigma_k, sigma_k, sigma_k, sigma_k, sigma_k),
  pi_k = list(pi_k, pi_k, pi_k, pi_k, matrix(c(0.7, 0.8, 0.5, 0.2, 0.1, 0.3, 0.2, 0.5, 0.8, 0.9), nrow = 2, byrow = T), pi_k),
  m_b = list(m_b, rep(0, 5), m_b, c(-1.5, -0.5, 0.0, 0.5, 1.5), m_b, m_b),
  S_b = list(s_b, rep(1, 5), s_b, c(1.0, 1.25, 1.5, 1.75, 2.25), s_b, s_b),
  batch_weights = list(batch_weights, batch_weights, c(1 / 2, 1 / 4, 1 / 8, 1 / 16, 1 / 16), batch_weights, batch_weights, batch_weights),
  frac_known = c(frac_known, frac_known, frac_known, frac_known, frac_known, frac_known),
  dof = list(dof, dof, dof, dof, dof, c(3, 5)),
  type = c(type, type, type, type, type, "MVT"),
)

saveRDS(scn_df, file = paste0(save_path, "scenario_descriptions.rds"))

n_scn <- nrow(scn_df)
n_sim <- 10

for (ii in 1:n_scn) {
  
  scn <- scn_df$Scenario[ii]
  N <- scn_df$N[ii]
  P <- scn_df$P[ii]
  cluster_means <- scn_df$mu_k[[ii]]
  std_dev <- scn_df$sigma_k[[ii]]
  batch_shift <- scn_df$m_b[[ii]]
  batch_var <- scn_df$S_b[[ii]]
  cluster_weights <- scn_df$pi_k[[ii]]
  batch_weights <- scn_df$batch_weights[[ii]]
  frac_known <- scn_df$frac_known[ii]
  type <- scn_df$type[ii]
  df <- scn_df$dof[[ii]]
  
  # The directory for the current scenario
  curr_dir <- paste0(save_path, scn)
  dir.create(curr_dir, showWarnings = F)
  
  for (jj in 1:n_sim) {
    set.seed(jj)
    generated <- FALSE

    if (scn == "varying_batch_representation") {
      my_data <- generateBatchDataVaryingRepresentation(
        N,
        P,
        cluster_means,
        std_dev,
        batch_shift,
        batch_var,
        cluster_weights,
        batch_weights,
        frac_known
      )
      generated <- TRUE
    }
    if (type == "MVT") {
      my_data <- generateBatchDataMVT(N,
        P,
        cluster_means,
        std_dev,
        batch_shift,
        batch_var,
        cluster_weights,
        batch_weights,
        frac_known,
        df = df
      )

      generated <- TRUE
    }
    if (!generated) {
      my_data <- BatchMixtureModel::generateBatchData(
        N,
        P,
        cluster_means,
        std_dev,
        batch_shift,
        batch_var,
        cluster_weights,
        batch_weights,
        frac_known
      )
    }
    
    my_data <- tibble(
      labels = factor(my_data$group_IDs),
      batch = factor(my_data$batch_IDs),
      fixed = my_data$fixed,
      X_observed = my_data$observed_data[, 1],
      Y_observed = my_data$observed_data[, 2],
      X_true = my_data$corrected_data[, 1],
      Y_true = my_data$corrected_data[, 2]
    )
    
    write.csv(my_data, file = paste0(curr_dir, "/seed_", jj, ".csv"), row.names = F)
    saveRDS(my_data, file = paste0(curr_dir, "/seed_", jj, ".rds"))
  }
}
