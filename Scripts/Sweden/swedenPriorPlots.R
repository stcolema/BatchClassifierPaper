#
# Plots for the ELISA data from Dopico et al., 2021.
# We plot the
#   * the observed data and the controls
#   * complete log-likelihood to check within/across chain convergence
#   * the inferred and observed data
#   * allocation probability proportion over time
#
# === Setup ====================================================================
library(BatchMixtureModel)
library(mdiHelpR)
library(ggplot2)
library(magrittr)
library(tidyr)
library(mcclust)
library(tibble)
library(dplyr)
library(stringr)
library(tibble)
library(patchwork)

calcSigmaBar <- function(X) {
  mean(diag(cov(X)))
}

empiricalPsi <- function(X, K) {
  P <- ncol(X)

  sigma_bar <- calcSigmaBar(X)

  psi <- diag(sigma_bar / (K^(2 / P)), nrow = P)
  psi
}

groupHyperparameters <- function(X, K, kappa = 0.01) {

  # kappa <- kappa
  nu <- ncol(X) + 2
  mu_0 <- colMeans(X)
  psi <- empiricalPsi(X, K)

  list(
    kappa = kappa,
    nu = nu,
    psi = psi,
    mu_0 = mu_0
  )
}

batchScale <- function(X, m_scale = 1.0) {
  sigma_bar <- calcSigmaBar(X)
  batch_prior_sd <- sigma_bar * m_scale
  batch_prior_sd
}


confidenceEllipse <- function(mu = c(0, 0), Sigma = matrix(c(1, 0, 0, 1), 2, 2), confidenceLevel = 0.95) {
  radius <- sqrt(2 * stats::qf(confidenceLevel, 2, Inf))
  chol_decomp <- chol(Sigma)
  angles <- (0:100) * 2 * pi / 100
  unit.circle <- cbind(cos(angles), sin(angles))
  ellipse <- t(mu + radius * t(unit.circle %*% chol_decomp))
  colnames(ellipse) <- c("X1", "X2")
  as.data.frame(ellipse)
}

setMyTheme()
set.seed(1)

# Data lives in the repo:
# "https://github.com/chr1swallace/seroprevalence-paper/blob/master/adjusted-data.RData"
# Please download before proceeding if it is not already in the Data directory.
# Read in the ELISA data
data_file <- "./Data/Sweden/adjusted-data.RData"

load(data_file)

save_dir <- "./Analysis/Sweden/Outputs/Plots/"

# These samples are poorly behaved, drop
drop_sample_in_12 <- which((m$Sample.ID) %in% 1:2)
m <- m[-drop_sample_in_12, ]

# Drop the "patient 4" samples
patient_4 <- which(m$type == "Patient 4")
m <- m[-patient_4, ]

# Find the controls
negative_controls <- which(m$type == "Historical controls")
positive_controls <- which((m$type == "COVID") | (m$type == "Patient 4"))
controls <- c(negative_controls, positive_controls)

N <- nrow(m)
N_positive <- length(positive_controls)
N_negative <- length(negative_controls)
N_controls <- length(controls)

# The non-control data is used when estimating seroprevalence
non_control_data <- m[-controls, ]

# Find which observations have a known label
fixed <- rep(0, N)
fixed[c(positive_controls, negative_controls)] <- 1

# Label
truth <- m$type %>%
  factor(levels = unique(m$type)) %>%
  as.numeric()

truth <- rep(2, N)
truth[negative_controls] <- 0
truth[positive_controls] <- 1

# Our batch variable
batches_present <- m$group %>%
  unique()

batch_vec <- m$group %>%
  as.factor() %>%
  as.numeric()

# The numeric data to model is contained in the features ``SPIKE`` and ``RBD``.
cols_of_interest <- c("SPIKE", "RBD")
cols_used <- which(colnames(m) %in% cols_of_interest)

# Log transform the data when in matrix form
X <- log(as.matrix(m[, ..cols_used])) %>%
  set_rownames(paste0("Person_", 1:N)) %>%
  set_colnames(cols_of_interest)

# Let's plot the original dataset
plot_df <- X %>%
  as_tibble() %>%
  add_column(
    Label = truth,
    Batch = batch_vec,
    Fixed = fixed
  )

K <- 2


n_draws <- 1e6


P <- ncol(X)

group_hypers <- groupHyperparameters(X, K)

expected_cov <- group_hypers$psi / (group_hypers$nu - P - 1)

expected_mean_distn <- mvtnorm::rmvnorm(n_draws,
  mean = group_hypers$mu_0,
  sigma = expected_cov / group_hypers$kappa
) %>%
  as.data.frame() %>%
  set_colnames(colnames(X))

prior_ellipse <- confidenceEllipse(
  mu = group_hypers$mu_0,
  Sigma = expected_cov / group_hypers$kappa
)

plot_df %>%
  ggplot(aes(x = SPIKE, y = RBD)) +
  geom_point() +
  geom_path(data = prior_ellipse, mapping = aes(x = X1, y = X2))
# +
#   stat_density_2d(data = expected_mean_distn, aes(fill = ..level..), geom = "polygon")


# === Prior hyperparameters ====================================================

# This sets the scale of the batch mean effect
m_scale <- c(1, 0.1, 0.01)
n_m_scales <- length(m_scale)

# Choice of rho and theta to have a common expectation
rho <- c(2, 20, 10) + 1
theta <- c(1, 10, 5)

# Quick check that the rhos and thetas align
n_theta <- length(theta)
n_rho <- length(rho)
rhos_match_thetas <- n_rho == n_theta
if (!rhos_match_thetas) {
  stop("Number of rhos to try and number of thetas to try not matching.")
}

hyper_params <- data.frame(
  m_scale = rep(m_scale, n_theta),
  rho = rep(rho, each = n_m_scales),
  theta = rep(theta, each = n_m_scales)
)

# How many combinations
n_hyperparam_combinations <- nrow(hyper_params)


batch_bias_prior_sd <- m_scale %>%
  sapply(function(x) {
    batchScale(X, x)
  })

prior_distns <- batch_bias_prior_sd %>%
  lapply(function(x) {
    rnorm(n_draws, mean = 0, sd = x)
  })

do.call(cbind, prior_distns) %>%
  data.frame() %>%
  set_colnames(c("1e-2", "1e-1", "1e-0")) %>%
  pivot_longer(everything()) %>%
  ggplot((aes(x = value))) +
  geom_density() +
  facet_wrap(~name, ncol = 1, scales = "free_y")


batch_prior_scales <- list()
for (i in seq.int(n_theta)) {
  batch_prior_scales[[i]] <- 1 + 1.0 / rgamma(n_draws, shape = rho[i], rate = theta[i])
}

p_S <- do.call(cbind, batch_prior_scales) %>%
  data.frame() %>%
  set_colnames(c(
    "shape = 3\nscale = 1",
    "shape = 21\nscale = 10",
    "shape = 11\nscale = 5"
  )) %>%
  pivot_longer(everything()) %>%
  ggplot((aes(x = value))) +
  geom_histogram(bins = 75) +
  facet_wrap(~name, ncol = 1, labeller = label_bquote(rows = .(name))) + # , scales = "free") +
  xlim(c(1, 4)) +
  labs(
    x = "Batch scaling effect",
    y = "Count"
  )

p_m <- do.call(cbind, prior_distns) %>%
  data.frame() %>%
  set_colnames(c(
    "1.00",
    "0.10",
    "0.01"
  )) %>%
  pivot_longer(everything()) %>%
  ggplot((aes(x = value))) +
  geom_histogram(bins = 75) +
  facet_wrap(~name,
    ncol = 1,
    scales = "free_x",
    labeller = label_bquote(rows = lambda ~ "=" ~ .(name))
  ) +
  labs(
    x = "Batch bias-effect",
    y = "Count"
  )



p_S <- do.call(cbind, batch_prior_scales) %>%
  data.frame() %>%
  set_colnames(c(
    "3",
    "21",
    "11"
  )) %>%
  pivot_longer(everything(), names_to = "Shape") %>%
  mutate(Scale = ifelse(Shape == "3", 1, ifelse(Shape == "11", 5, 10))) %>%
  ggplot((aes(x = value))) +
  geom_histogram(bins = 75) +
  facet_wrap(~ Shape + Scale,
    ncol = 1,
    labeller = label_bquote(rows = alpha ~ "=" ~ .(Shape) ~ "," ~ beta ~ "=" ~ .(Scale))
  ) + # , scales = "free") +
  xlim(c(1, 4)) +
  labs(
    x = "Batch scaling effect",
    y = "Count"
  )

p_patchwork <- p_S + p_m

ggsave(paste0(save_dir, "priorDrawsBatchEffects.png"),
  plot = p_patchwork,
  height = 5,
  width = 5
)
