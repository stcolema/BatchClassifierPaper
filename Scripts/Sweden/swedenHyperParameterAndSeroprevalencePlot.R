#
# Plots for the hyperparameter choices for the  ELISA data from Dopico et al., 2021.
# We plot the
#   * the draws from the prior distributions for the batch effects with different
#     hyperparameter choices
#   * the seroprevalence estimate
#
# === Packages =================================================================
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
library(lubridate)

# === Functions ================================================================

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

#' @title Samples to data frame
#' @description Turns the output from the mixture models into a data.frame
#' @param samples Output from the ``batchSemiSupervisedMixtureModel`` or
#' ``batchMixtureModel``.
#' @param type The type of mixture model used; this changes which parameters
#' the function expects to find.
#' @param R The number of iterations run. Defaults to the number of slices in
#' the sampled batch mean array.
#' @param thin The thinning factor of the sampler. Defaults to 1.
#' @param keep_allocation A logical indicating if the final data frame should
#' include the sampled class/cluster membership variables.
#' @return A wide data.frame of all the sampled parameters and the iteration.
#' @examples
#' # Data in matrix format
#' X <- matrix(c(rnorm(100, 0, 1), rnorm(100, 3, 1)), ncol = 2, byrow = TRUE)
#'
#' # Observed batches represented by integers
#' batch_vec <- sample(seq(1, 5), size = 100, replace = TRUE)
#'
#' # MCMC iterations (this is too low for real use)
#' R <- 100
#' thin <- 5
#'
#' # MCMC samples
#' samples <- batchMixtureModel(X, R, thin, batch_vec, "MVN")
#'
#' samples_df <- samplesToDF(samples, "MVN", R = R, thin = thin)
#' @importFrom stringr str_match
#' @export
samplesToDF <- function(samples, type,
                        R = nrow(samples$samples),
                        thin = 1,
                        keep_allocation = TRUE) {
  
  # Number of classes and batches
  K <- ncol(samples$means[, , 1])
  B <- ncol(samples$batch_shift[, , 1])
  P <- nrow(samples$means[, , 1])
  N <- ncol(samples$samples)
  
  # Stack the sampled matrices on top of each other
  means_df <- data.frame(t(apply(samples$means, 3L, rbind)))
  batch_shift_df <- data.frame(t(apply(samples$batch_shift, 3L, rbind)))
  batch_scale_df <- data.frame(t(apply(samples$batch_scale, 3L, rbind)))
  mean_sums_df <- data.frame(t(apply(samples$mean_sum, 3L, rbind)))
  
  # Indices over columns and batches
  col_inds <- seq(1, P)
  batch_inds <- seq(1, B)
  group_inds <- seq(1, K)
  
  # Give sensible column names
  colnames(means_df) <- suppressWarnings(
    paste0(
      "Mu_",
      sort(as.numeric(levels(interaction(group_inds, col_inds, sep = ""))))
    )
  )
  
  colnames(batch_shift_df) <- suppressWarnings(
    paste0(
      "m_",
      sort(as.numeric(levels(interaction(batch_inds, col_inds, sep = ""))))
    )
  )
  
  colnames(batch_scale_df) <- suppressWarnings(
    paste0(
      "S_",
      sort(as.numeric(levels(interaction(batch_inds, col_inds, sep = ""))))
    )
  )
  
  # The combination objects are awkward to name correctly
  mean_sum_names <- suppressWarnings(levels(interaction(
    colnames(means_df),
    colnames(batch_shift_df)
  )))
  
  mean_sum_names <- as.data.frame(stringr::str_match(
    mean_sum_names,
    "Mu_([:digit:]*).m_([:digit:]*)"
  ))
  
  colnames(mean_sum_names) <- c("Comb", "Mu", "m")
  mean_sum_names$Mu <- as.numeric(mean_sum_names$Mu)
  mean_sum_names$m <- as.numeric(mean_sum_names$m)
  
  correct_comb <- which(mean_sum_names$Mu %% 10 == mean_sum_names$m %% 10)
  mean_sum_names <- mean_sum_names[correct_comb, ]
  
  inds <- matrix(seq(1, (P * B * K)), nrow = 4, byrow = T)
  colnames(mean_sums_df) <- mean_sum_names$Comb[order(mean_sum_names$Mu)][c(inds)]
  
  # The covariance is slightly more awkward
  cov_df <- data.frame(t(apply(samples$covariance, 3L, rbind)))
  colnames(cov_df) <- suppressWarnings(
    paste0(
      "Sigma_",
      sort(
        as.numeric(
          levels(
            interaction(
              list(group_inds, col_inds, col_inds),
              sep = ""
            )
          )
        )
      )
    )
  )
  
  # The combined batch and cluster covariances - this only effects the diagonal
  # entries, so let's keep only them
  cov_comb_df <- data.frame(t(apply(samples$cov_comb, 3L, rbind)))
  cov_comb_df <- cov_comb_df[, seq(1, ncol(cov_comb_df)) %% 4 %in% c(0, 1)]
  
  comb_cov_names <- c(
    suppressWarnings(
      levels(
        interaction(
          colnames(cov_df),
          colnames(batch_scale_df)
        )
      )
    )[(P**2) * seq(0, ((2 * K * B) - 1)) + 1][(seq(1, (2 * K * B)) %% 4) %in% c(1, 2)],
    suppressWarnings(
      levels(
        interaction(
          colnames(cov_df),
          colnames(batch_scale_df)
        )
      )
    )[(P**2) * seq(0, ((2 * K * B) - 1)) + 4][((seq(1, (2 * K * B)) %% 4) %% 4) %in% c(0, 3)]
  )
  
  inds <- matrix(c(
    seq(1, (B * K)),
    seq((B * K + 1), (2 * B * K))
  ),
  nrow = 2,
  byrow = TRUE
  )
  
  colnames(cov_comb_df) <- comb_cov_names[c(inds)]
  
  # The sampled weights
  weights_df <- as.data.frame(samples$weights)
  colnames(weights_df) <- c(paste0("pi_", group_inds))
  
  # Add a variable for the iteration the sample comes from
  Iteration <- seq(1, (R / thin)) * thin
  
  output_df <- as.data.frame(
    cbind(
      Iteration,
      weights_df,
      means_df,
      cov_df,
      batch_shift_df,
      batch_scale_df,
      mean_sums_df,
      cov_comb_df
    )
  )
  
  # Type dependent parameters
  if (type == "MVT") {
    df_df <- samples$t_df
    colnames(df_df) <- paste0("t_df_", group_inds)
    output_df <- cbind(output_df, df_df)
  }
  
  if (type == "MSN") {
    shape_df <- data.frame(t(apply(samples$shapes, 3L, rbind)))
    colnames(shape_df) <- suppressWarnings(
      paste0(
        "phi_",
        sort(as.numeric(levels(interaction(group_inds, col_inds, sep = ""))))
      )
    )
    output_df <- cbind(output_df, shape_df)
  }
  
  # Sampled observed likelihood and BIC
  output_df$Complete_likelihood <- samples$complete_likelihood
  output_df$Observed_likelihood <- samples$observed_likelihood
  output_df$BIC <- samples$BIC
  
  # The sampled allocations
  if (keep_allocation) {
    samples_df <- data.frame(samples$samples)
    colnames(samples_df) <- paste0("c_", seq(1, N))
    output_df <- cbind(output_df, samples_df)
  }
  
  output_df
}

# === Setup ====================================================================


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
positive_controls <- which((m$type == "COVID"))
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

# === MCMC =====================================================================

# MCMC parameters
R <- 50000
thin <- 100
burn <- 20000
eff_burn <- burn / thin

# === Extract data =============================================================

# Find the allocations across time
# The non-control data has the week of the sample collection included in the
# sample ID
weeks_str <- non_control_data$Sample.ID

# Extract the week
weeks <- weeks_str %>%
  str_match_all("Wk(\\d{1,2})") %>%
  do.call(rbind, .) %>%
  `[`(, 2) %>%
  as.numeric()

weeks_present <- unique(weeks)

chain_used <- 3

my_list <- list(1)

# Tibble to hold results
chain_tib <- hyper_params %>%
  as_tibble() %>%
  mutate(
    Chain = chain_used,
    R = R,
    thin = thin,
    burn = burn,
    MVT_samples = my_list,
    MVN_samples = my_list,
    MVN_prob = my_list,
    MVN_pred = my_list,
    MVT_prob = my_list,
    MVT_pred = my_list
  )

chain_dir <- "./Analysis/Sweden/Outputs/Chains/"

for (i in seq.int(n_hyperparam_combinations)) {
  .m_scale <- chain_tib$m_scale[i]
  .rho <- chain_tib$rho[i]
  .theta <- chain_tib$theta[i]

  m_scale_sci_notation <- formatC(.m_scale, format = "e", digits = 0)

  f_mvt <- paste0(
    chain_dir,
    "sweden_",
    "mvt",
    "_chain_",
    chain_used,
    "_m_scale_",
    m_scale_sci_notation,
    "_rho_",
    .rho,
    "_theta_",
    .theta,
    ".rds"
  )

  mvt_samples <- readRDS(f_mvt)

  f_mvn <- paste0(
    chain_dir,
    "/sweden_",
    "mvn",
    "_chain_",
    chain_used,
    "_m_scale_",
    m_scale_sci_notation,
    "_rho_",
    .rho,
    "_theta_",
    .theta,
    ".rds"
  )

  mvn_samples <- readRDS(f_mvn)


  # Allocations
  mvn_prob <- calcAllocProb(mvn_samples$alloc, eff_burn)
  mvn_pred <- predictClass(mvn_prob)

  mvt_prob <- calcAllocProb(mvt_samples$alloc, eff_burn)
  mvt_pred <- predictClass(mvt_prob)

  # 95% credible interval
  mvn_ci_seropositive <- mvn_samples$alloc[, 2, ] %>%
    apply(1, quantile, probs = c(0.025, 0.975)) %>%
    t()

  mvt_ci_seropositive <- mvt_samples$alloc[, 2, ] %>%
    apply(1, quantile, probs = c(0.025, 0.975)) %>%
    t()

  mvn_seropositive <- data.frame(
    "Estimate" = mvn_prob[-controls, 2],
    "CI_0.025" = mvn_ci_seropositive[-controls, 1],
    "CI_0.975" = mvn_ci_seropositive[-controls, 2],
    Week = weeks,
    "Model" = "MVN",
    "Chain" = chain_used,
    m_scale = .m_scale,
    rho = .rho,
    theta = .theta
  )

  mvt_seropositive <- data.frame(
    "Estimate" = mvt_prob[-controls, 2],
    "CI_0.025" = mvt_ci_seropositive[-controls, 1],
    "CI_0.975" = mvt_ci_seropositive[-controls, 2],
    Week = weeks,
    "Model" = "MVT",
    "Chain" = chain_used,
    m_scale = .m_scale,
    rho = .rho,
    theta = .theta
  )

  seropositive_df <- rbind(mvn_seropositive, mvt_seropositive)

  if (i == 1) {
    seroprevalence_prob_df <- seropositive_df
  } else {
    seroprevalence_prob_df <- rbind(seroprevalence_prob_df, seropositive_df)
  }

  # # === Seroprevalence v2
  # sampled_allocs <- mvt_samples$samples[-seq(1, eff_burn), -controls]
  #
  # reduced_batch_vec <- batch_vec[-controls]
  #
  # alloc_df <- sampled_allocs %>%
  #   as_tibble() %>%
  #   set_colnames(1:(N-N_controls)) %>%
  #   add_column(Iteration = seq(burn + thin, R, thin), .before = 1) %>%
  #   pivot_longer(-Iteration, names_to = "Sample_str", values_to = "Label") %>%
  #   mutate(Sample = as.numeric(Sample_str))
  #
  # alloc_df$Week <- weeks[alloc_df$Sample]
  #
  # alloc_df %>%
  #   group_by(Week, Iteration) %>%
  #   summarise(Seroprevalence = sum(Label) / n()) %>%
  #   group_by(Week) %>%
  #   summarise(Estimate = mean(Seroprevalence)) %>%
  #   ggplot(aes(x = Week, y = Estimate)) +
  #   geom_point()
  #
  # sampled_allocs %>%
  #   rowSums() %>%
  #   `/`(N - N_controls) %>%
  #   quantile(probs = c(0.025, 0.5, 0.975))

  # Inferred datasets
  mvn_inferred_data <- rowMeans(mvn_samples$batch_corrected_data[, , -c(1:eff_burn), drop = F], dims = 2) %>%
    as_tibble(.name_repair = "minimal") %>%
    set_colnames(colnames(X)) %>%
    add_column(
      "Label" = factor(mvn_pred),
      "Prob" = apply(mvn_prob, 1, function(x) {
        x[which.max(x)]
      }),
      "Batch" = factor(batch_vec),
      "Type" = "Inferred",
      "Fixed" = factor(fixed),
      "Model" = "MVN",
      "Chain" = chain_used,
      m_scale = .m_scale,
      rho = .rho,
      theta = .theta
    )

  mvt_inferred_data <- rowMeans(mvt_samples$batch_corrected_data[, , -c(1:eff_burn), drop = F], dims = 2) %>%
    as_tibble(.name_repair = "minimal") %>%
    set_colnames(colnames(X)) %>%
    add_column(
      "Label" = factor(mvt_pred),
      "Prob" = apply(mvt_prob, 1, function(x) {
        x[which.max(x)]
      }),
      "Batch" = factor(batch_vec),
      "Type" = "Inferred",
      "Fixed" = factor(fixed),
      "Model" = "MVT",
      "Chain" = chain_used,
      m_scale = .m_scale,
      rho = .rho,
      theta = .theta
    )

  .inf_df <- rbind(mvn_inferred_data, mvt_inferred_data)

  if (i == 1) {
    inferred_data <- .inf_df
  } else {
    inferred_data <- rbind(inferred_data, .inf_df)
  }

  # Associate the allocation probabilities with the week
  .week_df <- data.frame(
    Week = weeks,
    m_scale = .m_scale,
    rho = .rho,
    theta = .theta,
    chain = chain_used,
    MVN_prob = mvn_prob[-controls, 2],
    MVT_prob = mvt_prob[-controls, 2],
    MVN_pred = mvn_pred[-controls] - 1,
    MVT_pred = mvt_pred[-controls] - 1
  )

  if (i == 1) {
    week_prob_df <- .week_df
  } else {
    week_prob_df <- rbind(week_prob_df, .week_df)
  }

  # Sampled parameters
  mvn_sampled_data <- samplesToDF(mvn_samples, "MVN", R, thin,
    keep_allocation = FALSE
  )
  mvn_sampled_data$Model <- "MVN"
  mvn_sampled_data$m_scale <- .m_scale
  mvn_sampled_data$rho <- .rho
  mvn_sampled_data$theta <- .theta
  mvn_sampled_data$chain <- chain_used
  mvn_sampled_data$t_df_1 <- mvn_sampled_data$t_df_2 <- NA

  mvt_sampled_data <- samplesToDF(mvt_samples, "MVT", R, thin,
    keep_allocation = FALSE
  )

  mvt_sampled_data$Model <- "MVT"
  mvt_sampled_data$m_scale <- .m_scale
  mvt_sampled_data$rho <- .rho
  mvt_sampled_data$theta <- .theta
  mvt_sampled_data$chain <- chain_used

  # Bind the sampled parameters together and move the degrees of freedom into
  # the df
  .samples_df <- rbind(mvn_sampled_data, mvt_sampled_data)[, c(1:99, 108:109, 100:107)]

  if (i == 1) {
    sampled_params <- .samples_df
  } else {
    sampled_params <- rbind(sampled_params, .samples_df)
  }

  # Save everything to the tibble
  chain_tib$MVT_samples[[i]] <- mvt_samples
  chain_tib$MVN_samples[[i]] <- mvn_samples

  chain_tib$MVN_prob[[i]] <- mvn_prob
  chain_tib$MVN_pred[[i]] <- mvn_pred

  chain_tib$MVT_prob[[i]] <- mvt_prob
  chain_tib$MVT_pred[[i]] <- mvt_pred
}

# === Plotting =================================================================

## == Seroprevalence ===========================================================

dopico_lda_svm <- tibble::tribble(
  ~Week,       ~N,   ~estimate,   ~lower.ci,  ~upper.ci,
     14,       NA,          NA,          NA,         NA,
     17,     200L,     0.04262,    0.021821,   0.081591,
     18,     200L,     0.04334,    0.022241,   0.082759,
     19,     200L,     0.07868,    0.048365,   0.125491,
     20,     200L,     0.07355,    0.044336,   0.119604,
     21,     200L,    0.074855,    0.045649,    0.12039,
     22,     200L,     0.03541,    0.017197,    0.07151,
     23,     200L,    0.084055,    0.052379,   0.132213,
     24,     200L,    0.068995,    0.041077,   0.113638,
     25,     200L,     0.06435,    0.037035,   0.109519,
     30,     200L,    0.131995,    0.089848,    0.18979,
     31,     200L,     0.09163,    0.058081,   0.141644,
     32,     200L,     0.11301,    0.074876,   0.167059,
     33,     200L,    0.108385,    0.071907,   0.160175,
     34,     200L,    0.129665,     0.08875,     0.1856,
     45,     200L,     0.11786,    0.079438,   0.171406,
     46,     200L,    0.142025,    0.099063,   0.199494,
     47,     200L,     0.13373,    0.091923,   0.190561,
     48,     200L,     0.13202,    0.090481,   0.188674,
     49,     200L,     0.11522,    0.076947,   0.169044,
     50,     200L,    0.157285,    0.112197,   0.216082
  )


# dopico_lda_svm <- tibble::tribble(
#   ~Week,    ~N,          ~estimate,           ~lower.ci,          ~upper.ci,
#      14,  100L, 0.0105599999999998, 0.00184711460040343, 0.0579842820109458,
#      17,  200L, 0.0431150000000002,  0.0221223112761615, 0.0823508337277539,
#      18,  200L, 0.0434700000000001,   0.022299970205302, 0.0830306178613479,
#      19,  200L, 0.0798099999999993,  0.0493467575514476,  0.126575203056695,
#      20,  200L, 0.0689349999999998,  0.0407436413657359,  0.114308116426858,
#      21,  200L, 0.0750750000000008,  0.0457905956485781,  0.120718292288378,
#      22,  200L, 0.0355349999999998,  0.0172721597953066, 0.0716993107746621,
#      23,  200L, 0.0844649999999996,  0.0527466648833784,  0.132586944479438,
#      24,  200L, 0.0694549999999997,  0.0414740782661843,   0.11406705818783,
#      25,  200L,           0.063885,  0.0366526359995393,  0.109060011958329,
#      30,  200L,           0.129955,  0.0881953852243244,  0.187423042342957,
#      31,  200L, 0.0919100000000001,  0.0583328478485414,  0.141902101175647,
#      32,  200L,           0.113775,  0.0755409717516833,  0.167847126993757,
#      33,  200L,            0.10752,  0.0711645331155792,  0.159263610639403,
#      34,  200L,           0.128185,  0.0874676569420377,  0.184034176973193,
#      45,  200L,           0.116765,  0.0785313072613766,  0.170175435429028,
#      46,  200L,            0.14462,   0.101157728048079,  0.202547830993344,
#      47,  200L,            0.13018,  0.0891305166227496,  0.186268917111879,
#      48,  200L,            0.12857,  0.0874735682462468,  0.185058597470984,
#      49,  200L,  0.113939999999999,  0.0760462408645264,  0.167297201965587,
#      50,  200L,            0.15694,   0.111767424356249,  0.215931416896436
#   )

dopico_bayesian_learner <- tibble::tribble(
  ~Week,             ~lower.ci,              ~estimate,           ~upper.ci,
    14L,   0.00333905545967725,       0.02349393230225, 0.0595868242834863,
    15L,    0.0127122783076255,     0.0317909490715836, 0.0608104458643637,
    16L,    0.0208638579090816,     0.0402845670351009, 0.0620504174459884,
    17L,    0.0295963254883042,       0.04729171578509,  0.065716751371419,
    18L,    0.0387240932482158,     0.0534179746962467, 0.0687252026630859,
    19L,    0.0451850018388308,     0.0585891107438571, 0.0728465674055788,
    20L,    0.0492397615814623,   0.0628226635399526, 0.0773572868958386,
    21L,    0.0523082116011915,     0.0663870558115672,  0.081887400709906,
    22L,    0.0550393379747273,     0.0697372994755059, 0.0865371126838118,
    23L,    0.0583468603596189,     0.0730796425250501, 0.0911370012466751,
    24L,    0.0612692386173902,     0.0775280803104066, 0.0968369607150626,
    25L,    0.0642767054550179,      0.082979108525533,  0.103884811898326,
    26L,    0.0676258376814479,     0.0886911427698066,  0.110871494316406,
    27L,    0.0719117007456125,     0.0954836965451289,  0.119849512389092,
    28L,    0.0765892444085248,      0.101902992632051,  0.127624463365591,
    29L,    0.0824869906226206,      0.108193454256367,  0.133409910319569,
    30L,    0.0886654086243632,      0.113420683947117,  0.138571350550605,
    31L,    0.0936423920760741,    0.117905209831753,  0.142341877658316,
    32L,    0.0980017402422175,    0.121467421041104,  0.145819964702898,
    33L,     0.100923405204302,    0.124285872769697,   0.14821128065733,
    34L,     0.103404619396124,     0.12654678244523,  0.150801478887081,
    35L,       0.1049210553012,    0.128706173383708,  0.153196501978474,
    36L,     0.106554430491549,    0.130737925244621,  0.156675972084953,
    37L,      0.10741386384969,    0.132193277360702,  0.158895293237457,
    38L,     0.108914019894821,    0.133553364094392,  0.161266379826124,
    39L,     0.110222876696505,    0.134863182391319,  0.163364192773837,
    40L,     0.111157314459714,    0.136137034981484,  0.163875388530462,
    41L,     0.112004408424312,    0.137414908779515,  0.165346810565715,
    42L,      0.11339854680334,    0.138736090590799,  0.167220015975339,
    43L,     0.114562420918706,    0.140196585774159,  0.167968901233307,
    44L,     0.116374552205785,    0.141504168333636,  0.168706975719988,
    45L,     0.118294644304495,    0.142835866141086,  0.169448290033326,
    46L,     0.121103353362968,    0.144717471907636,  0.170899594285061,
    47L,     0.123457802105182,    0.147206357497034,   0.17363673056285,
    48L,     0.126231772075056,    0.149814712612457,  0.176899754114051,
    49L,     0.128411606492678,    0.152893960773629,  0.181648768886363,
    50L,     0.130739865592731,    0.156403885781036,  0.188660929301118
)




svm_lda_est <- dopico_lda_svm %>%
  mutate(
    Model = "SVM-LDA"
  ) %>%
  # group_by(Week, Model) %>%
  mutate(
    "Seroprevalence" = 1 * estimate,
    "97.5%" = 1 * upper.ci,
    "02.5%" = 1 * lower.ci,
    Name = "Dopico",
    Date = as.Date(paste0("2020", Week, "7"), "%Y%U%u"),
    m_scale = NA,
    rho = NA,
    theta = NA
  ) %>%
  select(-c(estimate, lower.ci, upper.ci))

bayesian_learner_est <- dopico_bayesian_learner %>%
  mutate(
    Model = "Bayesian learner"
  ) %>%
  # group_by(Week, Model) %>%
  mutate(
    "Seroprevalence" = 1 * estimate,
    "97.5%" = 1 * upper.ci,
    "02.5%" = 1 * lower.ci,
    Name = "Dopico",
    Date = as.Date(paste0("2020", Week, "7"), "%Y%U%u"),
    m_scale = NA,
    rho = NA,
    theta = NA
  ) %>%
  select(-c(estimate, lower.ci, upper.ci))


my_sero_df <- seroprevalence_prob_df %>%
  group_by(Week, Model, m_scale, rho, theta) %>%
  summarise(
    N = n(),
    Seroprevalence = sum(Estimate) / N,
    "97.5%" = Seroprevalence + 1.96 * sqrt(Seroprevalence * (1 - Seroprevalence) / N), # 100 * sum(CI_0.975) / n(),
    "02.5%" = max(0, Seroprevalence - 1.96 * sqrt(Seroprevalence * (1 - Seroprevalence) / N)), # 100 * sum(CI_0.025) / n(),
    Name = paste0(Model, "_", m_scale, "_", rho, "_", theta)
  ) %>%
  filter(Model == "MVT") %>%
  mutate(Date = as.Date(paste0("2020", Week, "3"), "%Y%U%u"))

p_sero_ci <- svm_lda_est %>%
  # mutate(Date = as.Date(paste0("2020", Week, "7"), "%Y%U%u")) %>%
  rbind(my_sero_df) %>%
  ggplot(aes(x = Date)) +
  # geom_line(aes(y = Seroprevalence,
  #               colour = factor(Model, labels = c("Dopico et al.", "MVT")),
  #               group = Name), lty = 2) +
  geom_point(aes(
    y = 100 * Seroprevalence,
    colour = factor(Model), #, labels = c("Dopico et al.", "MVT")),
    group = Name
  )) +
  geom_errorbar(aes(
    ymin = 100 * `02.5%`,
    ymax = 100 * `97.5%`,
    color = factor(Model), # , labels = c("Dopico et al.", "MVT")),
    group = Name
  ),
  width = 3
  ) +
  # geom_errorbar(data = my_sero_df,
  # mapping = aes(ymin=`02.5%`, ymax=`97.5%`), width=3, colour = "#00BFC4") +
  # geom_line(aes(y = `97.5%`,
  #               colour = factor(Model, labels = c("Dopico et al.", "MVT")),
  #               group = Name), lty = 2) +
  # geom_line(aes(y = `02.5%`,
  #               colour = factor(Model, labels = c("Dopico et al.", "MVT")),
  #               group = Name), lty = 2) +
  labs(
    y = "Seroprevalence (%)",
    colour = "Model:"
  ) +
  theme(legend.position = "bottom")

p_sero_ci <- p_sero_ci +
  geom_ribbon(data = bayesian_learner_est, 
    mapping = aes(ymin = 100 * `02.5%`, ymax = 100 * `97.5%`), 
    alpha = 0.3
  ) +
  geom_line(data = bayesian_learner_est, 
    mapping = aes(y = 100 * Seroprevalence, colour = Model)
  ) +
  ggthemes::scale_color_colorblind()
  
# seroprevalence_df <- week_prob_df %>%
#   select(-c(MVT_pred, MVN_pred)) %>%
#   pivot_longer(c(MVT_prob, MVN_prob),
#     names_sep = "_", #
#     names_to = c("Model", ".value")
#   ) %>%
#   group_by(Week, Model, m_scale, rho, theta) %>%
#   summarise(
#     Num_seropositive = sum(prob) / n(),
#     Seroprevalence = sum(prob) / n()
#   )
#
# sero_plot_df <- seroprevalence_df %>%
#   filter(Model == "MVT") %>%
#   select(-Num_seropositive) %>%
#   pivot_wider(
#     names_from = c(Model, m_scale, rho, theta),
#     values_from = Seroprevalence
#   ) %>%
#   add_column(Dopico = dopico_pred$Bayesian_est / 100) %>%
#   pivot_longer(-Week) %>%
#   mutate(Model = ifelse(grepl("MVN", name),
#     "MVN",
#     ifelse(grepl("MVT", name),
#       "MVT",
#       "Dopico et al"
#     )
#   ))
#
#
# p_sero <- sero_plot_df %>%
#   mutate(Date = as.Date(paste0("2020", Week, "7"), "%Y%U%u")) %>%
#   ggplot() +
#   geom_point(aes(x = Date, y = 100 * value, colour = Model, group = name)) +
#   geom_line(aes(x = Date, y = 100 * value, colour = Model, group = name), lty = 2) +
#   labs(
#     y = "Seroprevalence (%)",
#     colour = "Model:"
#   ) +
#   theme(legend.position = "bottom")
#
# ggsave(paste0(save_dir, "seroprevalence.png"),
#   height = 5,
#   width = 5,
#   plot = p_sero
# )
#
#
# # Add in national covid data
# # National population from https://www.worldometers.info/world-population/sweden-population/
# sweden_pop <- 10170753
#
# # Our world in data
# # Citation: https://ourworldindata.org/coronavirus/country/sweden#citation
# owid_data <- read.csv("./Data/Sweden/owid-sweden-covid-data.csv") %>%
#   select(date, total_cases, new_cases, new_deaths, total_deaths) %>%
#   mutate(cases_less_deaths = total_cases - total_deaths,
#          seroprevalence = cases_less_deaths / sweden_pop)
#
# # Dates covered by samples
# min_date <- dmy("30/03/2020")
# max_date <- dmy("13/12/2020")
#
# rel_owid_data <- owid_data %>%
#   mutate(date = dmy(date)) %>%
#   filter(date >= min_date, date <= max_date)
#
# p_sero_w_national_data <- p_sero +
#   geom_point(data = rel_owid_data, aes(x = date, y = 700 * seroprevalence)) +
#   scale_y_continuous(
#     # Add a second axis and specify its features
#     sec.axis = sec_axis(~./7, name="National confirmed cases less COVID-related deaths")
#   )
#
# ggsave(paste0(save_dir, "seroprevalenceWithNationalData.png"),
#        height = 5,
#        width = 5,
#        plot = p_sero_w_national_data
# )
#
# p_sero_ci_w_national_data <- p_sero_ci +
#   geom_point(data = rel_owid_data,
#     aes(x = date, y = 700 * seroprevalence),
#     colour = "black",
#     size = 1.2
#   ) +
#   # geom_line(data = rel_owid_data,
#   #            aes(x = date, y = 700 * seroprevalence),
#   #            colour = "black"
#   # ) +
#   scale_y_continuous(
#     # Add a second axis and specify its features
#     sec.axis = sec_axis(~./7, name="National confirmed cases less COVID-related deaths")
#   )
#
# ggsave(paste0(save_dir, "seroprevalenceWithCIAndNationalData.png"),
#        height = 5,
#        width = 5,
#        plot = p_sero_ci_w_national_data
# )

# === Prior distributions ======================================================

K <- 2
n_draws <- 1e6
P <- ncol(X)

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

batch_prior_scales <- list()
for (i in seq.int(n_theta)) {
  batch_prior_scales[[i]] <- 1 + 1.0 / rgamma(n_draws, shape = rho[i], rate = theta[i])
}

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
    # scales = "free_y",
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


# p_patchwork <- (p_S + p_m) / p_sero +
#   plot_annotation(tag_levels = "A")
#
# ggsave(paste0(save_dir, "hyperparameterPriorsAndSeroprevalence.png"),
#   height = 8,
#   width = 6,
#   plot = p_patchwork
# )

p_patchwork_ci <- (p_S + p_m) / p_sero_ci +
  plot_annotation(tag_levels = "A")

ggsave(paste0(save_dir, "hyperparameterPriorsAndSeroprevalenceCI.png"),
  height = 8,
  width = 8,
  plot = p_patchwork_ci
)

# p_patchwork_w_national_data <- (p_S + p_m) / p_sero_w_national_data +
#   plot_annotation(tag_levels = "A")
#
# ggsave(paste0(save_dir, "hyperparameterPriorsAndSeroprevalenceWNationalData.png"),
#        height = 8,
#        width = 6,
#        plot = p_patchwork_w_national_data
# )
#
# p_patchwork_w_ci_and_national_data <- (p_S + p_m) / p_sero_ci_w_national_data +
#   plot_annotation(tag_levels = "A")
#
# ggsave(paste0(save_dir, "hyperparameterPriorsAndSeroprevalenceWCIsNationalData.png"),
#        height = 8,
#        width = 6,
#        plot = p_patchwork_w_ci_and_national_data
# )
