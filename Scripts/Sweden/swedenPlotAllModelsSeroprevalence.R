#
# Plots for the hyperparameter choices for the  ELISA data from Dopico et al., 2021.
# We plot the
#   * the draws from the prior distributions for the batch effects with different
#     hyperparameter choices
#   * the seroprevalence estimate
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
library(lubridate)

setMyTheme()
set.seed(1)

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

my_wd <- "./"

output_df <- read.csv(paste0(my_wd, "/Analysis/Sweden/Outputs/otherModelsSeroprevalence.csv"))


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
    Name = "SVM-LDA",
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
    Name = "Bayesian learner",
    Date = as.Date(paste0("2020", Week, "7"), "%Y%U%u"),
    m_scale = NA,
    rho = NA,
    theta = NA
  ) %>%
  select(-c(estimate, lower.ci, upper.ci))
# 
# dopico_pred <- tibble::tribble(
#   ~type,      ~Week,    ~N,          ~estimate,           ~lower.ci,          ~upper.ci,
#   "Combined",     14,  100L, 0.0105599999999998, 0.00184711460040343, 0.0579842820109458,
#   "Combined",     17,  200L, 0.0431150000000002,  0.0221223112761615, 0.0823508337277539,
#   "Combined",     18,  200L, 0.0434700000000001,   0.022299970205302, 0.0830306178613479,
#   "Combined",     19,  200L, 0.0798099999999993,  0.0493467575514476,  0.126575203056695,
#   "Combined",     20,  200L, 0.0689349999999998,  0.0407436413657359,  0.114308116426858,
#   "Combined",     21,  200L, 0.0750750000000008,  0.0457905956485781,  0.120718292288378,
#   "Combined",     22,  200L, 0.0355349999999998,  0.0172721597953066, 0.0716993107746621,
#   "Combined",     23,  200L, 0.0844649999999996,  0.0527466648833784,  0.132586944479438,
#   "Combined",     24,  200L, 0.0694549999999997,  0.0414740782661843,   0.11406705818783,
#   "Combined",     25,  200L,           0.063885,  0.0366526359995393,  0.109060011958329,
#   "Combined",     30,  200L,           0.129955,  0.0881953852243244,  0.187423042342957,
#   "Combined",     31,  200L, 0.0919100000000001,  0.0583328478485414,  0.141902101175647,
#   "Combined",     32,  200L,           0.113775,  0.0755409717516833,  0.167847126993757,
#   "Combined",     33,  200L,            0.10752,  0.0711645331155792,  0.159263610639403,
#   "Combined",     34,  200L,           0.128185,  0.0874676569420377,  0.184034176973193,
#   "Combined",     45,  200L,           0.116765,  0.0785313072613766,  0.170175435429028,
#   "Combined",     46,  200L,            0.14462,   0.101157728048079,  0.202547830993344,
#   "Combined",     47,  200L,            0.13018,  0.0891305166227496,  0.186268917111879,
#   "Combined",     48,  200L,            0.12857,  0.0874735682462468,  0.185058597470984,
#   "Combined",     49,  200L,  0.113939999999999,  0.0760462408645264,  0.167297201965587,
#   "Combined",     50,  200L,            0.15694,   0.111767424356249,  0.215931416896436
# )
# 

# dopico_pred <- tribble(
#   ~Week,  ~Bayesian_est,  ~SPIKE_3SD,  ~RBD_3SD,  ~SPIKE_6SD,  ~RBD_6SD,
#   14,            2.4,         2.5,       0.5,         0.5,       0.5,
#   17,            4.8,           6,         5,         4.5,         4,
#   18,            5.4,           8,         5,           5,         4,
#   19,            5.9,        11.5,         8,           9,         8,
#   20,            6.3,           8,         9,           7,         8,
#   21,            6.7,          12,        10,           8,       7.5,
#   22,              7,         4.5,       3.5,           4,       3.5,
#   23,            7.3,        10.5,         9,         9.5,         7,
#   24,            7.8,           8,         8,         6.5,         7,
#   25,            8.3,         8.5,         7,         7.5,         7,
#   30,           11.4,          22,      20.5,          17,      14.5,
#   31,           11.9,        14.5,      13.5,        10.5,         9,
#   32,           12.3,          18,        16,        13.5,      10.5,
#   33,           12.5,        15.5,      12.5,        12.5,        10,
#   34,           12.7,          20,        20,        16.5,       9.5,
#   45,           14.1,        15.5,        14,        12.5,        11,
#   46,           14.2,          21,        17,        16.5,      13.5,
#   47,           14.3,        17.5,        16,          17,        12,
#   48,           14.5,        20.5,        18,          16,      12.5,
#   49,           14.7,        16.5,      14.5,        12.5,        12,
#   50,           14.8,          20,      17.5,          18,      15.5
# )
# 
# 
# dopico_est <- dopico_pred %>%
#   mutate(
#     Model = "Dopico"
#   ) %>%
#   # group_by(Week, Model) %>%
#   mutate(
#     "Seroprevalence" = 1 * estimate,
#     "97.5%" = 1 * upper.ci,
#     "02.5%" = 1 * lower.ci,
#     Name = "Dopico",
#     Date = as.Date(paste0("2020", Week, "7"), "%Y%U%u"),
#     m_scale = NA,
#     rho = NA,
#     theta = NA
#   ) %>%
#   select(-c(type, estimate, lower.ci, upper.ci))

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


# output_df$Week <- output_df$Time
output_df$Date <- as.Date(paste0("2020", output_df$Week, "7"), "%Y%U%u")

output_df_new <- output_df %>% 
  mutate(`97.5%` = 0.0,
         `02.5%` = 0.0,
         Name = Model,
         m_scale = NA,
         rho = NA, 
         theta = NA,
         Seroprevalence = Seroprevalence / 100 )


sero_all <- svm_lda_est %>%
  rbind(distinct(my_sero_df), output_df_new) %>% 
  select(-N) %>% 
  rbind(bayesian_learner_est)

p_sero_all <- sero_all %>% 
  ggplot(aes(x = Date)) +
  geom_point(aes(
    y = 100 * Seroprevalence,
    colour = factor(Model), # , labels = c("Dopico et al.", "MVT")),
    group = Name
  ))  +
  labs(
    y = "Seroprevalence (%)",
    colour = "Model:"
  ) +
  theme(legend.position = "bottom") +
  geom_line(aes(
    y = 100 * Seroprevalence,
    colour = factor(Model), # , labels = c("Dopico et al.", "MVT")),
    group = Name
  )) +
  ggthemes::scale_color_colorblind()


dopico_est %>%
  rbind(my_sero_df, output_df_new) %>% 
  group_by(Name) %>% 
  summarise(S_n = sum(Seroprevalence))


dopico_est %>%
  rbind(my_sero_df, output_df_new) %>%
  ggplot(aes(x = Name)) +
  geom_boxplot(aes(
    y = 100 * Seroprevalence,
    colour = factor(Model), # , labels = c("Dopico et al.", "MVT")),
    group = Name
  ))  +
labs(
  y = "Seroprevalence (%)",
  colour = "Model:"
) +
  theme(legend.position = "bottom") +
  geom_line(aes(
    y = 100 * Seroprevalence,
    colour = factor(Model), # , labels = c("Dopico et al.", "MVT")),
    group = Name
  )) +
  ggthemes::scale_color_colorblind()


dopico_est %>%
  rbind(distinct(my_sero_df), output_df_new) %>% 
  mutate(Total_estimated = Seroprevalence * N) %>% 
  group_by(Name) %>% 
  summarise(N_total = sum(Total_estimated))


ggsave(paste0(my_wd, "Analysis/Sweden/Outputs/Plots/seroprevalenceAllModels.png"),
  plot = p_sero_all,
  width = 6,
  height = 5
)

rep_model_df <- my_sero_df %>% 
  filter(Name == "MVT_0.1_11_5") %>% 
  mutate(Name = "MVT",
         Date = as.Date(paste0("2020", Week, "7"), "%Y%U%u")
         )


bayesian_learner_est_reduced <- bayesian_learner_est %>%
  filter(Week %in% weeks_present)

all_models_sero_df <- svm_lda_est %>%
  rbind(distinct(rep_model_df), output_df_new) %>% 
  select(-N) %>% 
  rbind(bayesian_learner_est_reduced)
  

comparison_sero_table <- all_models_sero_df %>% 
  select(Date, Seroprevalence, Name) %>% 
  mutate(Seroprevalence = Seroprevalence * 100) %>% 
  pivot_wider(values_from = Seroprevalence, names_from = Name) 

## Standard xtable use
xtable::xtable(comparison_sero_table)

# \begin{table}[ht]
#  \centering
#  \begin{tabular}{c|cccccc}
#  \hline
#    Date & SVM-LDA  & Bayesian learner & MVT & RF & SVM & LR & LR (BC) \\
#  \hline
#    2020/04/05 &  NA & 2.35 & 1.60 & 1.00 & 1.00 & 1.00 & 2.00  \\ 
#    2020/04/26 & 4.26 & 4.73 & 4.92 & 4.50 & 5.00 & 5.00 & 5.00  \\ 
#    2020/05/03 & 4.33 & 5.34 & 5.74 & 4.50 & 4.50 & 5.00 & 5.50  \\ 
#    2020/05/10 & 7.87 & 5.86 & 8.87 & 8.50 & 8.50 & 8.50 & 10.50  \\ 
#    2020/05/17 & 7.36 & 6.28 & 8.05 & 7.50 & 4.50 & 8.00 & 8.50  \\ 
#    2020/05/24 & 7.49 & 6.64 & 8.34 & 8.00 & 6.50 & 7.50 & 9.50  \\ 
#    2020/05/31 & 3.54 & 6.97 & 4.15 & 4.00 & 3.50 & 4.00 & 5.00  \\ 
#    2020/06/07 & 8.41 & 7.31 & 9.83 & 9.50 & 8.00 & 10.00 & 10.50  \\ 
#    2020/06/14 & 6.90 & 7.75 & 7.55 & 7.00 & 7.00 & 7.00 & 8.50  \\ 
#    2020/06/21 & 6.44 & 8.30 & 7.80 & 7.00 & 6.50 & 7.50 & 8.00  \\ 
#    2020/07/26 & 13.20 & 11.34 & 16.72 & 15.00 & 11.50 & 16.50 & 16.00  \\ 
#    2020/08/02 & 9.16 & 11.79 & 10.48 & 9.00 & 8.50 & 10.00 & 10.00  \\ 
#    2020/08/09 & 11.30 & 12.15 & 13.43 & 12.00 & 8.00 & 13.00 & 13.50  \\ 
#    2020/08/16 & 10.84 & 12.43 & 12.15 & 12.00 & 10.50 & 11.50 & 11.00  \\ 
#    2020/08/23 & 12.97 & 12.65 & 15.55 & 15.50 & 11.00 & 15.00 & 15.00  \\ 
#    2020/11/08 & 11.79 & 14.28 & 13.72 & 12.00 & 11.50 & 12.50 & 14.00  \\ 
#    2020/11/15 & 14.20 & 14.47 & 18.85 & 15.50 & 14.50 & 17.00 & 18.50  \\ 
#    2020/11/22 & 13.37 & 14.72 & 16.44 & 15.50 & 15.50 & 15.50 & 16.00  \\ 
#    2020/11/29 & 13.20 & 14.98 & 17.36 & 15.00 & 15.00 & 16.00 & 16.50  \\ 
#    2020/12/06 & 11.52 & 15.29 & 14.47 & 12.50 & 11.00 & 13.00 & 13.50  \\ 
#    2020/12/13 & 15.73 & 15.64 & 19.65 & 18.00 & 16.00 & 18.00 & 19.00  \\
#  \hline
#  \end{tabular}
#  \end{table}


