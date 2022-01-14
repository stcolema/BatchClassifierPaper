# /usr/bin/Rscript
#
# Analysis of the Seattle hosptial data from
# Dingens, A.S., Crawford, K.H.D., Adler, A. et al. (2020).
#
# === Setup ====================================================================
library(BatchMixtureModel)
library(mdiHelpR)
library(ggplot2)
library(magrittr)
library(tibble)
library(dplyr)
library(tidyr)
library(data.table)


setMyTheme()
set.seed(1)

my_wd <- "/home/sdc56/rds/hpc-work/BatchMixtureModel"

rbd_screen_results <- read.csv(paste0(my_wd, "/Data/Seattle/Original/RBD_screen_results.csv"))
save_dir <- paste0(my_wd, "/Analysis/Seattle/Plots/")

dropped_samples <- which(rbd_screen_results$category == "no sera")
rbd_screen_results <- rbd_screen_results[-dropped_samples, ]


cutoffs <- list(
  "3SD" = unique(rbd_screen_results$Cutoff_3SD),
  "5SD" = unique(rbd_screen_results$Cutoff_5SD)
)

X <- log(matrix(rbd_screen_results$Screen_OD450, ncol = 1))
colnames(X) <- "RBD"
N <- nrow(X)
hist(X)

batch_vec <- rbd_screen_results$screen_batch
types <- factor(rbd_screen_results$category)

initial_labels <- rep(0, N)

neg_controls <- which(types == "pre-2020 sera pool")
pos_controls <- which(types == "CR3022 antibody")
control_inds <- c(neg_controls, pos_controls)

N_controls <- length(control_inds)

hist(X[control_inds])
hist(X[neg_controls])
hist(X[pos_controls])

initial_labels[neg_controls] <- 0
initial_labels[pos_controls] <- 1
initial_labels[-control_inds] <- sample(c(0, 1), size = N - N_controls, replace = T, prob = c(0.5, 0.5))

fixed <- rep(0, N)
fixed[control_inds] <- 1

observed_labels <- rep("Unknown", N)
observed_labels[pos_controls] <- "Seropositive"
observed_labels[neg_controls] <- "Seronegative"

# Let's plot the original dataset
plot_df <- X %>%
  as_tibble() %>%
  add_column(
    Label = factor(observed_labels),
    Batch = factor(batch_vec)
  )

p_raw_data <- plot_df %>%
  # mutate(Plot_label = as.numeric(factor(Label)) * fixed) %>%
  ggplot(aes(x = Label, y = RBD, colour = Label)) +
  geom_jitter(size = 0.7) +
  labs(
    title = "Seattle ELISA data",
    colour = "Group"
  ) +
  facet_wrap(~Batch, labeller = label_both) +
  ggthemes::scale_color_colorblind() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1))

# ggsave("./Data/ELISA_observed.png", plot = p_raw_data, width = 5, height = 3.75)

ggsave(paste0(save_dir, "ELISA_observed.png"),
  plot = p_raw_data,
  width = 8,
  height = 5
)



p_controls <- plot_df[fixed == 1, ] %>%
  # mutate(Plot_label = as.numeric(factor(Label, levels = c("Historical controls", "COVID")))) %>%
  ggplot(aes(
    x = Label,
    y = RBD,
    colour = Label
  )) +
  geom_jitter(size = 0.7) +
  labs(
    title = "ELISA data",
    subtitle = "Controls",
    colour = "Class"
  ) +
  facet_wrap(~Batch, labeller = label_both) +
  ggthemes::scale_color_colorblind() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1))

ggsave(paste0(save_dir, "ELISA_controls.png"),
  plot = p_controls,
  width = 8,
  height = 5
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
R <- 15000
thin <- 25
burn <- 5000
eff_burn <- burn / thin

# Number of chains to run
n_chains <- 7L
mvt_samples <- mvn_samples <- vector("list", n_chains * n_hyperparam_combinations)

# What chain are we currently saving
chain_count <- 0


for (j in seq(1, n_hyperparam_combinations)) {

  # Hyperparameters for the prior distribution of the batch effects
  .m_scale <- hyper_params$m_scale[j]
  .rho <- hyper_params$rho[j]
  .theta <- hyper_params$theta[j]

  # Generate a Markov chain for each model and record acceptance rates separately
  for (i in 1:n_chains) {
    set.seed(i)

    # Modelling
    t0 <- Sys.time()
    mvt_samples[[i]] <- .mvt <- batchSemiSupervisedMixtureModel(
      X,
      R,
      thin,
      initial_labels,
      fixed,
      batch_vec,
      "MVT",
      alpha = 1,
      mu_proposal_window = 0.16**2,
      cov_proposal_window = 100,
      m_proposal_window = 0.17**2,
      S_proposal_window = 50,
      t_df_proposal_window = 20,
      m_scale = .m_scale,
      rho = .rho,
      theta = .theta
    )

    t1 <- Sys.time()

    # Modelling MVN
    set.seed(i)
    t2 <- Sys.time()
    mvn_samples[[i]] <- .mvn <- batchSemiSupervisedMixtureModel(
      X,
      R,
      thin,
      initial_labels,
      fixed,
      batch_vec,
      "MVN",
      alpha = 1,
      mu_proposal_window = 0.16**2,
      cov_proposal_window = 125,
      m_proposal_window = 0.17**2,
      S_proposal_window = 50,
      m_scale = .m_scale,
      rho = .rho,
      theta = .theta
    )

    t3 <- Sys.time()

    # Model likelihood
    likelihood_df_local <- data.frame(
      "mvt_complete" = .mvt$complete_likelihood,
      "mvt_observed" = .mvt$observed_likelihood,
      "mvn_complete" = .mvn$complete_likelihood,
      "mvn_observed" = .mvn$observed_likelihood,
      "mvt_bic" = .mvt$BIC,
      "mvn_bic" = .mvn$BIC,
      "Chain" = i,
      "Iteration" = (1:(R / thin)) * thin,
      "m_scale" = .m_scale,
      "rho" = .rho,
      "theta" = .theta
    ) %>%
      pivot_longer(c(
        mvt_complete,
        mvt_observed,
        mvn_complete,
        mvn_observed,
        mvt_bic,
        mvn_bic
      ),
      names_sep = "_",
      names_to = c("model", "param")
      )

    # Acceptance rates
    mvt_cluster_acceptance_df <- data.frame(
      Mu = .mvt$mu_acceptance_rate,
      Cov = .mvt$cov_acceptance_rate,
      DF = .mvt$t_df_acceptance_rate,
      K = seq(1, length(.mvt$mu_acceptance_rate)),
      Model = "MVT",
      Chain = i,
      m_scale = .m_scale,
      rho = .rho,
      theta = .theta
    )

    mvt_batch_acceptance_df <- data.frame(
      m = .mvt$m_acceptance_rate,
      S = .mvt$S_acceptance_rate,
      B = seq(1, length(.mvt$m_acceptance_rate)),
      Model = "MVT",
      Chain = i,
      m_scale = .m_scale,
      rho = .rho,
      theta = .theta
    )

    mvn_cluster_acceptance_df <- data.frame(
      Mu = .mvn$mu_acceptance_rate,
      Cov = .mvn$cov_acceptance_rate,
      DF = NA,
      K = seq(1, length(.mvn$mu_acceptance_rate)),
      Model = "MVN",
      Chain = i,
      m_scale = .m_scale,
      rho = .rho,
      theta = .theta
    )

    mvn_batch_acceptance_df <- data.frame(
      m = .mvn$m_acceptance_rate,
      S = .mvn$S_acceptance_rate,
      B = seq(1, length(.mvn$m_acceptance_rate)),
      Model = "MVN",
      Chain = i,
      m_scale = .m_scale,
      rho = .rho,
      theta = .theta
    )


    if (i == 1 & j == 1) {
      batch_acceptance_df <- rbind(mvn_batch_acceptance_df, mvt_batch_acceptance_df)
      cluster_acceptance_df <- rbind(mvn_cluster_acceptance_df, mvt_cluster_acceptance_df)
      likelihood_df <- likelihood_df_local
    } else {
      cluster_acceptance_df <- rbind(
        cluster_acceptance_df,
        mvn_cluster_acceptance_df,
        mvt_cluster_acceptance_df
      )
      batch_acceptance_df <- rbind(
        batch_acceptance_df,
        mvn_batch_acceptance_df,
        mvt_batch_acceptance_df
      )

      likelihood_df <- rbind(likelihood_df, likelihood_df_local)
    }

    # For saving numbers with a decimal place
    m_scale_sci_notation <- formatC(.m_scale, format = "e", digits = 0)

    # Save the samples for the current chain
    saveRDS(.mvt, file = paste0(
      my_wd,
      "/Analysis/Seattle/Outputs/seattle_mvt_chain_",
      i,
      "_m_scale_",
      m_scale_sci_notation,
      "_rho_",
      .rho,
      "_theta_",
      .theta,
      ".rds"
    ))

    saveRDS(.mvn, file = paste0(
      my_wd,
      "/Analysis/Seattle/Outputs/seattle_mvn_chain_",
      i,
      "_m_scale_",
      m_scale_sci_notation,
      "_rho_",
      .rho,
      "_theta_",
      .theta,
      ".rds"
    ))
  }
}



# Save the various data frames that we use to check model behaviour
write.csv(cluster_acceptance_df,
  paste0(my_wd, "/Analysis/Seattle/ModelChecks/groupParamAcceptance.csv"),
  row.names = FALSE
)

write.csv(batch_acceptance_df,
  paste0(my_wd, "/Analysis/Seattle/ModelChecks/batchParamAcceptance.csv"),
  row.names = FALSE
)

write.csv(likelihood_df,
  paste0(my_wd, "/Analysis/Seattle/ModelChecks/likelihoods.csv"),
  row.names = FALSE
)

# # Check that the acceptance rates are reasonable (ideally in [0.1, 0.5] but
# # avoid really high or low values).
# batch_acceptance_df %>%
#   pivot_longer(c(S, m), names_to = "Parameter", values_to = "Acceptance_rate") %>%
#   ggplot(aes(x = Chain, y = Acceptance_rate, colour = Parameter)) +
#   geom_point() +
#   facet_grid(Parameter ~ Model)
#
# cluster_acceptance_df %>%
#   pivot_longer(c(Mu, Cov, DF), names_to = "Parameter", values_to = "Acceptance_rate") %>%
#   ggplot(aes(x = Chain, y = Acceptance_rate, colour = Parameter)) +
#   geom_point() +
#   facet_grid(Parameter ~ Model, scales = "free_y") +
#   geom_hline(yintercept = c(0.1, 0.5))
#
# # Check out the model likelihoods
# likelihood_df %>%
#   filter(Iteration > burn) %>%
#   ggplot(aes(x = Iteration, y = value, colour = model, group = interaction(model, Chain))) +
#   geom_point() +
#   facet_wrap(~param, ncol = 1, scales = "free_y") +
#   geom_line() +
#   labs(title = "Complete log-likelihood", y = "Log-likelihood") +
#   ggthemes::scale_color_colorblind()
#
# ggsave(paste0(save_dir, "Seattle_model_fit_mvn_mvt.png"),
#   height = 6,
#   width = 4
# )
#
#
# likelihood_df %>%
#   filter(Iteration > burn, param == "complete") %>%
#   ggplot(aes(x = Iteration, y = value, colour = factor(Chain), group = interaction(model, Chain))) +
#   geom_point() +
#   facet_wrap(~model, ncol = 1, scales = "free_y") +
#   geom_line() +
#   labs(title = "Complete log-likelihood", y = "Log-likelihood") +
#   ggthemes::scale_color_colorblind()
#
#
# ggsave(paste0(save_dir, "Seattle_model_complete_ll.png"),
#   height = 5,
#   width = 4
# )
#
#
# chain_used <- 7
#
# likelihood_df %>%
#   filter(Iteration > burn, param == "complete", Chain == chain_used) %>%
#   ggplot(aes(x = Iteration, y = value, colour = model, group = interaction(model, Chain))) +
#   geom_point() +
#   # facet_wrap(~model, ncol = 1, scales = "free_y") +
#   geom_line() +
#   labs(title = "Complete log-likelihood", y = "Log-likelihood") +
#   ggthemes::scale_color_colorblind()
#
# saveRDS(mvt_samples[[chain_used]],
#   file = paste0("./Analysis/Seattle/seattle_mvt_chain_", chain_used, ".rds")
# )
# saveRDS(mvn_samples[[chain_used]],
#   file = paste0("./Analysis/Seattle/seattle_mvn_chain_", chain_used, ".rds")
# )
#
#
# # === Sampled parameters =======================================================
#
#
# plotSampledBatchScales(mvt_samples[[1]], R = R, thin = thin, burn_in = burn)
# plotSampledBatchScales(mvt_samples[[2]], R = R, thin = thin, burn_in = burn)
# plotSampledBatchScales(mvt_samples[[3]], R = R, thin = thin, burn_in = burn)
# plotSampledBatchScales(mvt_samples[[4]], R = R, thin = thin, burn_in = burn)
# plotSampledBatchScales(mvt_samples[[5]], R = R, thin = thin, burn_in = burn)
# plotSampledBatchScales(mvt_samples[[6]], R = R, thin = thin, burn_in = burn)
#
# plotSampledBatchMeans(mvt_samples[[1]], R = R, thin = thin, burn_in = burn)
# plotSampledBatchMeans(mvt_samples[[2]], R = R, thin = thin, burn_in = burn)
# plotSampledBatchMeans(mvt_samples[[3]], R = R, thin = thin, burn_in = burn)
# plotSampledBatchMeans(mvt_samples[[4]], R = R, thin = thin, burn_in = burn)
# plotSampledBatchMeans(mvt_samples[[5]], R = R, thin = thin, burn_in = burn)
# plotSampledBatchMeans(mvt_samples[[6]], R = R, thin = thin, burn_in = burn)
# plotSampledBatchMeans(mvt_samples[[7]], R = R, thin = thin, burn_in = burn)
#
# plotSampledClusterMeans(mvt_samples[[1]], R = R, thin = thin, burn_in = burn)
# plotSampledClusterMeans(mvt_samples[[2]], R = R, thin = thin, burn_in = burn)
# plotSampledClusterMeans(mvt_samples[[3]], R = R, thin = thin, burn_in = burn)
# plotSampledClusterMeans(mvt_samples[[4]], R = R, thin = thin, burn_in = burn)
# plotSampledClusterMeans(mvt_samples[[5]], R = R, thin = thin, burn_in = burn)
# plotSampledClusterMeans(mvt_samples[[6]], R = R, thin = thin, burn_in = burn)
# plotSampledClusterMeans(mvt_samples[[7]], R = R, thin = thin, burn_in = burn)
#
# plotSampledBatchScales(mvn_samples[[1]], R = R, thin = thin, burn_in = burn)
# plotSampledBatchScales(mvn_samples[[2]], R = R, thin = thin, burn_in = burn)
# plotSampledBatchScales(mvn_samples[[3]], R = R, thin = thin, burn_in = burn)
# plotSampledBatchScales(mvn_samples[[4]], R = R, thin = thin, burn_in = burn)
# plotSampledBatchScales(mvn_samples[[5]], R = R, thin = thin, burn_in = burn)
# plotSampledBatchScales(mvn_samples[[6]], R = R, thin = thin, burn_in = burn)
#
# plotSampledBatchMeans(mvn_samples[[1]], R = R, thin = thin, burn_in = burn)
# plotSampledBatchMeans(mvn_samples[[2]], R = R, thin = thin, burn_in = burn)
# plotSampledBatchMeans(mvn_samples[[3]], R = R, thin = thin, burn_in = burn)
# plotSampledBatchMeans(mvn_samples[[4]], R = R, thin = thin, burn_in = burn)
# plotSampledBatchMeans(mvn_samples[[5]], R = R, thin = thin, burn_in = burn)
# plotSampledBatchMeans(mvn_samples[[6]], R = R, thin = thin, burn_in = burn)
# plotSampledBatchMeans(mvn_samples[[7]], R = R, thin = thin, burn_in = burn)
#
# plotSampledClusterMeans(mvn_samples[[1]], R = R, thin = thin, burn_in = burn)
# plotSampledClusterMeans(mvn_samples[[2]], R = R, thin = thin, burn_in = burn)
# plotSampledClusterMeans(mvn_samples[[3]], R = R, thin = thin, burn_in = burn)
# plotSampledClusterMeans(mvn_samples[[4]], R = R, thin = thin, burn_in = burn)
# plotSampledClusterMeans(mvn_samples[[5]], R = R, thin = thin, burn_in = burn)
# plotSampledClusterMeans(mvn_samples[[6]], R = R, thin = thin, burn_in = burn)
# plotSampledClusterMeans(mvn_samples[[7]], R = R, thin = thin, burn_in = burn)
#
# # === Allocation ===============================================================
#
# # Allocations
# mvn_prob <- calcAllocProb(mvn_samples[[chain_used]]$alloc, eff_burn)
# mvn_pred <- predictClass(mvn_prob)
#
# mvt_prob <- calcAllocProb(mvt_samples[[chain_used]]$alloc, eff_burn)
# mvt_pred <- predictClass(mvt_prob)
#
# # Inferred datasets
# mvt_inferred_data <- rowMeans(mvt_samples[[chain_used]]$batch_corrected_data[, , -c(1:eff_burn), drop = FALSE], dims = 2) %>%
#   as_tibble(.name_repair = "minimal") %>%
#   set_colnames(colnames(X)) %>%
#   add_column(
#     "Label" = factor(mvt_pred),
#     "Prob" = apply(mvt_prob, 1, function(x) {
#       x[which.max(x)]
#     }),
#     "Batch" = factor(batch_vec),
#     "Type" = "Inferred",
#     "Fixed" = factor(fixed),
#     "Model" = "MVT"
#   )
#
# p_mvt <- mvt_inferred_data %>%
#   ggplot(aes(y = RBD, x = Label)) +
#   # geom_point(aes(shape = Fixed, colour = Label, alpha = Prob)) +
#   geom_jitter(aes(shape = Fixed, colour = Label, alpha = Prob)) +
#   scale_alpha_continuous(range = c(0.4, 1.0)) +
#   # geom_density_2d(aes(colour = Predicted)) +
#   scale_color_viridis_d() +
#   labs(
#     title = "Seattle ELISA data - MVT",
#     subtitle = "Predicted labels",
#     y = "Inferred RBD",
#     caption = "log transformed"
#   )
#
# mvn_inferred_data <- rowMeans(mvn_samples[[chain_used]]$batch_corrected_data[, , -c(1:eff_burn), drop = F], dims = 2) %>%
#   as_tibble(.name_repair = "minimal") %>%
#   set_colnames(colnames(X)) %>%
#   add_column(
#     "Label" = factor(mvn_pred),
#     "Prob" = apply(mvn_prob, 1, function(x) {
#       x[which.max(x)]
#     }),
#     "Batch" = factor(batch_vec),
#     "Type" = "Inferred",
#     "Fixed" = factor(fixed),
#     "Model" = "MVN"
#   )
#
# p_mvn <- mvn_inferred_data %>%
#   ggplot(aes(y = RBD, x = Label)) +
#   # geom_point(aes(shape = Fixed, colour = Label, alpha = Prob)) +
#   geom_jitter(aes(shape = Fixed, colour = Label, alpha = Prob)) +
#   scale_alpha_continuous(range = c(0.4, 1.0)) +
#   # geom_density_2d(aes(colour = Predicted)) +
#   scale_color_viridis_d() +
#   labs(
#     title = "Seattle ELISA data - MVN",
#     subtitle = "Predicted labels",
#     y = "Inferred RBD",
#     caption = "log transformed"
#   )
#
# observed_labels <- initial_labels + 1
# observed_labels[fixed == 0] <- 3
#
# observed_data <- X %>%
#   as_tibble() %>%
#   add_column(
#     "Label" = observed_labels,
#     "Prob" = 1,
#     "Batch" = factor(batch_vec),
#     "Type" = "Observed",
#     "Fixed" = factor(fixed),
#     "Model" = "Observed"
#   )
#
# p_all_datasets <- rbind(observed_data, mvt_inferred_data, mvn_inferred_data) %>%
#   mutate(Label = factor(Label, levels = c(1, 2, 3), labels = c("Seronegative", "Seropositive", "Unknown"))) %>%
#   ggplot(aes(y = RBD, x = Batch)) +
#   geom_jitter(aes(shape = Fixed, colour = Label, alpha = Prob)) +
#   scale_alpha_continuous(range = c(0.4, 1.0)) +
#   facet_wrap(~Model, ncol = 1) +
#   # geom_density_2d(aes(colour = Predicted)) +
#   ggthemes::scale_color_colorblind() +
#   labs(
#     title = "Elisa data",
#     subtitle = "Observed and inferred datasets and labels",
#     caption = "log transformed",
#     colour = "Group",
#     alpha = "Probability"
#   )
#
# # rbind(observed_data, mvt_inferred_data, mvn_inferred_data) %>%
# #   mutate(Label = factor(Label, levels = c(1, 2, 3), labels = c("Seronegative", "Seropositive", "Unknown"))) %>%
# #   ggplot(aes(x = SPIKE, y = RBD)) +
# #   geom_point(aes(shape = Fixed, colour = Label, alpha = Prob), size = 0.4) +
# #   scale_alpha_continuous(range = c(0.4, 1.0)) +
# #   facet_wrap(~Model, ncol = 1) +
# #   # geom_density_2d(aes(colour = Predicted)) +
# #   ggthemes::scale_color_colorblind() +
# #   labs(
# #     title = "Elisa data",
# #     subtitle = "Observed and inferred datasets and labels",
# #     caption = "log transformed",
# #     colour = "Class",
# #     alpha = "Probability"
# #   )
#
# ggsave(paste0(save_dir, "Seattle_full_chain_", chain_used, "_mvn_v_mvt_inferred_norm_prob.png"),
#   plot = p_all_datasets,
#   height = 7,
#   width = 6
# )
#
#
# # #   # Modelling
# # # t0 <- Sys.time()
# # # mvt_samples <- batchSemiSupervisedMixtureModel(
# # #   X,
# # #   R,
# # #   thin,
# # #   initial_labels,
# # #   fixed,
# # #   batch_vec,
# # #   "MVT",
# # #   alpha = 1,
# # #   mu_proposal_window = 0.10**2,
# # #   cov_proposal_window = 200,
# # #   m_proposal_window = 0.10**2,
# # #   S_proposal_window = 75,
# # #   t_df_proposal_window = 25
# # # )
# # #
# # # t1 <- Sys.time()
# # #
# # # t1 - t0
# #
# # mvt_samples$cov_acceptance_rate
# # mvt_samples$mu_acceptance_rate
# # mvt_samples$t_df_acceptance_rate
# # mvt_samples$m_acceptance_rate
# # mvt_samples$S_acceptance_rate
# #
# # plot(mvt_samples$t_df[, 1])
# # plot(mvt_samples$t_df[, 2])
# # plot(mvt_samples$complete_likelihood)
# # plot(mvt_samples$complete_likelihood[-c(1:100)])
# #
# # plotSampledBatchScales(mvt_samples, R = R, thin = thin, burn_in = 0)
# # plotSampledBatchMeans(mvt_samples, R = R, thin = thin, burn_in = 0)
# # plotSampledClusterMeans(mvt_samples, R = R, thin = thin, burn_in = 0)
# #
# # plotSampledBatchScales(mvt_samples, R = R, thin = thin, burn_in = burn)
# # plotSampledBatchMeans(mvt_samples, R = R, thin = thin, burn_in = burn)
# # plotSampledClusterMeans(mvt_samples, R = R, thin = thin, burn_in = burn)
# #
# #
# # # Allocations
# # mvt_prob <- calcAllocProb(mvt_samples$alloc, eff_burn)
# # mvt_pred <- predictClass(mvt_prob)
# #
# # # Inferred datasets
# # mvt_inferred_data <- rowMeans(mvt_samples$batch_corrected_data[, , -c(1:eff_burn), drop = FALSE], dims = 2) %>%
# #   as_tibble(.name_repair = "minimal") %>%
# #   set_colnames(colnames(X)) %>%
# #   add_column(
# #     "Label" = factor(mvt_pred),
# #     "Prob" = apply(mvt_prob, 1, function(x) {
# #       x[which.max(x)]
# #     }),
# #     "Batch" = factor(batch_vec),
# #     "Type" = "Inferred",
# #     "Fixed" = factor(fixed),
# #     "Model" = "MVT"
# #   )
# #
# # mvt_inferred_data %>%
# #   ggplot(aes(y = RBD, x = Label)) +
# #   # geom_point(aes(shape = Fixed, colour = Label, alpha = Prob)) +
# #   geom_jitter(aes(shape = Fixed, colour = Label, alpha = Prob)) +
# #   scale_alpha_continuous(range = c(0.4, 1.0)) +
# #   # geom_density_2d(aes(colour = Predicted)) +
# #   scale_color_viridis_d() +
# #   labs(
# #     title = "Elisa data - batch adjusted",
# #     subtitle = "Predicted labels",
# #     caption = "log transformed"
# #   )
#
#
# # === Allocations across time ==================================================
#
# non_control_data <- rbd_screen_results[-control_inds, ]
#
# # The non-control data has the week of the sample collection included in the
# # sample ID
# weeks_str <- non_control_data$screen_batch
# library(stringr)
# # Extract the week
# periods <- weeks_str %>%
#   str_extract_all("(\\d{1})") %>%
#   do.call(rbind, .) %>%
#   as.numeric()
#
# periods_present <- unique(periods)
#
# # Associate the allocation probabilities with the week
# week_prob_df <- data.frame(
#   Week = periods,
#   Seronegative_prob = mvt_prob[-control_inds, 1],
#   Seropositive_prob = mvt_prob[-control_inds, 2],
#   Allocation = mvt_pred[-control_inds]
# )
#
# # Plot the number of seropositive allocations as a function of time
# p_alloc_time <- week_prob_df %>%
#   group_by(Week) %>%
#   summarise(
#     Seronegative = sum(Allocation == 1),
#     Seropositive = sum(Allocation == 2)
#   ) %>%
#   ggplot(aes(x = Week, y = Seropositive)) +
#   geom_point() +
#   labs(
#     title = "Seropositive allocations across time",
#     y = "Number of serorpositive allocations"
#   )
#
# ggsave(paste0(save_dir, "Positive_allocations_across_time.png"),
#   height = 5,
#   width = 5,
#   plot = p_alloc_time
# )
#
# # === Other stuff ==============================================================
#
# mvt_pred %>% table()
# mvn_pred %>% table()
#
#
# papers_pred <- which(rbd_screen_results$Hit_5SD == "True")
# mvt_pred[papers_pred]
# mvn_pred[papers_pred]
#
# plot_df[papers_pred, ] %>%
#   # mutate(Plot_label = as.numeric(factor(Label)) * fixed) %>%
#   ggplot(aes(x = Label, y = RBD, colour = Label)) +
#   geom_jitter(size = 0.7) +
#   labs(
#     title = "Seattle ELISA data",
#     colour = "Group"
#   ) +
#   facet_wrap(~Batch, labeller = label_both) +
#   ggthemes::scale_color_colorblind() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1))
#
# serology_results <- read.csv("./Data/Seattle/Original/serology_results.csv")
#
# rbd_order <- order(rbd_screen_results[papers_pred, ]$Screen_OD450)
# reduced_df <- rbd_screen_results[papers_pred, ]
#
# reduced_further <- reduced_df[reduced_df$Screen_OD450 %in% serology_results[serology_results$positive_at_5_sd == "True", ]$RBD.IgG.screen, ]
#
# preds_of_interest <- as.numeric(row.names(reduced_further[order(reduced_further$Screen_OD450), ]))
#
# mvt_pred[preds_of_interest]
# mvn_pred[preds_of_interest]
#
# serology_results[serology_results$positive_at_5_sd == "True", ]
#
#
# serology_results[serology_results$Abbott.Interp != "", ]
#
# dim(rbd_screen_results)
#
# rbd_kept <- which(rbd_screen_results$category == "population samples" &
#   !rbd_screen_results$sample %in% c("COV_0001", "COV_0004"))
#
# rbd_pop <- rbd_screen_results[rbd_kept, ]
# mvt_pred_pop <- mvt_pred[rbd_kept]
# mvn_pred_pop <- mvn_pred[rbd_kept]
#
# serology_kept <- serology_results[-c(1:2), ]
#
# serology_matching_order <- serology_kept[match(serology_kept$RBD.IgG.screen, rbd_pop$Screen_OD450), ]
#
#
#
# map_to_serology <- match(rbd_screen_results$sample, serology_kept$sample_id)
# map_to_serology <- map_to_serology[!is.na(map_to_serology)]
#
# serology_kept$MVT_pred <- as.character(mvt_pred[map_to_serology] == 2)
# serology_kept$MVN_pred <- as.character(mvn_pred[map_to_serology] == 2)
#
# serology_kept[serology_kept == "True"] <- TRUE
# serology_kept[serology_kept == "False"] <- FALSE
# serology_kept[serology_kept == ""] <- NA
#
# serology_kept$Abbott.Interp[serology_kept$Abbott.Interp == "Positive"] <- "TRUE"
# serology_kept$Abbott.Interp[serology_kept$Abbott.Interp == "Negative"] <- "FALSE"
#
# serology_kept %>%
#   select(-c(Abbott.index, n_assays_tested, sample_id, n_assay_hits_3sd, n_assay_hits_5sd)) %>%
#   pivot_longer(c(
#     RBD.IgG.screen.hit.3.sd,
#     RBD.IgG.screen.hit.5.sd,
#     RBD.IgG.titration.hit.3.sd,
#     Spike.IgG.titration.hit.3.sd,
#     RBD.IgG.titration.hit.5.sd,
#     Spike.IgG.titration.hit.5.sd,
#     Abbott.Interp,
#     Abbott.index.hit.5.sd,
#     Abbott.index.hit.3.sd,
#     positive_at_3_sd,
#     positive_at_5_sd,
#     MVT_pred,
#     MVN_pred
#   )) %>%
#   ggplot(aes(x = RBD.IgG.titration, y = Spike.IgG.titration, colour = value)) +
#   geom_point() +
#   facet_wrap(~name)
