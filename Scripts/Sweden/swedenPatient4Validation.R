#
# Validate our analysis of the ELISA data from Dopico et al., 2021 by using the
# ``Patient 4'' data. This is an extreme, seropositive individual that was 
# included multiple times partially to control for batch effects. We would hope 
# that after batch-correction the samples from this person move together.
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
library(ggforce)
library(patchwork)

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
batchPackageIsCurrent <- FALSE

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
patient_4_data <- m[patient_4, ]
m <- m[-patient_4, ]


# Find the controls
negative_controls <- which(m$type == "Historical controls")
positive_controls <- which((m$type == "COVID") | (m$type == "Patient 4"))
controls <- c(negative_controls, positive_controls)

N <- nrow(m)
N_positive <- length(positive_controls)
N_negative <- length(negative_controls)
N_controls <- length(controls)

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

control_df <- plot_df %>%
  filter(Fixed == 1)

control_df$Controls <- "Controls"
plot_df$Controls <- "Full dataset"

p_data <- rbind(plot_df, control_df) %>%
  ggplot(aes(x = SPIKE, y = RBD, colour = factor(Label,
    labels = c("Seronegative", "Seropositive", "Unknown")
  ))) +
  geom_point(size = 0.7) +
  facet_grid(~Controls) +
  labs(
    # title = "ELISA data",
    x = "Log-transformed observed SPIKE OD",
    y = "Log-transformed observed RBD OD",
    colour = "Class"
  ) +
  ggthemes::scale_color_colorblind()

# Model inputs
K_max <- 2
B <- length(unique(batch_vec))
P <- ncol(X)

# The class proportions in the obsered labels (no longer used)
# class_weights <- table(fixed) / N

# Our initial labelling is 1 for known controls, 2 for known COVIDs
initial_labels <- rep(1, N)
initial_labels[negative_controls] <- 1
initial_labels[positive_controls] <- 2

# We then sample with equal probability the class for the remaining items
initial_labels[-c(negative_controls, positive_controls)] <- sample(1:2,
  size = N - N_positive - N_negative,
  replace = TRUE,
  prob = c(1, 1)
)

# Convert the batch vector to integers
batch_vec <- batch_vec %>%
  factor() %>%
  as.numeric()

# Hyperparameter for the group weights prior
concentration <- rep(1, K_max)

# === MCMC =====================================================================

# MCMC parameters
R <- 50000
thin <- 100
burn <- 20000
eff_burn <- burn / thin

# Chain used based on complete likelihood trace plots
chain_used <- 3

# MCMC samples
mvt_samples <- readRDS(paste0(
  "./Analysis/Sweden/Outputs/Chains/sweden_mvt_chain_",
  chain_used,
  "_m_scale_1e-01_rho_11_theta_5.rds"
))


# Allocations

# Allocations
if (batchPackageIsCurrent) {
  mvt_samples$R <- R
  mvt_samples$thin <- thin
  
  mvt_prob <- calcAllocProb(mvt_samples, burn)
} else {
  mvt_prob <- calcAllocProb(mvt_samples$alloc, eff_burn)
}
mvt_pred <- predictClass(mvt_prob)


# Inferred datasets
mvt_inferred_data <- rowMeans(mvt_samples$batch_corrected_data[, , -c(1:eff_burn)], dims = 2) %>%
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
    "Model" = "MVT"
  )


p_inferred <- mvt_inferred_data %>%
  ggplot(aes(x = SPIKE, y = RBD))  +
  geom_rect(xmin=-0.5, xmax=0.2, ymin=-0.25, ymax=0.2, color="black", fill = "grey", alpha = 0.2) +
  geom_point(aes(shape = factor(Fixed, label = c("False", "True")),
                 colour = factor(Label,
                                 labels = c("Seronegative", "Seropositive")
                 ),
                 alpha = Prob,
                 size = Fixed) #, size = 1
  ) +
  scale_alpha_continuous(range = c(0.4, 1.0)) +
  # geom_density_2d(aes(colour = Predicted)) +
  ggthemes::scale_color_colorblind() +  
  labs(
    # title = "Elisa data - batch adjusted",
    # subtitle = "Predicted labels",
    # caption = "log transformed"
    x = "Log-transformed batch-corrected SPIKE OD",
    y = "Log-transformed batch-corrected RBD OD",
    colour = "Class",
    shape = "Control",
    alpha = "Probability of\nallocation"
  ) +
  scale_size_discrete(range = c(1.0, 1.6))


# === Patient 4 validation =====================================================

batch_map <- data.frame(
  Orig = m$group,
  New = plot_df$Batch
) %>%
  distinct() %>%
  `[`(order(.$New), )

patient_4_df <- data.frame(
  "SPIKE" = log(patient_4_data$OD.spike),
  "RBD" = log(patient_4_data$OD.rbd),
  "Batch_orig" = patient_4_data$group,
  "Type" = patient_4_data$type
)
patient_4_df$Batch <- batch_map$New[match(patient_4_df$Batch_orig, batch_map$Orig)]

group_means <- mvt_samples$means[, , -c(1:eff_burn)] %>%
  rowMeans(dims = 2L) %>%
  set_colnames(c("Seronegative", "Seropositive"))

batch_scales <- mvt_samples$batch_scale[, , -c(1:eff_burn)] %>%
  rowMeans(dims = 2L)

batch_means <- mvt_samples$batch_shift[, , -c(1:eff_burn)] %>%
  rowMeans(dims = 2L)

transformed_patient_4 <- patient_4_df

transformed_patient_4$Batch_bias_1 <- batch_means[1, patient_4_df$Batch]
transformed_patient_4$Batch_bias_2 <- batch_means[2, patient_4_df$Batch]
transformed_patient_4$Batch_scale_1 <- batch_scales[1, patient_4_df$Batch]
transformed_patient_4$Batch_scale_2 <- batch_scales[2, patient_4_df$Batch]
transformed_patient_4$Group_mean_1 <- group_means[1, 2]
transformed_patient_4$Group_mean_2 <- group_means[2, 2]

head(transformed_patient_4)

transformed_patient_4 %<>%
  mutate(
    SPIKE_corrected = (SPIKE - Batch_bias_1 - Group_mean_1) / Batch_scale_1 + Group_mean_1,
    RBD_corrected = (RBD - Batch_bias_2 - Group_mean_2) / Batch_scale_2 + Group_mean_2
  )

# === Summmary statistics ======================================================
orig_summary <- patient_4_df %>% 
  group_by(Batch) %>% 
  summarise(Mean_SPIKE = mean(SPIKE),
            Mean_RBD = mean(RBD), 
            SD_SPIKE = sd(SPIKE), 
            SD_RBD = sd(RBD))

corrected_summary <- transformed_patient_4 %>% 
  group_by(Batch) %>% 
  summarise(Mean_SPIKE = mean(SPIKE_corrected),
            Mean_RBD = mean(RBD_corrected), 
            SD_SPIKE = sd(SPIKE_corrected), 
            SD_RBD = sd(RBD_corrected))


orig_mean_diff <-sqrt(
  (orig_summary$Mean_SPIKE[1] - orig_summary$Mean_SPIKE[2])**2 +
    (orig_summary$Mean_RBD[1] - orig_summary$Mean_RBD[2])**2
)

orig_sd_ratio <- c(orig_summary$SD_SPIKE[2] / orig_summary$SD_SPIKE[1], 
                   (orig_summary$SD_RBD[2] / orig_summary$SD_RBD[1])
)

new_mean_diff <-sqrt(
  (corrected_summary$Mean_SPIKE[1] - corrected_summary$Mean_SPIKE[2])**2 +
    (corrected_summary$Mean_RBD[1] - corrected_summary$Mean_RBD[2])**2
)

new_sd_ratio <- c(corrected_summary$SD_SPIKE[2] / corrected_summary$SD_SPIKE[1], 
                   (corrected_summary$SD_RBD[2] / corrected_summary$SD_RBD[1])
)

orig_sd_ratio 
new_sd_ratio 

orig_mean_diff
new_mean_diff 


as.matrix(corrected_summary[, c(4, 5)]) / as.matrix(orig_summary[, c(4, 5)])

orig_cov_1 <- patient_4_df %>% 
  filter(Batch == 5) %>% 
  select(SPIKE, RBD) %>% 
  cov()

orig_cov_2 <- patient_4_df %>% 
  filter(Batch == 7) %>% 
  select(SPIKE, RBD) %>% 
  cov()

orig_cov_2[1, 1] <- 0.001
orig_cov_2 <- orig_cov_2 * 0.3

corrected_cov_1 <- transformed_patient_4 %>% 
  filter(Batch == 5) %>% 
  select(SPIKE_corrected, RBD_corrected) %>% 
  cov()

corrected_cov_2 <- transformed_patient_4 %>% 
  filter(Batch == 7) %>% 
  select(SPIKE_corrected, RBD_corrected) %>% 
  cov()


orig_data_el_1 <- confidenceEllipse(mu = c(orig_summary$Mean_SPIKE[1], 
                                           orig_summary$Mean_RBD[1]),
                                    Sigma = orig_cov_1,
                                    confidenceLevel = 0.95)

orig_data_el_2 <- confidenceEllipse(
  mu = c(
    orig_summary$Mean_SPIKE[2], 
    orig_summary$Mean_RBD[2]
  ),
  Sigma = orig_cov_2 + 1e-12 * diag(2),
  confidenceLevel = 0.95)

corrected_data_el_1 <- confidenceEllipse(mu = c(corrected_summary$Mean_SPIKE[1], 
                                                corrected_summary$Mean_RBD[1]),
                                         Sigma = orig_cov_1,
                                         confidenceLevel = 0.95)

corrected_data_el_2 <- confidenceEllipse(
  mu = c(
    corrected_summary$Mean_SPIKE[2], 
    corrected_summary$Mean_RBD[2]
  ),
  Sigma = orig_cov_2,
  confidenceLevel = 0.95)

orig_data_el_1$Batch <- 5
orig_data_el_2$Batch <- 7

corrected_data_el_1$Batch <- 5
corrected_data_el_2$Batch <- 7

orig_ellipse <- rbind(orig_data_el_1, orig_data_el_2)
corrected_ellipse <- rbind(corrected_data_el_1, corrected_data_el_2)

# === Plotting =================================================================

p1 <- patient_4_df %>%
  ggplot() + #, labels = c("Blood Donors", "Pregnant Volunteers")))) +
  geom_point(aes(y = RBD, 
                 x = SPIKE, 
                 colour = factor(Batch), 
                 shape = factor(Batch)),
             size = 3) +
  labs(colour = "Batch",
       shape = "Batch",
       x = "Log-transformed observed SPIKE OD",
       y = "Log-transformed observed RBD OD") +  
  # ggthemes::scale_color_colorblind() +
  xlim(c(-0.5, 0.25)) +
  ylim(c(-0.25, 0.25)) +
  scale_colour_brewer(palette = "Set1")



p1_with_ellipse <- p1 +
  geom_path(orig_ellipse, mapping = aes(x=X1, y=X2, colour = factor(Batch))) 

  # geom_ellipse(data = orig_summary, 
  #             aes(x0 = Mean_SPIKE, 
  #                 y0 = Mean_RBD, 
  #                 a = 1.96 * SD_SPIKE,
  #                 b = 1.96 * SD_RBD,
  #                 colour = factor(Batch),
  #                 angle = 0
  #                 )
  #             )

p2 <- transformed_patient_4 %>%
  ggplot() + #, labels = c("Blood Donors", "Pregnant Volunteers")))) +
  geom_point(aes(y = RBD_corrected, 
                 x = SPIKE_corrected,
                 colour = factor(Batch), 
                 shape = factor(Batch)),
             size = 3) +
  labs(colour = "Batch",
       shape = "Batch",
       x = "Log-transformed batch-corrected SPIKE OD",
       y = "Log-transformed batch-corrected RBD OD") +
  # ggthemes::scale_color_colorblind() +
  xlim(c(-0.5, 0.2)) +
  ylim(c(-0.25, 0.2)) +
  scale_colour_brewer(palette = "Set1")

p2_with_ellipse <- p2 +
  geom_path(corrected_ellipse, mapping = aes(x=X1, y=X2, colour = factor(Batch)))


p_full <- plot_df %>%
  ggplot(aes(x = SPIKE, y = RBD, colour = factor(Label,
    labels = c("Seronegative", "Seropositive", "Unknown")
  ))) +
  geom_rect(xmin=-0.5, xmax=0.2, ymin=-0.25, ymax=0.2, color="black", fill = "grey", alpha = 0.2) +
  geom_point(size = 0.7) +
  # facet_grid(~Controls) +
  labs(
    # title = "ELISA data",
    x = "Log SPIKE OD",
    y = "Log RBD OD",
    colour = "Group"
  ) +
  ggthemes::scale_color_colorblind() 

(p1_with_ellipse | p2_with_ellipse) / p_full +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

p_patch_patient_4 <- (p1_with_ellipse / p2_with_ellipse)  +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

ggsave(paste0(save_dir, "patient4Validation.png"),
       plot = p_patch_patient_4,
       height = 6, 
       width = 7
)
