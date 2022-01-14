


# === Setup ====================================================================
library(BatchMixtureModel)
library(mdiHelpR)
library(ggplot2)
library(magrittr)
library(tibble)
library(dplyr)
library(data.table)


setMyTheme()
set.seed(1)

generateELISAdata <- function(N,
                              P,
                              mu_k,
                              cov_k,
                              df_k,
                              m_b,
                              S_b,
                              class_weights,
                              batch_weights,
                              frac_known = 0.5) {
  
  
  # The number of clusters to generate
  K <- ncol(mu_k)
  
  # The number of batches to generate
  B <- ncol(m_b)
  
  if (ncol(class_weights) != B) {
    stop("Number of columns in class weight matrix does not match the number of batches.")
  }
  
  if (nrow(class_weights) != K) {
    stop("Number of rows in class weight matrix does not match the number of classes.")
  }
  
  # The membership vector for the N points, currently empty
  labels <- rep(0, N)
  
  # The batch labels for the N points
  batches <- sample(seq(1, B), N, replace = T, prob = batch_weights)
  
  # The fixed labels for the semi-supervised case
  fixed <- sample(seq(0, 1), N, replace = T, prob = c(1 - frac_known, frac_known))
  
  # The data matrices
  observed_data <- true_data <- matrix(0, nrow = N, ncol = P)
  
  # Iterate over the batches to sample appropriate labels
  for (b in seq(1, B)) {
    batch_ind <- which(batches == b)
    N_b <- length(batch_ind)
    labels[batch_ind] <- sample(1:K, N_b, replace = T, prob = class_weights[, b])
  }
  
  # Generate the data
  for (n in seq(1, N)) {
    
    # Find the class and batch
    b <- batches[n]
    k <- labels[n]
    
    # Group and batch parameters
    .cov <- cov_k[[k]]
    .mu <- mu_k[, k]
    .df <- df_k[k]
    
    .S <- S_b[, b]
    .m <- m_b[, b]
    
    # Random data
    y <- x <- mvtnorm::rmvt(1, sigma = .cov, delta = .mu, df = .df)
    
    for(p in seq(1, P)) {
      
      # The observed data includes the batch effects
      x[p] <- (x[p] - .mu[p]) * .S[p] + (.mu[p] + .m[p])
      
    }
    
    true_data[n,] <- y
    observed_data[n, ] <- x
    
  }
  
  positives <- which(labels == 2)
  
  fixed[positives] <- makeFixed(true_data[positives, ])
  
  # Return a list of the class labels, batch labels, fixed points and the
  # observed and true datasets.
  list(
    observed_data = observed_data,
    true_data = true_data,
    class_IDs = labels,
    batch_IDs = batches,
    fixed = fixed
  )
}

euclidsDist <- function(x, y) {
  sqrt(sum((x - y)**2))
}

makeFixed <- function(X, dist = "exp") {
  
  mean_x <- colMeans(X)
  max_x <- X %>% apply(2, max)
  
  Y <- apply(X, 1, euclidsDist, max_x)
  
  if(dist == "gamma") {
    probs <- dgamma(Y, shape = 1.0, rate = 0.5, log = F)
  }
  if(dist == "chisq") {
    probs <- dchisq(Y, df = 2, log = F)
  }
  if(dist == "exp") {
    probs <- dexp(Y, rate = 1, log = F)
  }
  
  fixed <- rep(0, nrow(X))
  for(n in 1:nrow(X)) {
    if(probs[n] > 1.0) {
      probs[n] <- 1.0
    }
    
    fixed[n] <- sample(c(0, 1), 1, prob = c(1 - probs[n], probs[n]))
  }
  fixed
}

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
P <- 2 # generate two measurements for each item
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

# === MCMC =====================================================================

# MCMC parameters
R <- 50000
thin <- 100
burn <- 20000
eff_burn <- burn / thin

# Chain used based on complete likelihood trace plots
chain_used <- 3

# MCMC samples
mvt_samples <- readRDS(paste0("./Analysis/Sweden/Outputs/Chains/sweden_mvt_chain_",
                              chain_used, 
                              "_m_scale_1e-01_rho_11_theta_5.rds")
)
# mvn_samples <- readRDS("./Analysis/Outputs/sweden_mvn_chain_4.rds")


mvt_cov <- mvt_samples$covariance[, , -c(1:eff_burn)] %>%
  rowMeans(dims = 2L)

mvt_mu <- mvt_samples$means[, , -c(1:eff_burn)] %>%
  rowMeans(dims = 2L)

mvt_S <- mvt_samples$batch_scale[, , -c(1:eff_burn)] %>%
  rowMeans(dims = 2L)

mvt_m <- mvt_samples$batch_shift[, , -c(1:eff_burn)] %>%
  rowMeans(dims = 2L)

mvt_df <- mvt_samples$t_df[-c(1:eff_burn), ] %>% colMeans()

cov_k <- list(
  mvt_cov[1:2, 1:2],
  mvt_cov[1:2, 3:4]
)

mvt_prob <- calcAllocProb(mvt_samples$alloc, eff_burn)
mvt_pred <- predictClass(mvt_prob)

N_b <- table(batch_vec)

pred_df <- data.frame(
  Group = mvt_pred,
  Batch = batch_vec,
  Fixed = fixed
) %>% group_by(Batch, Group) %>% 
  summarise(Count = n())

class_weights <- matrix(pred_df$Count / N_b[pred_df$Batch], nrow = 2)


batch_weights <- N_b / N

n_sims <- 10

save_dir <- "./Simulations/ELISA_like/Data/"
sims_lst <- vector("list", n_sims)

for(i in 1:n_sims) {

  sim_data <- generateELISAdata(N,
    P,
    mvt_mu,
    cov_k,
    mvt_df,
    mvt_m,
    mvt_S,
    class_weights,
    batch_weights,
    frac_known = 0.33
  )

  new_data_df <- data.frame(
    Labels = factor(sim_data$class_IDs),
    Batch = factor(sim_data$batch_IDs),
    Fixed = sim_data$fixed,
    X_observed = sim_data$observed_data[, 1],
    Y_observed = sim_data$observed_data[, 2],
    X_true = sim_data$true_data[, 1],
    Y_true = sim_data$true_data[, 2]
  )

# new_data_df %>% 
#   filter(Fixed == 1) %>% 
#   ggplot(aes(x = X_observed, y = Y_observed, colour = Labels)) +
#   geom_point() +
#   xlim(c(-6, 2)) +
#   ylim(c(-6, 2))
# 
# new_data_df %>% 
#   # filter(Fixed == 1) %>% 
#   ggplot(aes(x = X_observed, y = Y_observed, colour = Labels)) +
#   geom_point() +
#   xlim(c(-6, 2)) +
#   ylim(c(-6, 2))
  
  fwrite(new_data_df, paste0(save_dir, "elisaLikeSeed", i, ".csv"))
  
  new_data_df$Seed <- i
  
  sims_lst[[i]] <- new_data_df

}


sims_df <- do.call(rbind, sims_lst)


p_all_obs <- sims_df %>% 
  ggplot(aes(x = X_observed, y = Y_observed, colour = factor(Labels))) +
  geom_point(alpha = 0.3) +
  facet_wrap(~Seed) +
  ggthemes::scale_color_colorblind() +
  labs(
    title = "ELISA-like simulations",
    subtitle = "All labels, observed measurements",
    x = "SPIKE (observed)",
    y = "RBD (observed)"
    )

p_fixed_obs <- sims_df %>% 
  filter(Fixed == 1) %>% 
  ggplot(aes(x = X_observed, y = Y_observed, colour = factor(Labels))) +
  geom_point() +
  facet_wrap(~Seed) +
  ggthemes::scale_color_colorblind() +
  labs(
    title = "ELISA-like simulations",
    subtitle = "Observed labels, observed measurements",
    x = "SPIKE (observed)",
    y = "RBD (observed)"
  )

p_all_true <- sims_df %>% 
  ggplot(aes(x = X_true, y = Y_true, colour = factor(Labels))) +
  geom_point(alpha = 0.3) +
  facet_wrap(~Seed) +
  ggthemes::scale_color_colorblind() +
  labs(
    title = "ELISA-like simulations",
    subtitle = "All labels, true measurements",
    x = "SPIKE (batch free)",
    y = "RBD (batch free)"
  )

p_fixed_true <- sims_df %>% 
  filter(Fixed == 1) %>% 
  ggplot(aes(x = X_true, y = Y_true, colour = factor(Labels))) +
  geom_point() +
  facet_wrap(~Seed) +
  ggthemes::scale_color_colorblind() +
  labs(
    title = "ELISA-like simulations",
    subtitle = "Observed labels, true measurements",
    x = "SPIKE (batch free)",
    y = "RBD (batch free)"
  )

ggsave(paste0(save_dir, "observed_data_and_all_labels.png"), 
       plot = p_all_obs,
       height = 10, 
       width = 12)

ggsave(paste0(save_dir, "observed_data_and_observed_labels.png"), 
       plot = p_fixed_obs,
       height = 10, 
       width = 12)

ggsave(paste0(save_dir, "true_data_and_all_labels.png"), 
       plot = p_all_true,
       height = 10, 
       width = 12)

ggsave(paste0(save_dir, "true_data_and_observed_labels.png"), 
       plot = p_fixed_true,
       height = 10, 
       width = 12)
