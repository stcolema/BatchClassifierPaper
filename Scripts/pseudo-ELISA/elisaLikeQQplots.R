
library(ggplot2)
library(magrittr)
library(tidyr)
library(dplyr)
library(tibble)
library(BatchMixtureModel)
library(patchwork)

mdiHelpR::setMyTheme()

set.seed(1)
my_wd <- "./"

# "https://github.com/chr1swallace/seroprevalence-paper/blob/master/adjusted-data.RData"
# Please download before proceeding if it is not already in the Data directory.
# Read in the ELISA data
data_file <- paste0(my_wd, "/Data/Sweden/adjusted-data.RData")
load(data_file)

save_dir <- paste0(my_wd, "Simulations/ELISA_like/Output/QQplots/")

# These samples are poorly behaved, drop
drop_sample_in_12 <- which((m$Sample.ID) %in% 1:2)
m <- m[-drop_sample_in_12, ]

# Drop the "patient 4" samples
patient_4 <- which(m$type == "Patient 4")
m <- m[-patient_4, ]

# Find the controls
negative_controls <- which(m$type == "Historical controls")
positive_controls <- which(m$type == "COVID")
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
Sweden_OD <- log(as.matrix(m[, ..cols_used])) %>%
  set_rownames(paste0("Person_", 1:N)) %>%
  set_colnames(cols_of_interest)



sim_dir <- "./Simulations/ELISA_like/Data/"
data_files <- list.files(sim_dir, pattern = "*.csv") %>%
  stringr::str_sort(numeric = T)

n_sims <- length(data_files) #

for (ii in seq(n_sims)) {
  sim_data <- data.table::fread(paste0(sim_dir, data_files[ii]))
  sim_OD <- as.matrix(sim_data[, c(4, 5)])

  # qqplot(sim_OD[,1], Sweden_OD[,1])
  # qqplot(sim_OD[,2], Sweden_OD[,2])

  # qqplot(sim_OD[,1], Sweden_OD[,1])
  # qqplot(sim_OD[,2], Sweden_OD[,2])

  # qq_df <- data.frame(Sim_SPIKE = sort(sim_OD[,1]),
  #                     Sim_RBD = sort(sim_OD[,2]),
  #                     Sweden_SPIKE = sort(Sweden_OD[,1]),
  #                     Sweden_RBD = sort(Sweden_OD[,2]),
  #                     Simulation = ii)

  Spike_df <- data.frame(
    Sim = sort(sim_OD[, 1]),
    Sweden = sort(Sweden_OD[, 1]),
    Simulation = ii,
    Antigen = "SPIKE"
  )

  RBD_df <- data.frame(
    Sim = sort(sim_OD[, 2]),
    Sweden = sort(Sweden_OD[, 2]),
    Simulation = ii,
    Antigen = "RBD"
  )


  my_df <- rbind(Spike_df, RBD_df)

  if (ii == 1) {
    qq_df <- my_df
  } else {
    qq_df <- rbind(qq_df, my_df)
  }

  # p1 <- qq_df %>%
  #   ggplot(aes(x = Sim_x, y = Sweden_x)) +
  #   geom_point() +
  #   geom_abline(slope = 1.0, intercept = 0.0)
}




p <- ggplot(as.data.frame(sim_OD), aes(x = X_observed, y = Y_observed))
p + stat_qq() + stat_qq_line()


p_all <- qq_df %>%
  ggplot(aes(x = Sim, y = Sweden)) +
  geom_point() +
  geom_abline(slope = 1.0, intercept = 0.0) +
  facet_grid(Simulation ~ Antigen)

p_spike <- qq_df %>%
  filter(Antigen == "SPIKE") %>%
  ggplot(aes(x = Sim, y = Sweden)) +
  geom_point() +
  geom_abline(slope = 1.0, intercept = 0.0) +
  facet_wrap(~Simulation) +
  labs(
    title = "SPIKE",
    x = "Simulation"
  )

p_rbd <- qq_df %>%
  filter(Antigen == "RBD") %>%
  ggplot(aes(x = Sim, y = Sweden)) +
  geom_point() +
  geom_abline(slope = 1.0, intercept = 0.0) +
  facet_wrap(~Simulation) +
  labs(
    title = "RBD",
    x = "Simulation"
  )


ggsave(
  file = paste0(save_dir, "SPIKE.png"),
  plot = p_spike,
  height = 5,
  width = 7
)

ggsave(
  file = paste0(save_dir, "RBD.png"),
  plot = p_rbd,
  height = 5,
  width = 7
)
