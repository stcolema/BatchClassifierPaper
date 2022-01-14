
library(tidyverse)
library(mdiHelpR)

# Set the gglot2 theme
setMyTheme()

# Set the random seed
set.seed(1)

# Save path for plots
save_path <- "./Article/Supplementary_materials/Plots/Simulations/"
dir.create(save_path)

# Where the generated files are saved
file_path <- "./Simulations/Generated_data/"

n_sim <- 10
sim_used <- sample(1:n_sim, size = 1)

scn_df <- readRDS(paste0(file_path, "scenario_descriptions.rds"))
scns <- scn_df$Scenario

data_lst <- vector("list", 6) %>% 
  set_names(scns)

for(scn in scns) {
  curr_dir <- paste0(file_path, scn)
  f <- paste0(curr_dir, "/seed_", sim_used, ".csv")
  my_df <- read.csv(f)
  
  plot_df <- my_df %>%
    pivot_longer(c(X_observed, X_true, Y_observed, Y_true),
                 names_to = c(".value", "Version"),
                 names_sep = "_",
                 names_ptypes = list("Version" = factor())
    )
  
  scn_str <- scn %>% 
    stringr::str_replace_all("_", " ") %>% 
    stringr::str_to_sentence(scn)
  
  # New facet label names for dose variable
  version_labs <- c("Observed", "Batch-corrected")
  names(version_labs) <- c("observed", "true")
  
  
  # Create the plot
  p <- plot_df %>% 
    ggplot(aes(x = X, y= Y, colour = factor(labels), shape = factor(fixed, labels = c("False", "True")))) +
    geom_point(size = 0.8) +
    facet_wrap(~Version, 
               ncol = 1,
               labeller = labeller(Version = version_labs)) +
    labs(
      # title = scn_str,
      # subtitle = "Observed and batch corrected datasets",
      colour = "Class",
      shape = "Observed label"
    ) +
    ggthemes::scale_color_colorblind()
  
  plot_name <- paste0(save_path, scn, "_sim_", sim_used, ".png")
  
  ggsave(plot_name, plot = p, height = 4, width = 6)
  
  plot_df$Scenario <- scn_str
  data_lst[[scn]] <- plot_df
}
