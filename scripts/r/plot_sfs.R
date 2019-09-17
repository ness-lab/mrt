# Load required packages
library(tidyverse)

# Load in data
sfs_dat <- read_csv('Google Drive/SLiM_Theta_Rho/N100_bot0.0_gen400_sfs.csv')

sfs_plot <- sfs_dat %>% 
  filter(X1 != "0") %>% 
  ggplot(., aes(x = X1, y = sfs)) + 
  geom_bar(stat = 'identity') + 
  theme_bw()
sfs_plot
