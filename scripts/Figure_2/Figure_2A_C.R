# -- Figure 2 A & C ------------------------------------------------------------
# ------------------------------------------------------------------------------
# Code for plotting DILP and DMS neuron counts
# ------------------------------------------------------------------------------

# load packages-----------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(dplyr)

### general setup --------------------------------------------------------------
options(scipen=999)

PATH_input = "./input"
PATH_output = "./output"

# Check if folders exist, if not, make them
if(!dir.exists(file.path(PATH_input))){
  dir.create(file.path(PATH_input))
  cat("Created input folder at:", normalizePath(PATH_input), "\n")
}
if(!dir.exists(file.path(PATH_output))){
  dir.create(file.path(PATH_output))
  cat("Created output folder at:", normalizePath(PATH_output), "\n")
}

# Folder name for saving figs + check if exists
fig_folder_name = "Figure_2"
if(!dir.exists(file.path(PATH_output, fig_folder_name))){
  dir.create(file.path(PATH_output, fig_folder_name))
  cat("Created figure folder at:", normalizePath(file.path(PATH_output, 
                                                           fig_folder_name)), "\n")
}

# ------------------------------------------------------------------------------
# Options for saving plots & csv files
# ------------------------------------------------------------------------------

write_plots = TRUE           # TRUE - save/replicate figure plots
                             # FALSE - plots not saved outside of R

# load data --------------------------------------------------------------------

DILP = read_delim(file.path(PATH_input,"Figure2A_DILP3mcherry_DILP2AB.csv"), delim = ",")
DMS = read_delim(file.path(PATH_input,"Figure2C_DMSmcherry.csv"), delim = ",")

# ------------------------------------------------------------------------------
# Figure 2A - DILP 
# ------------------------------------------------------------------------------
panel_name = "A"
filename_dilp = paste0(fig_folder_name, panel_name)

# Reshape the data for plotting
df_long <- DILP %>%
  mutate(condition = row_number()) %>%
  pivot_longer(cols = c(mcherry, DILP2), names_to = "value_label", 
               values_to = "value") %>%
  arrange(condition) %>%
  mutate(id = row_number()) %>%
  select(id, condition, value_label, value)

# Calculate means per category
means <- df_long %>%
  group_by(value_label) %>%
  summarise(mean_value = mean(value, na.rm = TRUE))

n_val <- length(unique(df_long$condition)) # sample size

# Line Plot
lp <- ggplot(df_long, aes(x = value_label, y = value, group = condition)) +
  geom_line(color = "black", linewidth = 0.8) +
  geom_point(size = 3) +
  scale_x_discrete(limits = c("mcherry", "DILP2"), 
                   labels = c("mcherry" = "DILP3 > mCherry", "DILP2" = "DILP2")) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5), expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8),
    axis.text.x = element_text(angle = 45, hjust = 1, color = c("magenta", "green"),
                               face = "bold", size = 12),
    axis.text.y = element_text(color = "black", face = "bold", size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 12)) +
  ylab("Number of neurons") +
  annotate("text", x = 1.5, y = 4, label = paste0("(n = ", n_val, ")"),
           size = 4, fontface = "bold", color = "black")

print(lp)

if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename_dilp,".pdf"))
  ggsave(pdf_path, plot = lp, width = 6, height = 10, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 2C - DMS 
# ------------------------------------------------------------------------------
panel_name = "C"
filename_dms = paste0(fig_folder_name, panel_name)

# info
summary(DMS$DMS_mcherry)

n_val <- nrow(DMS)  # sample size

plot2C <- ggplot(DMS, aes(x = "DMS_mcherry", y = DMS_mcherry)) +
  geom_jitter(width = 0.1, size = 2, color = "magenta") +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, fatten = 3, color = "black") +
  scale_x_discrete(labels = c("DMS_mcherry" = "DMS > mCherry")) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 2), expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8),
    axis.text.x = element_text(color = c("magenta"), face = "bold", size = 12),
    axis.text.y = element_text(color = "black", face = "bold", size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 12)
  ) +
  ylab("Number of neurons") + 
  annotate("text", x = 1, y = 2, label = paste0("(n = ", n_val, ")"), 
           size = 4, fontface = "bold", color = "black")

print(plot2C)

if (write_plots) {
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename_dms, ".pdf"))
  ggsave(pdf_path,plot = plot2C, width = 5, height = 7, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}
# ------------------------------------------------------------------------------
