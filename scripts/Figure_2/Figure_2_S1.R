# -- Figure 2 S1 A-C -----------------------------------------------------------
# ------------------------------------------------------------------------------
# This is the code for plotting DILP counts
# ------------------------------------------------------------------------------

# load packages ----------------------------------------------------------------
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
# When running script for the first time, use settings provided:
write_plots = TRUE           # TRUE - save/replicate figure plots
                             # FALSE - plots not saved outside of R

# load data --------------------------------------------------------------------

DILP2mcherry_DILP2AB = read_delim(file.path(PATH_input,"Figure2_S1A_DILP2mcherry_DILP2AB.csv"), delim = ",")
DILP2mcherry_DILP3AB = read_delim(file.path(PATH_input,"Figure2_S1B_DILP2mcherry_DILP3AB.csv"), delim = ",")
DILP5mcherry_DILP2AB = read_delim(file.path(PATH_input,"Figure2_S1C_DILP5mcherry_DILP2AB.csv"), delim = ",")

# ------------------------------------------------------------------------------
# Figure 2 S1A - DILP2 mcherry and DILP2 antibody 
# ------------------------------------------------------------------------------
panel_name = "_S1A"
filename_dilp2 = paste0(fig_folder_name, panel_name)

# Reshape the data for plotting
df_long <- DILP2mcherry_DILP2AB %>%
  mutate(condition = row_number()) %>%
  pivot_longer(cols = c(mCherry, DILP2), names_to = "value_label", 
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
  scale_x_discrete(limits = c("mCherry", "DILP2"), 
                   labels = c("mCherry" = "DILP2 > mCherry", "DILP2" = "DILP2")) +
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, 5), expand = c(0, 0)) +
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
  annotate("text", x = 1.5, y = 5, label = paste0("(n = ", n_val, ")"),
           size = 4, fontface = "bold", color = "black")

print(lp)

if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename_dilp2,".pdf"))
  ggsave(pdf_path, plot = lp, width = 6, height = 10, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 2 S1B - DILP2 mcherry and DILP3 antibody 
# ------------------------------------------------------------------------------
panel_name = "_S1B"
filename_s1b = paste0(fig_folder_name, panel_name)

# Reshape the data for plotting
df_long <- DILP2mcherry_DILP3AB %>%
  mutate(condition = row_number()) %>%
  pivot_longer(cols = c(mCherry, DILP3), names_to = "value_label", 
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
s1b <- ggplot(df_long, aes(x = value_label, y = value, group = condition)) +
  geom_line(color = "black", linewidth = 0.8) +
  geom_point(size = 3) +
  scale_x_discrete(limits = c("mCherry", "DILP3"), 
                   labels = c("mCherry" = "DILP2 > mCherry", "DILP3" = "DILP3")) +
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, 5), expand = c(0, 0)) +
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
  annotate("text", x = 1.5, y = 5, label = paste0("(n = ", n_val, ")"),
           size = 4, fontface = "bold", color = "black")

print(s1b)

if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename_s1b,".pdf"))
  ggsave(pdf_path, plot = s1b, width = 6, height = 10, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 2 S1C - DILP5 mcherry and DILP2 antibody 
# ------------------------------------------------------------------------------
panel_name = "_S1C"
filename_s1c = paste0(fig_folder_name, panel_name)

# Reshape the data for plotting
df_long <- DILP5mcherry_DILP2AB %>%
  mutate(condition = row_number()) %>%
  pivot_longer(cols = c(mCherry, DILP2), names_to = "value_label", 
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
s1c <- ggplot(df_long, aes(x = value_label, y = value, group = condition)) +
  geom_line(color = "black", linewidth = 0.8) +
  geom_point(size = 3) +
  scale_x_discrete(limits = c("mCherry", "DILP2"), 
                   labels = c("mCherry" = "DILP5 > mCherry", "DILP2" = "DILP2")) +
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
  annotate("text", x = 1.5, y = 5, label = paste0("(n = ", n_val, ")"),
           size = 4, fontface = "bold", color = "black")

print(s1c)

if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename_s1c,".pdf"))
  ggsave(pdf_path, plot = s1c, width = 6, height = 10, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}
# ------------------------------------------------------------------------------
