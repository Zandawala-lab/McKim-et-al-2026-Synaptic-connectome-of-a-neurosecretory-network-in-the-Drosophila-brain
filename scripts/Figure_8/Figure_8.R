# -- Figure 8 ------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Code for plotting and analysis of total feeding (CAFE), survival, & egg-laying
# ------------------------------------------------------------------------------

# load packages ----------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(dplyr)
library(survival)
library(survminer)
library(rstatix)
library(ggpubr)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggsignif)
library(broom)


### general setup --------------------------------------------------------------
options(scipen=999)

PATH_input = "./input"
PATH_output = "./output"

input_files = list.files(path = PATH_input, full.names = FALSE, recursive = FALSE)
# Check if 'tmp' folder exists inside PATH_input, if not, make it
tmp_path <- file.path(PATH_input, "tmp")
if(!dir.exists(tmp_path)){
  dir.create(tmp_path)
  cat("Created tmp folder at:", normalizePath(tmp_path), "\n")
}
input_files_tmp = list.files(file.path(path = paste0(PATH_input,"/tmp/")),
                             full.names = FALSE, recursive = FALSE)
input_files = c(input_files,input_files_tmp)

v = read_delim(file.path(PATH_input,"version.csv"),
               col_types  =  cols(version  =  col_character()),delim  =  ";")
v = v$version[1]

# Folder name for saving figs + check if exists
fig_folder_name = "Figure_8"
if(!dir.exists(file.path(PATH_output, fig_folder_name))){
  dir.create(file.path(PATH_output, fig_folder_name))
  cat("Created figure folder at:", normalizePath(file.path(PATH_output, 
                                                           fig_folder_name)), "\n")
}

# ------------------------------------------------------------------------------
# Options for saving plots & csv files
# ------------------------------------------------------------------------------
# When running script for the first time, set both to TRUE:
write_plots = TRUE           # TRUE - save/replicate figure plots
                             # FALSE - plots not saved outside of R
write_csv   = TRUE           # TRUE - save processed data associated w/figures
                             # FALSE - data not saved outside of R


# load data --------------------------------------------------------------------

DNg27_feed = read_delim(file.path(PATH_input, "Figure8_CAFEassay_DNg27.csv"), delim = ",")
CRZ_feed = read_delim(file.path(PATH_input, "Figure8_CAFEassay_CRZ.csv"), delim = ",")

DNg27_survival = read_delim(file.path(PATH_input, "Figure8_Survival_DNg27.csv"), delim = ",")
CRZ_survival = read_delim(file.path(PATH_input, "Figure8_Survival_CRZ.csv"), delim = ",")

CRZ_eggs = read_delim(file.path(PATH_input, "Figure8_Egglaying_CRZ.csv"), delim = ",")
DNg27_eggs = read_delim(file.path(PATH_input, "Figure8_Egglaying_DNg27.csv"), delim = ",")
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Fig 8A - CRZ - CAFE assay - feeding
#-------------------------------------------------------------------------------
panel_name = "A"
filename = paste0(fig_folder_name, panel_name)

# pivot df
CRZ_long <- pivot_longer(CRZ_feed, 
                           cols = everything(),
                           names_to = "Genotype", 
                           values_to = "Feeding")

# add factors for order of plot
CRZ_long$Genotype <- factor(CRZ_long$Genotype, 
                              levels = c("Crz > w1118",
                                         "w1118 > Kir2.1",
                                         "Crz > Kir2.1"))

# normality check
CRZ_long %>%
  group_by(Genotype) %>%
  shapiro_test(Feeding)

ggqqplot(CRZ_long, x = "Feeding", facet.by = "Genotype")

CRZ_long %>%
  levene_test(Feeding ~ Genotype)

# one-way anova
anova_result <- aov(Feeding ~ Genotype, data = CRZ_long)
summary(anova_result)

# post hoc tests 
stat.test <- CRZ_long %>%
  tukey_hsd(Feeding ~ Genotype) %>%
  add_xy_position(x = "Genotype")

# plot
crz <- ggplot(CRZ_long, aes(x = Genotype, y = Feeding, color = Genotype)) +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.5) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.8) +
  scale_color_manual(
    values = c(
      "Crz > w1118"  = "black",
      "w1118 > Kir2.1" = "gray40",
      "Crz > Kir2.1" = "blue"
    ),
    labels = c(
      "Crz > w1118"  = expression(italic("Crz > w"^"1118")),
      "w1118 > Kir2.1" = expression(italic("w"^"1118"*" > Kir2.1")),
      "Crz > Kir2.1" = expression(italic("Crz > Kir2.1"))
    )
  ) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 1),
                     expand = c(0, 0)) +
  labs(y = "Total Feeding (µL/fly)", color = NULL) +
  guides(color = guide_legend(override.aes = list(linetype = 0,
                                                  shape = 16, size = 3))) +
  theme_classic() +
  theme(
    panel.background  = element_rect(fill = "transparent", color = NA),
    plot.background   = element_rect(fill = "transparent", color = NA),
    panel.grid        = element_blank(),
    axis.line.x       = element_blank(),
    axis.line.y       = element_line(color = "black"),
    axis.ticks.x      = element_blank(),
    axis.text.x       = element_blank(),
    axis.title.x      = element_blank(),
    axis.text.y       = element_text(size = 10, color = "black"),
    axis.title.y      = element_text(size = 10, color = "black"),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.text       = element_text(size = 8, color = "black"),
    legend.key.size    = unit(0.4, "cm"),
    legend.spacing.y   = unit(0.2, "cm"),
    legend.title      = element_blank()
  ) + 

  stat_pvalue_manual(
    stat.test,
    label = "p.adj.signif",
    tip.length = 0.01,
    hide.ns = TRUE
  )
print(crz)

if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename, ".pdf"))
  ggsave(pdf_path, plot = crz, width = 12, height = 6, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
  
  png_path <- file.path(PATH_output, fig_folder_name, paste0(filename, ".png"))
  ggsave(png_path, plot = crz, width = 12, height = 6, units = "cm",
         bg = "transparent")
  cat("Saved PNG to:", normalizePath(png_path), "\n")
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Fig 8D - DNg27 - CAFE assay - feeding
#-------------------------------------------------------------------------------
panel_name = "D"
filename_d = paste0(fig_folder_name, panel_name)

# pivot df
DNg27_long <- pivot_longer(DNg27_feed, 
                           cols = everything(),
                           names_to = "Genotype", 
                           values_to = "Feeding")

# add factors for order of plot
DNg27_long$Genotype <- factor(DNg27_long$Genotype, 
                              levels = c("DNg27 > w1118",
                                         "w1118 > Kir2.1",
                                         "DNg27 > Kir2.1"))

# normality check
DNg27_long %>%
  group_by(Genotype) %>%
  shapiro_test(Feeding)

ggqqplot(DNg27_long, x = "Feeding", facet.by = "Genotype")

DNg27_long %>%
  levene_test(Feeding ~ Genotype)

# one-way anova
anova_result <- aov(Feeding ~ Genotype, data = DNg27_long)
# not statistically different
summary(anova_result)

# plot
dng27 <- ggplot(DNg27_long, aes(x = Genotype, y = Feeding, color = Genotype)) +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.5) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.8) +
  scale_color_manual(
    values = c(
      "DNg27 > w1118"  = "black",
      "w1118 > Kir2.1" = "gray40",
      "DNg27 > Kir2.1" = "red"
    ),
    labels = c(
      "DNg27 > w1118"  = expression(italic("DNg27 > w"^"1118")),
      "w1118 > Kir2.1" = expression(italic("w"^"1118"*" > Kir2.1")),
      "DNg27 > Kir2.1" = expression(italic("DNg27 > Kir2.1"))
    )
  ) +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 1),
                     expand = c(0, 0)) +
  labs(y = "Total Feeding (µL/fly)", color = NULL) +
  guides(color = guide_legend(override.aes = list(linetype = 0,
                                                  shape = 16, size = 3))) +
  theme_classic() +
  theme(
    panel.background  = element_rect(fill = "transparent", color = NA),
    plot.background   = element_rect(fill = "transparent", color = NA),
    panel.grid        = element_blank(),
    axis.line.x       = element_blank(),
    axis.line.y       = element_line(color = "black"),
    axis.ticks.x      = element_blank(),
    axis.text.x       = element_blank(),
    axis.title.x      = element_blank(),
    axis.text.y       = element_text(size = 10, color = "black"),
    axis.title.y      = element_text(size = 10, color = "black"),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.text       = element_text(size = 8, color = "black"),
    legend.key.size    = unit(0.4, "cm"),
    legend.spacing.y   = unit(0.2, "cm"),
    legend.title      = element_blank()
  )
print(dng27)

if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename_d, ".pdf"))
  ggsave(pdf_path, plot = dng27, width = 12, height = 6, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
  
  png_path <- file.path(PATH_output, fig_folder_name, paste0(filename_d, ".png"))
  ggsave(png_path, plot = dng27, width = 12, height = 6, units = "cm",
         bg = "transparent")
  cat("Saved PNG to:", normalizePath(png_path), "\n")
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Fig 8B - CRZ - Survival Visualization
#-------------------------------------------------------------------------------
panel_name = "B"
filename_b = paste0(fig_folder_name, panel_name)

CRZ_surv_long <- CRZ_survival %>%
  pivot_longer(cols = starts_with("hr_"),
               names_to = "hour",
               values_to = "dead") %>%
  mutate(hour = as.numeric(gsub("hr_", "", hour)),
         alive = 20 - dead,
         pct_survival = (alive / 20) * 100)

CRZ_summary <- CRZ_surv_long %>%
  group_by(Genotype, hour) %>%
  summarise(
    mean_survival = mean(pct_survival),
    se = sd(pct_survival) / sqrt(n()),
    .groups = "drop"
  )

# add starting 0 time
zero_time <- data.frame(
  Genotype = c("Crz > w1118",
               "w1118 > Kir2.1",
               "Crz > Kir2.1"),
  hour = 0,
  mean_survival = 100,
  se = 0
)

# combine dfs
CRZ_summary <- bind_rows(zero_time, CRZ_summary)

# set factor order for plotting
CRZ_summary$Genotype <- factor(CRZ_summary$Genotype, 
                               levels = c("Crz > w1118",
                                          "w1118 > Kir2.1",
                                          "Crz > Kir2.1"))

# ggplot survival curve
crz_surv <- ggplot(CRZ_summary, aes(x = hour, y = mean_survival,
                                    color = Genotype)) +
  geom_line(size = 0.8) +
  geom_ribbon(aes(ymin = mean_survival - se, ymax = mean_survival + se,
                  fill = Genotype), alpha = 0.15, color = NA) +
  annotate("text", 
           x = 65, 
           y = 50, 
           label = "****", 
           color = "black", 
           size = 3.5, 
           hjust = 0.5) +
  scale_color_manual(values = c(
    "Crz > w1118"    = "black",
    "w1118 > Kir2.1" = "gray40",
    "Crz > Kir2.1"   = "blue"
  ),
  labels = c(
    "Crz > w1118"  = expression(italic("Crz > w"^"1118")),
    "w1118 > Kir2.1" = expression(italic("w"^"1118"*" > Kir2.1")),
    "Crz > Kir2.1" = expression(italic("Crz > Kir2.1"))
  )) +
  scale_fill_manual(values = c(
    "Crz > w1118"    = "black",
    "w1118 > Kir2.1" = "gray40",
    "Crz > Kir2.1"   = "blue"
  ),
  labels = c(
    "Crz > w1118"  = expression(italic("Crz > w"^"1118")),
    "w1118 > Kir2.1" = expression(italic("w"^"1118"*" > Kir2.1")),
    "Crz > Kir2.1" = expression(italic("Crz > Kir2.1"))
  )) +
  scale_x_continuous(
    limits = c(0, 108),
    breaks = sort(unique(c(0, seq(12, 108, by = 12), 108))),
    expand = c(0, 0)) +  
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(x = "Time (hours)", y = "Survival (%)", color = NULL, fill = NULL) +
  theme_classic() +
  theme(
    panel.background  = element_rect(fill = "transparent", color = NA),
    plot.background   = element_rect(fill = "transparent", color = NA),
    panel.grid        = element_blank(),
    axis.line         = element_line(color = "black"),
    axis.text         = element_text(size = 8, color = "black"),
    axis.title        = element_text(size = 8, color = "black"),
    legend.text       = element_text(size = 8, color = "black"),
    legend.title      = element_blank(),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key        = element_blank(),
    legend.key.size    = unit(0.4, "cm"),
    legend.spacing.y   = unit(0.2, "cm")
  )
print(crz_surv)

if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename_b, ".pdf"))
  ggsave(pdf_path, plot = crz_surv, width = 14, height = 6, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
  
  png_path <- file.path(PATH_output, fig_folder_name, paste0(filename_b, ".png"))
  ggsave(png_path, plot = crz_surv, width = 14, height = 6, units = "cm",
         bg = "transparent")
  cat("Saved PNG to:", normalizePath(png_path), "\n")
}

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Fig 8E - DNg27 - Survival Visualization
#-------------------------------------------------------------------------------
panel_name = "E"
filename_e = paste0(fig_folder_name, panel_name)

# reshape the data
DNg27_surv_long <- DNg27_survival %>%
  pivot_longer(cols = starts_with("hr_"),
               names_to = "hour",
               values_to = "dead") %>%
  mutate(hour = as.numeric(gsub("hr_", "", hour)),
         alive = 20 - dead,
         pct_survival = (alive / 20) * 100)

# summarize
DNg27_summary <- DNg27_surv_long %>%
  group_by(Genotype, hour) %>%
  summarise(
    mean_survival = mean(pct_survival),
    se = sd(pct_survival) / sqrt(n()),
    .groups = "drop"
  )

# add starting 0 time
zero_time <- data.frame(
  Genotype = c("DNg27 > w1118",
               "w1118 > Kir2.1",
               "DNg27 > Kir2.1"),
  hour = 0,
  mean_survival = 100,
  se = 0
)

# combine dfs
DNg27_summary <- bind_rows(zero_time, DNg27_summary)

DNg27_summary$hour <- as.numeric(DNg27_summary$hour)

# set factor order for plotting
DNg27_summary$Genotype <- factor(DNg27_summary$Genotype, 
                               levels = c("DNg27 > w1118",
                                          "w1118 > Kir2.1",
                                          "DNg27 > Kir2.1"))

# ggplot survival curve
dng27_surv <- ggplot(DNg27_summary, aes(x = hour, y = mean_survival,
                                        color = Genotype)) +
  geom_line(size = 0.8) +
  geom_ribbon(aes(ymin = mean_survival - se, ymax = mean_survival + se,
                  fill = Genotype), alpha = 0.15, color = NA) +
  scale_color_manual(values = c(
    "DNg27 > w1118"    = "black",
    "w1118 > Kir2.1" = "gray40",
    "DNg27 > Kir2.1"   = "red"
  ),
  labels = c(
    "DNg27 > w1118"  = expression(italic("DNg27 > w"^"1118")),
    "w1118 > Kir2.1" = expression(italic("w"^"1118"*" > Kir2.1")),
    "DNg27 > Kir2.1" = expression(italic("DNg27 > Kir2.1"))
  )) +
  scale_fill_manual(
    values = c(
      "DNg27 > w1118"  = "black",
      "w1118 > Kir2.1" = "gray40",
      "DNg27 > Kir2.1" = "red"
    ),
    labels = c(
      "DNg27 > w1118"  = expression(italic("DNg27 > w"^"1118")),
      "w1118 > Kir2.1" = expression(italic("w"^"1118"*" > Kir2.1")),
      "DNg27 > Kir2.1" = expression(italic("DNg27 > Kir2.1"))
    )
  ) +
  #scale_x_continuous(breaks = c(0, seq(24, 100, by = 8))) +
  scale_x_continuous(
    limits = c(0, 108),
    breaks = sort(unique(c(0, seq(12, 108, by = 12), 108))),
    expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(x = "Time (hours)", y = "Survival (%)", color = NULL, fill = NULL) +
  theme_classic() +
  theme(
    panel.background  = element_rect(fill = "transparent", color = NA),
    plot.background   = element_rect(fill = "transparent", color = NA),
    panel.grid        = element_blank(),
    axis.line         = element_line(color = "black"),
    axis.text         = element_text(size = 8, color = "black"),
    axis.title        = element_text(size = 8, color = "black"),
    legend.text       = element_text(size = 9, color = "black"),
    legend.title      = element_blank(),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key        = element_blank(),
    legend.key.size    = unit(0.4, "cm"),
    legend.spacing.y   = unit(0.2, "cm")
  )

print(dng27_surv)

if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename_e, ".pdf"))
  ggsave(pdf_path, plot = dng27_surv, width = 14, height = 6, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
  
  png_path <- file.path(PATH_output, fig_folder_name, paste0(filename_e, ".png"))
  ggsave(png_path, plot = dng27_surv, width = 14, height = 6, units = "cm",
         bg = "transparent")
  cat("Saved PNG to:", normalizePath(png_path), "\n")
}

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#                Survminer package for statistics for survival plots
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Fig 8B - CRZ - Survival Analysis
#-------------------------------------------------------------------------------
# convert cumulative dead to individual-level event data
# each replicate needs one row per fly with time of death
surv_individual <- CRZ_surv_long %>%
  arrange(Genotype, Replicate, hour) %>%
  group_by(Genotype, Replicate) %>%
  mutate(new_deaths = dead - lag(dead, default = 0)) %>%
  ungroup()

# expand to one row per fly
events <- surv_individual %>%
  filter(new_deaths > 0) %>%
  uncount(new_deaths) %>%
  mutate(status = 1)  # 1 = died

# add censored flies (still alive at hr_100)
censored <- CRZ_surv_long %>%
  group_by(Genotype, Replicate) %>%
  slice_max(hour) %>%
  mutate(
    alive = 20 - dead,
    status = 0   # 0 = censored
  ) %>%
  filter(alive > 0) %>%
  uncount(alive) %>%
  ungroup()

# combine
surv_data <- bind_rows(
  events %>% select(Genotype, Replicate, hour, status),
  censored %>% select(Genotype, Replicate, hour, status)
)

surv_data$Genotype <- factor(surv_data$Genotype,
                             levels = c("Crz > w1118",
                                        "w1118 > Kir2.1",
                                        "Crz > Kir2.1"))

# fit KM curves
fit <- survfit(Surv(hour, status) ~ Genotype, data = surv_data)

# plot
ggsurvplot(
  fit,
  data        = surv_data,
  palette     = c("black", "gray40", "blue"), 
  size        = 0.8,
  censor      = TRUE,
  censor.size = 3,
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.y.text.col = T, # colour risk table text annotations
  risk.table.y.text = FALSE,
  xlab        = "Hour",
  ylab        = "Survival probability",
  xlim = c(0,100), 
  break.time.by = 12, 
  legend.title = "",
  ggtheme     = theme_classic()
)

# overall log-rank test
surv_diff <- survdiff(Surv(hour, status) ~ Genotype, data = surv_data)

# pairwise
surv_pairwise <- pairwise_survdiff(Surv(hour, status) ~ Genotype, 
                                   data = surv_data, 
                                   p.adjust.method = "bonferroni")

# median survival table
surv_table <- summary(fit)$table

# save
txt_path <- file.path(PATH_output, fig_folder_name, 
                      paste0(filename_b, "_CRZ_survival.txt"))
sink(txt_path)

cat("=== Overall Log-Rank Test ===\n")
print(surv_diff)

cat("\n=== Pairwise Log-Rank Tests (Bonferroni-adjusted) ===\n")
print(surv_pairwise)

cat("\n=== Median Survival Table ===\n")
print(surv_table)

sink()
cat("Saved TXT to:", normalizePath(txt_path), "\n")

# csv
# pairwise p-values
pairwise_csv_path <- file.path(PATH_output, fig_folder_name,
                               "Figure_8B_CRZ_survival_pairwise.csv")
write.csv(
  as.data.frame(surv_pairwise$p.value),
  pairwise_csv_path
)
cat("Saved CSV to:", normalizePath(pairwise_csv_path), "\n")

# median survival table
median_csv_path <- file.path(PATH_output, fig_folder_name,
                             "Figure_8B_CRZ_survival_median.csv")
write.csv(
  as.data.frame(surv_table),
  median_csv_path
)
cat("Saved CSV to:", normalizePath(median_csv_path), "\n")
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Fig 8E - DNg27 - Survival Analysis
#-------------------------------------------------------------------------------

# convert cumulative dead to individual-level event data
# each replicate needs one row per fly with time of death
surv_individual <- DNg27_surv_long %>%
  arrange(Genotype, Replicate, hour) %>%
  group_by(Genotype, Replicate) %>%
  mutate(new_deaths = dead - lag(dead, default = 0)) %>%
  ungroup()

# expand to one row per fly
events <- surv_individual %>%
  filter(new_deaths > 0) %>%
  uncount(new_deaths) %>%
  mutate(status = 1)  # 1 = died

# add censored flies (still alive at hr_100)
censored <- DNg27_surv_long %>%
  group_by(Genotype, Replicate) %>%
  slice_max(hour) %>%
  mutate(
    alive = 20 - dead,
    status = 0   # 0 = censored
  ) %>%
  filter(alive > 0) %>%
  uncount(alive) %>%
  ungroup()

# combine
surv_data <- bind_rows(
  events %>% select(Genotype, Replicate, hour, status),
  censored %>% select(Genotype, Replicate, hour, status)
)

surv_data$Genotype <- factor(surv_data$Genotype,
                             levels = c("DNg27 > w1118",
                                        "w1118 > Kir2.1",
                                        "DNg27 > Kir2.1"))

# fit KM curves
fit <- survfit(Surv(hour, status) ~ Genotype, data = surv_data)

# plot
ggsurvplot(
  fit,
  data        = surv_data,
  palette     = c("black", "gray40", "red"), #check order
  size        = 0.8,
  censor      = TRUE,
  censor.size = 3,
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE,
  xlab        = "Hour",
  ylab        = "Survival probability",
  xlim = c(0,100), 
  break.time.by = 12, 
  legend.title = "",
  ggtheme     = theme_classic()
)

surv_diff <- survdiff(Surv(hour, status) ~ Genotype, data = surv_data)

# pairwise comparisons
surv_pairwise <- pairwise_survdiff(Surv(hour, status) ~ Genotype,
                                   data = surv_data, p.adjust.method = "bonferroni")

surv_table <- summary(fit)$table

# save
txt_path <- file.path(PATH_output, fig_folder_name, paste0(filename_e,
                                                           "_DNg27_survival.txt"))
sink(txt_path)

cat("=== Overall Log-Rank Test ===\n")
print(surv_diff)

cat("\n=== Pairwise Log-Rank Tests (Bonferroni-adjusted) ===\n")
print(surv_pairwise)

cat("\n=== Median Survival Table ===\n")
print(surv_table)

sink()
cat("Saved TXT to:", normalizePath(txt_path), "\n")

# csv
# pairwise p-values
pairwise_csv_path <- file.path(PATH_output, fig_folder_name,
                               "Figure_8E_DNg27_survival_pairwise.csv")
write.csv(
  as.data.frame(surv_pairwise$p.value),
  pairwise_csv_path
)
cat("Saved CSV to:", normalizePath(pairwise_csv_path), "\n")

# median survival table
median_csv_path <- file.path(PATH_output, fig_folder_name,
                             "Figure_8E_DNg27_survival_median.csv")
write.csv(
  as.data.frame(surv_table),
  median_csv_path
)
cat("Saved CSV to:", normalizePath(median_csv_path), "\n")

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Fig 8C - CRZ Egg-laying
#-------------------------------------------------------------------------------
panel_name = "C"
filename_c = paste0(fig_folder_name, panel_name)

# Pivot to long format
CRZ_eggs_long <- CRZ_eggs %>%
  pivot_longer(
    cols = starts_with("day_"),
    names_to = "day",
    values_to = "eggs"
  ) %>%
  mutate(day = as.integer(str_remove(day, "day_")))

CRZ_eggs_long <- CRZ_eggs_long %>%
  mutate(UniqueRep = interaction(Genotype, Replicate))

# add factors for order of plot
CRZ_eggs_long$Genotype <- factor(CRZ_eggs_long$Genotype, 
                                 levels = c("Crz > w1118",
                                            "w1118 > Kir2.1",
                                            "Crz > Kir2.1"))

# rm anova stats
m <- lmer(eggs ~ Genotype * factor(day) + (1 | UniqueRep), data = CRZ_eggs_long)
summary(m)
anova(m)            # Type III F-tests for fixed effects

# by time effects
emm <- emmeans(m, ~ Genotype | factor(day))   # genotype differences at each day
pairs_by_day <- pairs(emm, adjust = "tukey")
print(pairs_by_day)

# save
# full text summary (model + anova + pairwise)
txt_path <- file.path(PATH_output, fig_folder_name, 
                      paste0(filename_c, "_CRZ_eggs_statistics.txt"))
sink(txt_path)

cat("=== Linear Mixed Model Summary ===\n")
print(summary(m))

cat("\n=== ANOVA (Type III F-tests) ===\n")
print(anova(m))

cat("\n=== Time effects - Tukey-adjusted) ===\n")
print(pairs_by_day)

sink()
cat("Saved TXT to:", normalizePath(txt_path), "\n")

# csv tables
anova_csv_path <- file.path(PATH_output, fig_folder_name, 
                            paste0(filename_c, "_CRZ_eggs_anova.csv"))
write.csv(
  as.data.frame(anova(m)),
  anova_csv_path
)
cat("Saved CSV to:", normalizePath(anova_csv_path), "\n")

pairwise_csv_path <- file.path(PATH_output, fig_folder_name, 
                               paste0(filename_c, "_CRZ_eggs_pairwise.csv"))
write.csv(
  as.data.frame(pairs_by_day),
  pairwise_csv_path
)
cat("Saved CSV to:", normalizePath(pairwise_csv_path), "\n")

# mean ± SEM
mean_se_custom <- function(x) {
  m  <- mean(x, na.rm = TRUE)
  se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  data.frame(y = m, ymin = m - se, ymax = m + se)
}
# color map
geno_colors <- c(
  "Crz > w1118"    = "black",
  "w1118 > Kir2.1" = "gray40",
  "Crz > Kir2.1"   = "blue"
)

# Plot
crz_eggs <- ggplot(
  CRZ_eggs_long,
  aes(x = day, y = eggs, colour = Genotype, group = Genotype)
) +
  
  # SEM error bars at each timepoint
  stat_summary(
    fun.data = mean_se,
    geom     = "errorbar",
    width    = 0.2,
    linewidth = 0.6,
    na.rm    = TRUE
  ) +
  
  # mean dot at each timepoint
  stat_summary(
    fun   = mean,
    geom  = "point",
    size  = 2.5,
    shape = 16,
    na.rm = TRUE
  ) +
  
  # line connecting the means
  stat_summary(
    fun      = mean,
    geom     = "line",
    linewidth = 0.8,
    na.rm    = TRUE
  ) +
  scale_x_continuous(breaks = 1:14, expand = expansion(add = c(0.5, 0.5))) +
  scale_y_continuous(breaks = c(0, 40, 80, 120, 160), expand = c(0, 0)) +
  coord_cartesian(xlim = c(1, 14), ylim = c(0, 160), clip = "off") + 
  scale_colour_manual(values = c(
    "Crz > w1118"    = "black",
    "w1118 > Kir2.1" = "gray40",
    "Crz > Kir2.1"   = "blue"
  ),
  labels = c(
    "Crz > w1118"  = expression(italic("Crz > w"^"1118")),
    "w1118 > Kir2.1" = expression(italic("w"^"1118"*" > Kir2.1")),
    "Crz > Kir2.1" = expression(italic("Crz > Kir2.1"))
  )) +
  labs(x = "Day", y = "Egg count") +
  theme_classic() +
  theme(
    panel.background  = element_rect(fill = "transparent", color = NA),
    plot.background   = element_rect(fill = "transparent", color = NA),
    panel.grid        = element_blank(),
    axis.line         = element_line(color = "black"),   # both axes
    axis.text.x       = element_text(size = 10, color = "black"),
    axis.title.x      = element_text(size = 10, color = "black"),
    axis.ticks.x      = element_line(color = "black"),
    axis.text.y       = element_text(size = 10, color = "black"),
    axis.title.y      = element_text(size = 10, color = "black"),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.text       = element_text(size = 8, color = "black"),
    legend.key.size   = unit(0.4, "cm"),
    legend.spacing.y  = unit(0.2, "cm"),
    legend.title      = element_blank()
  )

# add statistics
crz_eggs_stats <- crz_eggs +
  annotate("segment", x = 4, xend = 14, y = 25, yend = 25, color = "black",
           linewidth = 0.6) +
  annotate("text", x = 9, y = 20, label = "***", size = 4, color = "black",
           hjust = 0.5)

print(crz_eggs_stats)

if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename_c, ".pdf"))
  ggsave(pdf_path, plot = crz_eggs_stats, width = 16, height = 8, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
  
  png_path <- file.path(PATH_output, fig_folder_name, paste0(filename_c, ".png"))
  ggsave(png_path, plot = crz_eggs_stats, width = 16, height = 8, units = "cm", 
         bg = "transparent")
  cat("Saved PNG to:", normalizePath(png_path), "\n")
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Fig 8F - DNg27 Egg-laying
#-------------------------------------------------------------------------------
panel_name = "F"
filename_f = paste0(fig_folder_name, panel_name)

# Pivot to long format
DNg27_eggs_long <- DNg27_eggs %>%
  pivot_longer(
    cols = starts_with("day_"),
    names_to = "day",
    values_to = "eggs"
  ) %>%
  mutate(day = as.integer(str_remove(day, "day_")))

DNg27_eggs_long <- DNg27_eggs_long %>%
  mutate(UniqueRep = interaction(Genotype, Replicate))

# add factors for order of plot
DNg27_eggs_long$Genotype <- factor(DNg27_eggs_long$Genotype, 
                                   levels = c("DNg27 > w1118", 
                                              "w1118 > Kir2.1", 
                                              "DNg27 > Kir2.1"))

# rm anova stats
e <- lmer(eggs ~ Genotype * factor(day) + (1 | UniqueRep), data = DNg27_eggs_long)
summary(e)
anova(e)            # Type III F-tests for fixed effects

# by time effects
emm_e <- emmeans(e, ~ Genotype | factor(day))   # genotype differences at each day
pairs_by_day <- pairs(emm_e, adjust = "tukey")   
print(pairs_by_day)

# save
# full text summary (model + anova + pairwise)
txt_path <- file.path(PATH_output, fig_folder_name, 
                      paste0(filename_f, "_DNg27_eggs_statistics.txt"))
sink(txt_path)

cat("=== Linear Mixed Model Summary ===\n")
print(summary(e))

cat("\n=== ANOVA (Type III F-tests) ===\n")
print(anova(e))

cat("\n=== Time effects - Tukey-adjusted) ===\n")
print(pairs_by_day)

sink()
cat("Saved TXT to:", normalizePath(txt_path), "\n")


# csv tables
anova_csv_path <- file.path(PATH_output, fig_folder_name, 
                            paste0(filename_f, "_DNg27_eggs_anova.csv"))
write.csv(
  as.data.frame(anova(m)),
  anova_csv_path
)
cat("Saved CSV to:", normalizePath(anova_csv_path), "\n")

pairwise_csv_path <- file.path(PATH_output, fig_folder_name, 
                               paste0(filename_f, "_DNg27_eggs_pairwise.csv"))
write.csv(
  as.data.frame(pairs_by_day),
  pairwise_csv_path
)
cat("Saved CSV to:", normalizePath(pairwise_csv_path), "\n")


# color map
geno_colors <- c(
  "DNg27 > w1118"    = "black",
  "w1118 > Kir2.1" = "gray40",
  "DNg27 > Kir2.1"   = "red"
)

# Plot
dng27_eggs <- ggplot(
  DNg27_eggs_long,
  aes(x = day, y = eggs, colour = Genotype, group = Genotype)
) +
  
  # SEM error bars at each timepoint
  stat_summary(
    fun.data = mean_se,
    geom     = "errorbar",
    width    = 0.2,
    linewidth = 0.6,
    na.rm    = TRUE
  ) +
  
  # mean dot at each timepoint
  stat_summary(
    fun   = mean,
    geom  = "point",
    size  = 2.5,
    shape = 16,
    na.rm = TRUE
  ) +
  
  # line connecting the means
  stat_summary(
    fun      = mean,
    geom     = "line",
    linewidth = 0.8,
    na.rm    = TRUE
  ) +
  scale_x_continuous(breaks = 1:14, expand = expansion(add = c(0.5, 0.5))) +
  scale_y_continuous(breaks = c(0, 40, 80, 120, 180), expand = c(0, 0)) +
  coord_cartesian(xlim = c(1, 14), ylim = c(0, 180), clip = "off") + 
  scale_colour_manual(values = c(
    "DNg27 > w1118"    = "black",
    "w1118 > Kir2.1" = "gray40",
    "DNg27 > Kir2.1"   = "red"
  ),
  labels = c(
    "DNg27 > w1118"  = expression(italic("DNg27 > w"^"1118")),
    "w1118 > Kir2.1" = expression(italic("w"^"1118"*" > Kir2.1")),
    "DNg27 > Kir2.1" = expression(italic("DNg27 > Kir2.1"))
  )) +
  labs(x = "Day", y = "Egg count") +
  theme_classic() +
  theme(
    panel.background  = element_rect(fill = "transparent", color = NA),
    plot.background   = element_rect(fill = "transparent", color = NA),
    panel.grid        = element_blank(),
    axis.line         = element_line(color = "black"),   # both axes
    axis.text.x       = element_text(size = 10, color = "black"),
    axis.title.x      = element_text(size = 10, color = "black"),
    axis.ticks.x      = element_line(color = "black"),
    axis.text.y       = element_text(size = 10, color = "black"),
    axis.title.y      = element_text(size = 10, color = "black"),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.text       = element_text(size = 8, color = "black"),
    legend.key.size   = unit(0.4, "cm"),
    legend.spacing.y  = unit(0.2, "cm"),
    legend.title      = element_blank()
  )

# add statistics
dng27_eggs_stats <- dng27_eggs +
  annotate("segment", x = 12, xend = 14, y = 40, yend = 40, color = "black",
           linewidth = 0.6) +
  annotate("text", x = 13, y = 35, label = "***", size = 4, color = "black",
           hjust = 0.5)

print(dng27_eggs_stats)

if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename_f, ".pdf"))
  ggsave(pdf_path, plot = dng27_eggs_stats, width = 16, height = 8, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
  
  png_path <- file.path(PATH_output, fig_folder_name, paste0(filename_f, ".png"))
  ggsave(png_path, plot = dng27_eggs_stats, width = 16, height = 8, units = "cm",
         bg = "transparent")
  cat("Saved PNG to:", normalizePath(png_path), "\n")
}
#-------------------------------------------------------------------------------