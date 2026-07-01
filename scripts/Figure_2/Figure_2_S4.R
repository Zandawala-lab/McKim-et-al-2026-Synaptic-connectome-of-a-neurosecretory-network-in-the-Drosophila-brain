# -- Figure 2 S4 ---------------------------------------------------------------
# ------------------------------------------------------------------------------
#  Code for plotting NSC morphological characteristics 
# ------------------------------------------------------------------------------

# load packages ----------------------------------------------------------------
library(tidyverse)
library(patchwork)

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

input_files = list.files(path = PATH_input, full.names = FALSE, recursive = FALSE)
# Check if 'tmp' folder exists inside PATH_input, if not, make it
tmp_path <- file.path(PATH_input, "tmp")
if(!dir.exists(tmp_path)){
  dir.create(tmp_path)
  cat("Created tmp folder at:", normalizePath(tmp_path), "\n")
}
input_files_tmp = list.files(path = paste0(PATH_input,"/tmp/"),
                             full.names = FALSE, recursive = FALSE)
input_files = c(input_files,input_files_tmp)

v = read_delim(file.path(PATH_input,"version.csv"),
               col_types  =  cols(version  =  col_character()),delim  =  ";")
v = v$version[1]

# Folder name for saving figs + check if exists
fig_folder_name = "Figure_2"
if(!dir.exists(file.path(PATH_output, fig_folder_name))){
  dir.create(file.path(PATH_output, fig_folder_name))
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Options for saving plots & csv files
# ------------------------------------------------------------------------------

# When running script for the first time, set to TRUE:
write_plots = TRUE           # TRUE - save/replicate figure plots
                             # FALSE - plots not saved outside of R
# set colors -------------------------------------------------------------------
m_NSC_unknown = "#BFD739"
l_NSC_unknown = "#009817"
m_NSC_DILP = "#BD0023"
m_NSC_DH44 = "#ffe200"
m_NSC_DMS = "#FE7E00"
l_NSC_CRZ = "#8100FF"
l_NSC_ITP = "#838383"
l_NSC_DH31 = "#00B4FF"
SEZ_NSC_CAPA = "#003BBD"
SEZ_NSC_Hugin = "#BD00B0"

syn_col <- c("m_NSC_unknown" = m_NSC_unknown,"l_NSC_unknown" = l_NSC_unknown,
             "m_NSC_DILP" = m_NSC_DILP,"m_NSC_DH44" = m_NSC_DH44,
             "m_NSC_DMS" = m_NSC_DMS,"l_NSC_CRZ" = l_NSC_CRZ,
             "l_NSC_ITP" = l_NSC_ITP,"l_NSC_DH31" = l_NSC_DH31,
             "SEZ_NSC_CAPA" = SEZ_NSC_CAPA,"SEZ_NSC_Hugin" = SEZ_NSC_Hugin)

# load data --------------------------------------------------------------------
NSC <- read_delim(file.path(PATH_input,paste0("NSC_v",v,".csv")),
                 col_types = cols(NSC_id = col_character()),delim = ",")

# Check for nuclei volume file
if (!paste0("nuclei_volume_v",v,".csv")  %in%  input_files) {
  reticulate::source_python("./scripts/Figure_2/Figure_2_S4_preparation.py")
}
nuclei <- read_delim(file.path(PATH_input, "tmp", paste0("nuclei_volume_v", v, ".csv")),
                     col_types = cols(pt_root_id = col_character()), delim = ",")

# Check for cell stats file
if (!paste0("cell_stats_v", v, ".csv") %in% list.files(PATH_input)) {
  stop("please go to https://codex.flywire.ai/api/download and download the cell size measurements file for the current version and save it in './input'.")
}
cell_stats <- read_delim(file.path(PATH_input, paste0("cell_stats_v", v, ".csv")),
                         col_types = cols(root_id = col_character()), delim = ",")

# prep for plotting ------------------------------------------------------------
colnames(NSC)[2] <- 'root_id'
colnames(nuclei)[8] <- 'root_id'

NSC <- left_join(NSC, cell_stats, by = "root_id")
NSC <- left_join(NSC, nuclei, by = "root_id")

# Filter out nuclei fragments - 2 DH31 NSC & 1 Hugin NSC - don't have nuclei data
NSC <- subset(NSC, NSC$volume > 1)

NSC$length_um <- NSC$length_nm / 1000
NSC$area_um <- NSC$area_nm / 1000000
NSC$size_um <- NSC$size_nm / 1000000000

NSC$NSC_name <- factor(NSC$NSC_name,
                       levels = rev(c("m_NSC_unknown", "l_NSC_unknown", "m_NSC_DILP",
                                      "m_NSC_DH44", "m_NSC_DMS", "l_NSC_CRZ",
                                      "l_NSC_ITP", "l_NSC_DH31", "SEZ_NSC_CAPA",
                                      "SEZ_NSC_Hugin")))

# ------------------------------------------------------------------------------
# Figure 2 S2A - FAFB - cable length
# ------------------------------------------------------------------------------

length_plot <- ggplot(NSC, aes(x = NSC_name, y = length_um, fill = NSC_name)) +
  geom_violin(trim = FALSE, color = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.5) +
  facet_wrap(~ NSC_name, scales = "fixed") +
  facet_grid(scales = "free_x", space = "free_x") +
  scale_fill_manual(values = syn_col) +
  ylab (~ "Length (µm)") +
  xlab("") +
  theme(panel.background = element_rect(fill = NA, color = NA),
        strip.background = element_rect(colour = NA, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.position = "none",
        text = element_text(size = 10, colour = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 0),
        axis.ticks.x = element_blank())
print(length_plot)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 2 S2B - FAFB - surface area
# ------------------------------------------------------------------------------

area_plot <- ggplot(NSC, aes(x = NSC_name, y = area_um, fill = NSC_name)) +
  geom_violin(trim = FALSE, color = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.5) +
  facet_wrap(~ NSC_name, scales = "fixed") + 
  facet_grid(scales = "free_x", space = "free_x") +
  scale_fill_manual(values = syn_col) +
  ylab("Area (µm²)") +
  xlab("") +
  theme(panel.background = element_rect(fill = NA, color = NA),
        strip.background = element_rect(colour = NA, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.position = "none",
        text = element_text(size = 10, colour = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 0),
        axis.ticks.x = element_blank())
print(area_plot)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 2 S2C - FAFB - cell volume
# ------------------------------------------------------------------------------

size_plot <- ggplot(NSC, aes(x = NSC_name, y = size_um, fill = NSC_name)) +
  geom_violin(trim = FALSE, color = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.5) +
  facet_wrap(~ NSC_name, scales = "fixed") +  
  facet_grid(scales = "free_x", space = "free_x") +
  scale_fill_manual(values = syn_col) +
  ylab("Size (µm³)") +
  xlab("") +
  theme(panel.background = element_rect(fill = NA, color = NA),
        strip.background = element_rect(colour = NA, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.position = "none",
        text = element_text(size = 10, colour = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 0),
        axis.ticks.x = element_blank())
print(size_plot)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 2 S2D - FAFB - nuclei volume
# ------------------------------------------------------------------------------

nuclei_plot <- ggplot(NSC, aes(x = NSC_name, y = volume, fill = NSC_name)) +
  geom_violin(trim = FALSE, color = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.5) +
  facet_wrap(~ NSC_name, scales = "fixed") +  
  facet_grid(scales = "free_x", space = "free_x") +
  scale_fill_manual(values = syn_col) +
  ylab("Nuclei (µm³)") +
  xlab("") +
  theme(panel.background = element_rect(fill = NA, color = NA),
        strip.background = element_rect(colour = NA, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.position = "none",
        text = element_text(size = 10, colour = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 0),
        axis.ticks.x = element_blank())
print(nuclei_plot)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 2 S2E - FAFB - PCA analysis
# ------------------------------------------------------------------------------

# PCA analysis of NSC
pca_data <- NSC[, c("area_um", "length_um", "size_um", "volume")]
# Standardize the data
pca_data_scaled <- scale(pca_data)
# Perform PCA
pca_result <- prcomp(pca_data_scaled, center = TRUE, scale. = TRUE)
# View summary of PCA result
summary(pca_result)

# df with PCA results & plot
pca_df <- as.data.frame(pca_result$x)
pca_df$NSC_name <- NSC$NSC_name
pca_df$root_id <- NSC$root_id

# Plot first two principal components
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = NSC_name)) +
  geom_point() +
  scale_color_manual(values = syn_col) +
  labs(title = "", x = "PC 1", y = "PC 2") +
  theme(panel.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.position = "right")
print(pca_plot)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Patchwork - combine panels A-E
# ------------------------------------------------------------------------------

combined_plot <- (length_plot + area_plot) / 
  (size_plot + nuclei_plot) / 
  (pca_plot + plot_spacer()) +
  plot_layout(guides = "collect")
print(combined_plot)

if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, "/Figure_2_S4.pdf")
  ggsave(pdf_path, plot = combined_plot, width = 16, height = 20, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
  
}
# ------------------------------------------------------------------------------

