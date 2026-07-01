# -- Figure 3 & S5  ------------------------------------------------------------
# ------------------------------------------------------------------------------
# Figure 3: A-J - input to NSCs
# S5: input to NSCs - bar plot of individual neuron input to multiple NSCs
# ------------------------------------------------------------------------------

# To replicate paper figures, download version from zenodo: 10.5281/zenodo.20795211

# ------------------------------------------------------------------------------
# How files were originally downloaded from Codex before running this script:
# ------------------------------------------------------------------------------
# All csv files were downloaded from Codex on 06/12/2024
# https://codex.flywire.ai/api/download?dataset=fafb

# If files are *not* in the input folder, code below will prompt you to add:

# Files needed:
# 1.) Classification / Hierarchical Annotations: classification.csv
      # last updated 02/24/2026 on codex - to replicate paper figures, download 
      # previous version from zenodo: 10.5281/zenodo.20795211
# 2.) Neurotransmitter Type Predictions: neurotransmitters.csv
      # *Note*: rename the filename to be: neurotransmitters (not neurons)
      # last updated 12/05/2024 on codex - to replicate paper figures, download 
      # previous version from zenodo: 10.5281/zenodo.20795211
# 3.) Consolidated Cell Types: consolidated_cell_types.csv
      # last updated 02/24/2026 on codex - to replicate paper figures, download 
      # previous version from zenodo: 10.5281/zenodo.20795211

# Put files in the 'input' directory: NSC_Connectome/input/

# load packages ----------------------------------------------------------------
library(tidyverse)
library(fafbseg)
library(ggplot2)
library(natverse)
library(ggtext)
library(patchwork)
library(cowplot)
library(ggalluvial)

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
fig_folder_name = "Figure_3"
if(!dir.exists(file.path(PATH_output, fig_folder_name))){
  dir.create(file.path(PATH_output, fig_folder_name))
  cat("Created figure folder at:", normalizePath(file.path(PATH_output, fig_folder_name)), "\n")
}

# ------------------------------------------------------------------------------
# Options for saving plots & csv files
# ------------------------------------------------------------------------------
# When running script for the first time, use settings provided:
write_plots = TRUE           # TRUE - save/replicate figure plots
                             # FALSE - plots not saved outside of R
write_csv   = TRUE           # TRUE - save processed data associated w/figures
                             # FALSE - data not saved outside of R
write_mesh = FALSE           # FALSE - do not generate mesh plots by default (takes awhile to render)

# set colors -------------------------------------------------------------------
m_NSC_unknown = "#BFD739"
l_NSC_unknown = "#009817"
m_NSC_DILP    = "#BD0023"
m_NSC_DH44    = "#ffe200"
m_NSC_DMS     = "#FE7E00"
l_NSC_CRZ     = "#8100FF"
l_NSC_ITP     = "#838383"
l_NSC_DH31    = "#00B4FF"
SEZ_NSC_CAPA  = "#003BBD"
SEZ_NSC_Hugin = "#BD00B0"

syn_col <- c("m_NSC_unknown" = m_NSC_unknown,"l_NSC_unknown" = l_NSC_unknown,
             "m_NSC_DILP" = m_NSC_DILP,"m_NSC_DH44" = m_NSC_DH44,
             "m_NSC_DMS" = m_NSC_DMS,"l_NSC_CRZ" = l_NSC_CRZ,
             "l_NSC_ITP" = l_NSC_ITP,"l_NSC_DH31" = l_NSC_DH31,
             "SEZ_NSC_CAPA" = SEZ_NSC_CAPA,"SEZ_NSC_Hugin" = SEZ_NSC_Hugin)

super_class_colors <- c("others"             = "#5f50a1", 
                        "undefined"          = "white",
                        "optic"              = "#999999",
                        "descending"         = "#c3905f",  
                        "visual_centrifugal" = "#f46d43",
                        "central"            = "#3384b8", 
                        "visual_projection"  = "#c8d684", 
                        "sensory"            = "#b73545", 
                        "ascending"          = "#a9d5a3")

neurotransmitter_colors <- c("ACH"       = "#B9529F",
                             "GLUT"      = "#51B848",
                             "GABA"      = "orange", 
                             "uncertain" = "black",
                             "other"     = "lightgrey")

custom_labels_map <- c(
  SEZ_NSC_Hugin  = "SEZ-NSC<sup>Hugin</sup>",
  SEZ_NSC_CAPA   = "SEZ-NSC<sup>CAPA</sup>",
  l_NSC_DH31     = "l-NSC<sup>DH31</sup>",
  l_NSC_ITP      = "l-NSC<sup>ITP</sup>",
  l_NSC_CRZ      = "l-NSC<sup>CRZ</sup>",
  m_NSC_DMS      = "m-NSC<sup>DMS</sup>",
  m_NSC_DH44     = "m-NSC<sup>DH44</sup>",
  m_NSC_DILP     = "m-NSC<sup>DILP</sup>",
  l_NSC_unknown  = "l-NSC<sup>unknown</sup>",
  m_NSC_unknown  = "m-NSC<sup>unknown</sup>"
)
# load data --------------------------------------------------------------------

NSC <- read_delim(file.path(PATH_input, paste0("NSC_v", v, ".csv")),
                delim = ",",
                escape_double = FALSE,
                col_types = cols(NSC_id = col_character()),
                trim_ws = TRUE)

# Query caveclient w/python to get input to NSCs
if (!paste0("NSC_input_filtered_v",v,".csv")  %in%  input_files) {
  reticulate::source_python("./scripts/Figure_3/Figure_3_preparation.py", 
                            envir = NULL)
}
# Load in file saved from running python script above
# If reticulate not setup or errors, download file from zenodo and read in here
synapses_input <- read_csv(file.path(PATH_input, paste0("tmp/NSC_input_filtered_v",v,".csv")),
                          col_types = cols(pre_pt_supervoxel_id = col_character(),
                                           pre_pt_root_id = col_character(),
                                           post_pt_supervoxel_id = col_character(),
                                           post_pt_root_id = col_character()))

# Classification file
if (!paste0("classification_v", v, ".csv") %in% list.files(PATH_input)) {
  stop("please go to zenodo and download the classification file provided (see 
       note at top of script) and save it in './input'.")
}

classification <- read_delim(file.path(PATH_input, paste0("classification_v", v, ".csv")),
                             delim = ",",
                             escape_double = FALSE,
                             col_types = cols(root_id = col_character(), flow = col_character()),
                             trim_ws = TRUE)

# Update classification
classification[classification$hemibrain_type %in% c("HBeyelet"),]$super_class <- "sensory"

# Neurotransmitter classification
if (!("neurotransmitters.csv") %in% list.files(PATH_input)) {
  stop("please go to zenodo and download the classification file provided (see 
       note at top of script) and save it in './input'.")
}

neurotransmitters <- read_delim(file.path(PATH_input, "neurotransmitters.csv"),
                             delim = ",",
                             escape_double = FALSE,
                             col_types = cols(root_id = col_character()),
                             trim_ws = TRUE)

# Consolidated cell types
if (!paste0("consolidated_cell_types.csv") %in% list.files(PATH_input)) {
  stop("please go to zenodo and download the classification file provided (see 
       note at top of script) and save it in './input'.")
}

consolidated_cell_types <- read_delim(file.path(PATH_input, "consolidated_cell_types.csv"),
                                      delim = ",",
                                      escape_double = FALSE,
                                      col_types = cols(root_id = col_character(),
                                                       primary_type = col_character(),
                                                       `additional_type(s)` = col_character()),
                                      trim_ws = TRUE)

# rename primary_type to cell_type & additional_types(s) to cell_type_additional
consolidated_cell_types <- consolidated_cell_types %>%
  rename(cell_type = primary_type,
         cell_type_additional = `additional_type(s)`)

# ------------------------------------------------------------------------------
# Run python scripts - open .py in R (requires reticulate) or run elsewhere
# otherwise, download files from zenodo and put into: output/neurons_v783/brainmesh/
# ------------------------------------------------------------------------------
# Save brainmesh.obj for 3D plotting -------------------------------------------
# Figure_3_S1_brainmesh.py

# ------------------------------------------------------------------------------
# Figure 3A - Plot input synapses of all NSCs
# ------------------------------------------------------------------------------

synapses_input$pre_x = xyzmatrix(synapses_input$pre_pt_position)[,1]
synapses_input$pre_y = xyzmatrix(synapses_input$pre_pt_position)[,2]
synapses_input$pre_z = xyzmatrix(synapses_input$pre_pt_position)[,3]
synapses_input$post_x = xyzmatrix(synapses_input$post_pt_position)[,1]
synapses_input$post_y = xyzmatrix(synapses_input$post_pt_position)[,2]
synapses_input$post_z = xyzmatrix(synapses_input$post_pt_position)[,3]

if (write_mesh) {
  # Brain mesh w/synaptic coordinates ------------------------------------------
  panel_name = "A1"
  filename_input = paste0(fig_folder_name, panel_name)
  # Open a new 3D plot
  open3d()
  # Set high resolution for the plot
  par3d(windowRect = c(0, 0, 3840, 2160))
  # Iterate over unique NSC names and plot synapses
  for (i in unique(NSC$NSC_name)) {
    # Get the NSC_ids for the current NSC_name
    NSC_ids <- NSC %>% filter(NSC_name == i) %>% pull(NSC_id)
    # Filter the synapses for the current NSC_ids
    synapses_filtered <- synapses_input %>% filter(post_pt_root_id %in% NSC_ids)
    # Plot the pre-synaptic positions
    plot3d(
      synapses_filtered$post_x,
      synapses_filtered$post_y,
      synapses_filtered$post_z,
      col = syn_col[i],
      size = 0.5,
      type = "s",
      add = TRUE
    )
  }
  
  # Add the surface model
  brainmesh <- readOBJ(paste0(PATH_output, "/neurons_v783/brainmesh/brainmesh.obj"))
  
  plot3d(brainmesh,
         add = TRUE,
         alpha = 0.1,
         col = "grey")
  # Adjust the view
  view3d(userMatrix = rotationMatrix(90 * pi / 90, 1, 0, 0), zoom = 0.5)
  
  if (write_plots) {
    # Export as html
    p <- rglwidget(webgl = TRUE,
                   width = 1920,
                   height = 1080)
    html_filename <- file.path(PATH_output,
                               fig_folder_name,
                               paste0(filename_input, ".html"))
    htmltools::save_html(p, html_filename)
    cat("Saved HTML to:", normalizePath(html_filename), "\n")
    
    # Export as png
    png_filename <- file.path(PATH_output,
                              fig_folder_name,
                              paste0(filename_input, ".png"))
    rgl.snapshot(png_filename)
    cat("Saved PNG to:", normalizePath(png_filename), "\n")
  }
  
  # Close the 3D plot
  close3d()
}
# Combined density plot of input synapses along the x-axis by NSC type ---------
panel_name = "A2"
filename_input2 = paste0(fig_folder_name, panel_name)

x_synapse <- ggplot() +
labs(x = "x coordinate", y = "# synapses") +
  theme_classic() +
  xlim(208000, 852000) + # based on an approximate of mesh coordinates in codex
  ylim (0, 3000)         # fixed the ymax limit manually

for (i in unique(NSC$NSC_name)) {
  # Get the NSC_ids for the current NSC_name
  NSC_ids <- NSC %>% filter(NSC_name == i) %>% pull(NSC_id)
  # Filter the synapses for the current NSC_ids
  synapses_filtered <- synapses_input %>% filter(post_pt_root_id %in% NSC_ids)
  # Skip cell types that have 1 or less synapses
  if (nrow(synapses_filtered) <= 1) {
    next
  }
  # Plot the pre-synaptic density
  x_synapse <- x_synapse +
    geom_freqpoly(data = synapses_filtered, aes(x = post_x), linewidth = 1, 
                  bins = 100, color = c(syn_col[i]), na.rm = TRUE)
}
print(x_synapse)

if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename_input2, ".pdf"))
  ggsave(pdf_path, plot = x_synapse, width = 6, height = 1.5)
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}

# Combined density plot of input synapses along the y-axis by NSC type ---------
panel_name = "A3"
filename_input3 = paste0(fig_folder_name, panel_name)

y_synapse <- ggplot() +
  labs(x = "# synapses", y = "y coordinate") +
  theme_classic() +
  xlim (0, 1500) +
  ylim (393600, 81200) # based on an approximate of mesh coordinates in codex

for (i in unique(NSC$NSC_name)) {
  # Get the NSC_ids for the current NSC_name
  NSC_ids <- NSC %>% filter(NSC_name == i) %>% pull(NSC_id)
  # Filter the synapses for the current NSC_ids
  synapses_filtered <- synapses_input %>% filter(post_pt_root_id %in% NSC_ids)
  # Skip cell types that have 1 or less synapses
  if (nrow(synapses_filtered) <= 1) {
    next
  }
  # Plot the pre-synaptic density
  y_synapse <- y_synapse +
    geom_freqpoly(data = synapses_filtered, aes(y = post_y), linewidth = 1, 
                  bins = 50, color = c(syn_col[i]), na.rm=TRUE)
}
print(y_synapse)

if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename_input3, ".pdf"))
  ggsave(pdf_path, plot = y_synapse, width = 2, height = 3.3)
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Organize connectivity data for all figures
# ------------------------------------------------------------------------------

# Summarize synapses data
NSC_input_sum <- synapses_input %>%
  group_by(pre_pt_root_id, post_pt_root_id) %>%
  summarise(n_synapses = n(),
            .groups = "drop")

# Save this data in tmp
if (write_csv){
  csv_path <- file.path(PATH_input, "tmp", paste0("NSC_input_sum_filtered_v", 
                                                  v, ".csv"))
  write.csv(NSC_input_sum, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# Make a large classification table with neurotransmitters and other annotations
classification <- left_join(classification, neurotransmitters, by = "root_id")
# Reclassify neurotransmitters based on stringent cutoff
classification$nt_type_stringent <- classification$nt_type
classification <- classification %>%
  mutate(nt_type_stringent = ifelse(nt_type_score <= 0.62, "uncertain", 
                                    nt_type_stringent))

# Join data for further analysis
classification_join <- classification
colnames(classification_join) <- c("pre_pt_root_id", colnames(classification)[-1])
# Filter and keep >= 5 synapses
NSC_input <- left_join(NSC_input_sum[NSC_input_sum$n_synapses >= 5, ], 
                       classification_join, by = "pre_pt_root_id")

# Convert NA values to a detectable value and then mutate
NSC_input$nt_type_score[is.na(NSC_input$nt_type_score)] <- -1
NSC_input <- NSC_input %>%
  mutate(nt_type_stringent = ifelse(nt_type_score == -1, "uncertain", 
                                    nt_type_stringent))

# Collapse neurotransmitters into 5 categories (GABA, GLUT, ACH, uncertain, other)
NSC_input$nt_type_stringent_corr <- NSC_input$nt_type_stringent
NSC_input <- NSC_input %>%
  mutate(nt_type_stringent_corr = case_when(
    nt_type_stringent_corr %in% c("GABA", "GLUT", "ACH", "uncertain") ~ nt_type_stringent_corr,
    TRUE ~ "other" ))

# Join dfs and rename cols - add info from NSCs (post info such as name)
NSC_join <- NSC
colnames(NSC_join) <- c("name_post", "post_pt_root_id", "hemisphere_post")
NSC_input <- left_join(NSC_input, NSC_join, by = "post_pt_root_id")
colnames(NSC_join) <- c("name_pre", "pre_pt_root_id", "hemisphere_pre")
NSC_input <- left_join(NSC_input, NSC_join, by = "pre_pt_root_id")

# add name text to the id in a new column
NSC_input$post_pt_root_id_name <- paste(NSC_input$post_pt_root_id, 
                                        NSC_input$name_post, sep = "_")

# ------------------------------------------------------------------------------
# **Save data**
# This is the file that is filtered, organized & processed for all following 
# plots & analyses of inputs (Figure 3 & Supplementary Files)
# ------------------------------------------------------------------------------
if (write_csv){
  csv_path <- file.path(PATH_input, paste0("NSC_input_classification_v", v, ".csv"))
  write.csv(NSC_input, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 3B - # input neurons to NSCs by super class
# ------------------------------------------------------------------------------
cell_count <- NSC_input %>%
  group_by(super_class) %>%
  summarise(n_superclass = length(unique(pre_pt_root_id)),
            n_NSC = length(unique(post_pt_root_id))) %>%
  arrange(desc(n_superclass))

# Save data to csv
if (write_csv){
  csv_path <- file.path(PATH_output, fig_folder_name, 
                        paste0("Figure_3B_numNeurons_inputNSC_v", v, ".csv"))
  write.csv(cell_count, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}
# ------------------------------------------------------------------------------

# Summarize NSC input data - total synapses & how many different pre neurons to NSC
NSC_input_total <- NSC_input %>%
  group_by(name_post) %>%
  summarise(n_synapses_total = sum(n_synapses, na.rm = TRUE),
            n_pre_partners_total = length(unique(pre_pt_root_id)))

# Summarize NSC input data - by super class
NSC_input_sum_grouped <- NSC_input %>%
  group_by(super_class, name_post) %>%
  summarise(n_synapses_sum = sum(n_synapses, na.rm = TRUE),
            avrg_synapses = mean(n_synapses, na.rm = TRUE),
            n_pre_partners = length(unique(pre_pt_root_id)),
            n_post_partners = length(unique(post_pt_root_id)),
            .groups = "drop")

# Combine groupings from above^ by NSC (post)
NSC_input_sum_grouped <- left_join(NSC_input_sum_grouped, NSC_input_total, 
                                   by = "name_post")
# Calculate % of input
NSC_input_sum_grouped$perc_of_input <- NSC_input_sum_grouped$n_synapses_sum / 
                                       NSC_input_sum_grouped$n_synapses_total

NSC_input_sum_grouped$name_post <- factor(NSC_input_sum_grouped$name_post,
                                          levels = rev(c("m_NSC_unknown", "l_NSC_unknown", 
                                                         "m_NSC_DILP", "m_NSC_DH44", 
                                                         "m_NSC_DMS", "l_NSC_CRZ",
                                                         "l_NSC_ITP", "l_NSC_DH31", 
                                                         "SEZ_NSC_CAPA", "SEZ_NSC_Hugin")))

NSC_input_sum_grouped$fill <- paste(NSC_input_sum_grouped$super_class, sep = "_")

# Adjust classes for plotting
NSC_input_sum_grouped$super_class_corr <- NSC_input_sum_grouped$super_class
NSC_input_sum_grouped$super_class_corr[NSC_input_sum_grouped$perc_of_input < 0.025] <- "others"
NSC_input_sum_grouped$super_class_corr[is.na(NSC_input_sum_grouped$super_class)] <- "undefined"
NSC_input_sum_grouped$fill_corr <- paste(NSC_input_sum_grouped$super_class_corr)

NSC_input_sum_grouped <- NSC_input_sum_grouped %>%
  group_by(name_post, fill_corr) %>%
  summarise(perc_of_input = sum(perc_of_input),
            .groups = "drop")

NSC_input_sum_grouped$fill_corr <- factor(NSC_input_sum_grouped$fill_corr,
                                          levels = rev(c("others",
                                                         "undefined",
                                                         "optic",
                                                         "descending",
                                                         "visual_centrifugal",
                                                         "central",
                                                         "visual_projection",
                                                         "sensory",
                                                         "ascending")))

# Diagnostic print
# print(NSC_input_sum_grouped)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 3C - Proportion input synapses to NSCs by super class
# ------------------------------------------------------------------------------
panel_name = "C"
filename_input_sc = paste0(fig_folder_name, panel_name)

# no ITP data
NSC_input_sum_grouped$name_post <- droplevels(NSC_input_sum_grouped$name_post)

colored_labels <- setNames(
  paste0("<span style='color:", syn_col[levels(NSC_input_sum_grouped$name_post)], "'>", 
         custom_labels_map[levels(NSC_input_sum_grouped$name_post)], "</span>"),
  levels(NSC_input_sum_grouped$name_post)
)

p <- ggplot(NSC_input_sum_grouped, aes(x = name_post, y = perc_of_input, 
                                       fill = fill_corr)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  facet_grid(scales = "free_x", space = "free_x") +
  scale_fill_manual(values = super_class_colors,
                    guide = guide_legend(nrow = 1)) +
  scale_x_discrete(limits = levels(NSC_input_sum_grouped$name_post), 
                   labels = colored_labels) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = function(x) ifelse(x %in% c(0, 0.5, 1), x, "")
  ) +
  ylab("Proportion of input synapses") +
  xlab("") +
  theme(panel.background = element_rect(fill = NA, color = NA),
        strip.background = element_rect(colour = NA, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        legend.title = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 12, colour = "black"),
        axis.text.x = element_markdown(angle = 270, hjust = 0, vjust = 0.5, colour = "black"),
        axis.text.y = element_text(angle = 0, colour = "black"),
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 8),
        legend.spacing.x = unit(0.1, "cm"))
print(p)

if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename_input_sc, ".pdf"))
  ggsave(pdf_path, plot = p, width = 12, height = 10, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Figure 3D - Top inputs from super classes to NSCs
# ------------------------------------------------------------------------------
panel_name = 'D'
# Ensure output directory exists
nsc_input_folder <- paste0(fig_folder_name, panel_name, "_NSC_inputs")
dir.create(file.path(PATH_output, fig_folder_name, nsc_input_folder), 
           recursive = TRUE, showWarnings = FALSE)
# Ensure output directory exists
nsc_top_input_folder <- paste0(fig_folder_name, panel_name, "_Top_inputs_allNSC")
dir.create(file.path(PATH_output, fig_folder_name, nsc_top_input_folder), 
           recursive = TRUE, showWarnings = FALSE)

filename_top_input = file.path(fig_folder_name, nsc_input_folder,
                               paste0(fig_folder_name, panel_name, 
                                      "_NSC_input_filtered_5_syn_threshold_"))

filename_top10 = file.path(fig_folder_name, nsc_top_input_folder,
                           paste0(fig_folder_name, panel_name, "_input_to_NSC_from_"))


col_to_keep <- c("pre_pt_root_id","cell_type","class","hemibrain_type","name_pre",
               "name_post","hemisphere_post","post_pt_root_id","n_synapses", 
               "nt_type_score", "nt_type", "nt_type_stringent", "nt_type_stringent_corr")

for (i in unique(NSC_input$super_class)) {
  if (!is.na(i)) {
  # Filter the data for the current NSC_ids
  NSC_subset_superclass <- unique(NSC_input[NSC_input$super_class == i,][col_to_keep])
  # Name and save file
  filename <- paste0(filename_top_input, i)
  csv_path <- file.path(PATH_output, paste0(filename, "_v", v, ".csv"))
  write.csv(NSC_subset_superclass, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")  
  
  # Assign values to hemibrain_type from cell_type or class if they are NA
  NSC_subset_superclass$hemibrain_type <- ifelse(is.na(NSC_subset_superclass$hemibrain_type),
                                            ifelse(is.na(NSC_subset_superclass$cell_type),
                                                   NSC_subset_superclass$class,
                                                   NSC_subset_superclass$cell_type),
                                            NSC_subset_superclass$hemibrain_type)
  
  # Filter out rows where hemibrain_type is still NA
  NSC_subset_superclass <- NSC_subset_superclass[!is.na(NSC_subset_superclass$hemibrain_type), ]
  
  top_tmp <- NSC_subset_superclass %>%
    group_by(hemibrain_type) %>%
    summarise(total_syn = sum(n_synapses))
  top10_tmp <- top_tmp[order(-top_tmp$total_syn), ]$hemibrain_type[1:10]
  NSC_subset_superclass_top = NSC_subset_superclass[NSC_subset_superclass$hemibrain_type %in% top10_tmp, ]
  
  # Create a data frame with hemibrain_type and their respective order
  order_df <- data.frame(hemibrain_type = top10_tmp, order = 1:10)
  # Merge the order information into NSC_subset_superclass_top
  NSC_subset_superclass_top <- merge(NSC_subset_superclass_top, order_df, 
                                     by = "hemibrain_type")
  # Save to csv
  filename_top10_full <- paste0(filename_top10, i, "_top_10_v", v, ".csv")
  csv_path_top10 <- file.path(PATH_output, filename_top10_full)
  write.csv(NSC_subset_superclass_top, csv_path_top10)
  cat("Saved CSV to:", normalizePath(csv_path_top10), "\n")
  }
}

# pre_pt_root_ids from each file based on super class 
# (e.g. Figure_3D_input_to_NSC_from_ascending_top_10_v783.csv) are put into
# neuroglancer (https://edit.flywire.ai/) for visualization & screenshots
# neuroglancer links provided as rtf's on zenodo
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 3E - (1) # strong connections to NSCs & (2) # synapses per NSC type
# ------------------------------------------------------------------------------
panel_name = "E"
filename_strong_input = paste0(fig_folder_name, panel_name)

# Get top inputs (interneurons-pre_pt_root_id) to NSC based on n_synapses
top_synp_INT <- NSC_input %>% filter(n_synapses >=50)

# Unique IDs to save
top_synp_INT_unique_IDs <- as.data.frame(unique(top_synp_INT$pre_pt_root_id)) 
colnames(top_synp_INT_unique_IDs) <- "pre_pt_root_id"

# Save as csv
if (write_csv){
  filename <- paste0(filename_strong_input, "_strong_synp_input_IDs")
  csv_path <- file.path(PATH_output, fig_folder_name, 
                        paste0(filename, "_v", v, ".csv"))
  write.csv(top_synp_INT_unique_IDs, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 3E (1-L) - Bar plot for # strong connections to NSCs 
# ------------------------------------------------------------------------------
# Summarize cell counts
top_snyp_INT_sum <- top_synp_INT %>%
  group_by(name_post) %>%
  summarise(cell_count = n())

top_snyp_INT_sum$name_post <- factor(top_snyp_INT_sum$name_post,
                                     levels = rev(c("m_NSC_unknown", "l_NSC_unknown", 
                                                    "m_NSC_DILP", "m_NSC_DH44", 
                                                    "m_NSC_DMS", "l_NSC_CRZ",
                                                    "l_NSC_ITP", "l_NSC_DH31", 
                                                    "SEZ_NSC_CAPA","SEZ_NSC_Hugin")))

# no ITP data
top_snyp_INT_sum$name_post <- droplevels(top_snyp_INT_sum$name_post)

colored_labels <- setNames(
  paste0("<span style='color:", syn_col[levels(top_snyp_INT_sum$name_post)], "'>", 
         custom_labels_map[levels(top_snyp_INT_sum$name_post)], "</span>"),
  levels(top_snyp_INT_sum$name_post)
)

b <- ggplot(top_snyp_INT_sum, aes(x = name_post, y = cell_count, fill = name_post)) + 
  geom_col() + 
  labs(y = "# of strong connections") + 
  scale_fill_manual(values = syn_col) + 
  scale_x_discrete(limits = levels(top_snyp_INT_sum$name_post), 
                   labels = colored_labels) +
  theme_classic() + 
  theme(axis.text.x = element_markdown(angle = 270, hjust = 0, vjust = 0.5, 
                                       color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 8, color = "black"), 
        axis.line = element_line(color = "black"), 
        axis.ticks.x = element_blank(),
        legend.position = "none",
        text = element_text(size = 10, color = "black")) + 
  scale_y_continuous(breaks = seq(0, 20, by = 5), limits = c(0, 20), 
                     expand = c(0, 0)) 
print(b)

if (write_plots){
  # Save to pdf
  pdf_filename <- paste0(filename_strong_input, "1.pdf")
  pdf_path <- file.path(PATH_output, fig_folder_name, pdf_filename)
  ggsave(pdf_path, plot = b, width = 8, height = 5, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 3E (2-R) - Bar plot for # synapses to NSCs (for strong connections)
# ------------------------------------------------------------------------------
top_synp_INT$name_post <- factor(top_synp_INT$name_post,
                                 levels = rev(c("m_NSC_unknown", "l_NSC_unknown", 
                                                "m_NSC_DILP","m_NSC_DH44", 
                                                "m_NSC_DMS", "l_NSC_CRZ",
                                                "l_NSC_ITP", "l_NSC_DH31", 
                                                "SEZ_NSC_CAPA", "SEZ_NSC_Hugin")))
# no ITP data
top_synp_INT$name_post <- droplevels(top_synp_INT$name_post)

colored_labels <- setNames(
  paste0("<span style='color:", syn_col[levels(top_synp_INT$name_post)], "'>", 
         custom_labels_map[levels(top_synp_INT$name_post)], "</span>"),
  levels(top_synp_INT$name_post)
)

g <- ggplot(top_synp_INT, aes(x = name_post, y = n_synapses, fill = name_post)) +
  geom_col() +
  labs(y = "# of synapses") +
  scale_fill_manual(values = syn_col) +
  scale_x_discrete(limits = levels(top_synp_INT$name_post), 
                   labels = colored_labels) +
  theme_classic() +
  theme(axis.text.x = element_markdown(angle = 270, hjust = 0, vjust = 0.5, 
                                       color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 8, color = "black"), 
        axis.line = element_line(color = "black"), 
        axis.ticks.x = element_blank(),
        legend.position = "none",
        text = element_text(size = 10, color = "black")) + 
  scale_y_continuous(breaks = seq(0, 2000, by = 500), limits = c(0, 2000), 
                     expand = c(0,0))
print(g)

if (write_plots){
  # Save to pdf
  pdf_filename <- paste0(filename_strong_input, "2.pdf")
  pdf_path <- file.path(PATH_output, fig_folder_name, pdf_filename)
  ggsave(pdf_path, plot = g, width = 8, height = 5, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 3F - Top input (pre) neurons to NSCs (>= 50 synapses)
# ------------------------------------------------------------------------------
panel_name = "F"
filename_strong_neurons = paste0(fig_folder_name, panel_name)

# Summarize pre by name_post for Figure 3F plots
top_snyp_INT_ordered <- top_synp_INT[order(top_synp_INT$name_post),]

# Save as csv
if (write_csv){
  filename <- paste0(filename_strong_neurons, "_strong_inputs_to_NSCs")
  csv_path <- file.path(PATH_output, fig_folder_name, 
                        paste0(filename, "_v", v, ".csv"))
  write.csv(top_snyp_INT_ordered, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# IDs in plot for neuroglancer
nsc_types_of_interest <- c("SEZ_NSC_CAPA", "l_NSC_DH31", "m_NSC_DH44", "m_NSC_DILP")

for (nsc_type in nsc_types_of_interest) {
  ids <- top_snyp_INT_ordered %>%
    filter(name_post == nsc_type) %>%
    pull(pre_pt_root_id) %>%
    unique() %>%
    as.character()
  
  cat(nsc_type, paste0("(n = ", length(ids), "):\n"), paste(ids, collapse = ", "), "\n\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 3G thru J - Neurons providing inputs to multiple NSC types 
# ------------------------------------------------------------------------------
# Start by getting proportions

# Count the number of unique name_post for each pre_pt_root_id 
pre_pt_post_count <- NSC_input %>%
  group_by(pre_pt_root_id) %>%
  summarise(n_unique_post = n_distinct(name_post))

# Merge with NSC_input
NSC_input_join <- NSC_input
NSC_input_join <- left_join(NSC_input_join, pre_pt_post_count, by = "pre_pt_root_id")

# Identify pre_pt_root_id with more than one type of post synaptic NSC
pre_pt_with_multiple_posts <- pre_pt_post_count %>%
  filter(n_unique_post > 1) %>%
  select(pre_pt_root_id)
# Join back to the original data
NSC_input_multi <- NSC_input_join %>%
  filter(pre_pt_root_id %in% pre_pt_with_multiple_posts$pre_pt_root_id)

# Calculate proportion of synapses
proportion_data <- NSC_input_multi %>%
  group_by(pre_pt_root_id, name_post) %>%
  summarise(total_synapses = sum(n_synapses), .groups = "drop_last") %>%
  mutate(proportion = total_synapses / sum(total_synapses)) %>%
  ungroup()


# ------------------------------------------------------------------------------
# Figure G - Sankey plot - Input to multiple NSCs 
# ------------------------------------------------------------------------------
panel_name = "G"
panel_type = "_prop_input_multipleNSC"
filename_prop_celltype = paste0(fig_folder_name, panel_name, panel_type)
filename_prop = paste0(fig_folder_name, panel_name)

# Add cell type and super class info for pre ids to proportion data
proportion_data_details <- proportion_data %>%
  left_join(classification %>% select(root_id, super_class), 
            by = c("pre_pt_root_id" = "root_id"))

proportion_data_details <- proportion_data_details %>%
  left_join(consolidated_cell_types, by = c("pre_pt_root_id" = "root_id"))

proportion_data_details <- proportion_data_details %>%
  relocate(super_class, cell_type, .after = 1)

# check for unknown cell type and drop those rows (4 ids)
proportion_data_details <- subset(proportion_data_details, !is.na(cell_type))

# group by cell_type & name_post, then sum synapses and then calc proportion 
# (so it is all by cell type)
# Calculate proportion of synapses
prop_cell_type <- proportion_data_details %>%
  group_by(super_class, cell_type, name_post) %>%
  summarise(total_synapses = sum(total_synapses), .groups = "drop_last") %>%
  mutate(proportion = total_synapses / sum(total_synapses)) %>%
  ungroup()

# prepare for plotting
super_class_light <- c("ascending"  = "#dff0dc",
                       "central"    = "#cce4f3",
                       "descending" = "#f0dfd0")

prop_cell_type <- prop_cell_type %>%
  filter(!is.na(cell_type)) %>%
  mutate(super_class = factor(super_class, levels = c("ascending", "central", "descending"))) %>%
  arrange(super_class, cell_type) %>%
  mutate(
    cell_type = factor(cell_type, levels = unique(cell_type)),
    name_post = factor(as.character(name_post), 
                       levels = rev(names(syn_col))[rev(names(syn_col)) %in% unique(as.character(name_post))])
  )

if (write_csv){
  csv_path <- file.path(PATH_output, fig_folder_name, 
                                paste0(filename_prop_celltype, "_v", v, ".csv"))
  write.csv(prop_cell_type, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# pre-neuron background colors
cell_type_bg <- super_class_light[as.character(prop_cell_type$super_class)]
names(cell_type_bg) <- as.character(prop_cell_type$cell_type)
cell_type_bg <- cell_type_bg[!duplicated(names(cell_type_bg))]

# calculate geom_rect positions for pre stratum coloring
total_y <- sum(prop_cell_type %>%
                 group_by(cell_type) %>%
                 summarise(total = sum(proportion)) %>%
                 pull(total))

cell_type_positions <- prop_cell_type %>%
  group_by(cell_type) %>%
  summarise(total = sum(proportion), super_class = first(super_class)) %>%
  arrange(super_class, cell_type) %>%
  mutate(
    ymax         = cumsum(total),
    ymin         = lag(ymax, default = 0),
    ymin_flipped = total_y - ymax,
    ymax_flipped = total_y - ymin,
    fill_color   = super_class_light[as.character(super_class)]
  )

width <- 0.08  # control all widths from one place

pre_xmin <- 1 - width/2
pre_xmax <- 1 + width/2

s_plot <- ggplot(prop_cell_type,
       aes(axis1 = cell_type, axis2 = name_post, y = proportion)) +
  geom_alluvium(aes(fill = name_post), width = width) +
  geom_stratum(fill = NA, width = width, color = "black", size = 0.2) +
  geom_rect(
    data = cell_type_positions,
    aes(xmin = pre_xmin, xmax = pre_xmax,
        ymin = ymin_flipped, ymax = ymax_flipped),
    fill  = cell_type_positions$fill_color,
    color = "black", size = 0.2,
    inherit.aes = FALSE
  ) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 1.5, color = "black") +
  scale_fill_manual(
    values = syn_col,
    breaks = rev(names(syn_col)),
    labels = custom_labels_map[rev(names(syn_col))],
    name   = "NSC class"
  ) +
  scale_x_discrete(limits = c("Pre", "Post"), expand = c(0.05, 0.05)) +
  theme_minimal() +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid   = element_blank(),
    legend.text  = element_markdown(size = 8),
    axis.text.x  = element_text(color = "black", size = 10, margin = margin(t = -10))
  ) +
  labs(title = "Cell types with input to multiple NSCs")
print(s_plot)

if (write_plots){
  # Save to pdf
  pdf_filename <- paste0(filename_prop,".pdf")
  pdf_path <- file.path(PATH_output, fig_folder_name, pdf_filename)
  ggsave(pdf_path, plot = s_plot, width = 10, height = 8)
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figures H,I,J - Neurons providing inputs to 2(H), 3(I), or 4(J) NSC types 
# ------------------------------------------------------------------------------
panel_name = "HIJ"
panel_type = "_NSC_input_multi_v"
panel_type_sum = "_NSC_input_multi_summary_v"
filename_multi = paste0(fig_folder_name, panel_name, panel_type)
filename_multi_summary = paste0(fig_folder_name, panel_name, panel_type_sum)

# IDs were input into neuroglancer for plotting
if (write_csv){
  csv_path <- file.path(PATH_output, fig_folder_name, paste0(filename_multi, v, ".csv"))
  write.csv(NSC_input_join, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# # of pre neurons to multiple NSCs
pre_pt_post_distribution <- pre_pt_post_count %>%
  group_by(n_unique_post) %>%
  summarise(count_pre_pt_root_id = n())

if (write_csv){
  csv_path <- file.path(PATH_output, fig_folder_name, paste0(filename_multi_summary, v, ".csv"))
  write.csv(pre_pt_post_distribution, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# Loop through each value of n_unique_post (1, 2, 3, 4)
for (n in sort(unique(NSC_input_join$n_unique_post))) {
  
  subset_data <- NSC_input_join %>% filter(n_unique_post == n)
  
  unique_ids <- subset_data %>% pull(pre_pt_root_id) %>% unique()
  unique_name_post <- subset_data %>% pull(name_post) %>% unique()
  
  cat("n_unique_post =", n, "\n")
  cat("  Unique pre_pt_root_id count:", length(unique_ids), "\n")
  cat("  Unique pre_pt_root_id:", paste(unique_ids, collapse = ", "), "\n")
}
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Figure 3 S5 - Prop of input - bar plot - input to multiple NSC types 
# ------------------------------------------------------------------------------
panel_name = "_S5"
filename_prop_NSC = paste0(fig_folder_name, panel_name)

# Order the plotting of NSC proportions
proportion_data$name_post <- factor(proportion_data$name_post,
                                    levels = rev(c("m_NSC_unknown","l_NSC_unknown",
                                                   "m_NSC_DILP","m_NSC_DH44",
                                                   "m_NSC_DMS","l_NSC_CRZ",
                                                   "l_NSC_ITP","l_NSC_DH31",
                                                   "SEZ_NSC_CAPA","SEZ_NSC_Hugin")))

# Drop unused levels (l_NSC_unknown, l_NSC_ITP)
proportion_data$name_post <- droplevels(proportion_data$name_post)

# Determine the ordering of pre_pt_root_id based on name_post
# Order cells such that those providing inputs to similar NSC are grouped together
ordered_pre_pt_root_id <- proportion_data %>%
  arrange(name_post, pre_pt_root_id) %>%
  pull(pre_pt_root_id) %>%
  unique()

proportion_data <- proportion_data %>%
  mutate(pre_pt_root_id = factor(pre_pt_root_id, levels = ordered_pre_pt_root_id))

# Stacked bar chart
pplot <- ggplot(proportion_data, aes(x = pre_pt_root_id, y = proportion, fill = name_post)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1, alpha = 0.7) +
  facet_grid(scales = "free_x", space = "free_x") +
  scale_fill_manual(values = syn_col, guide = guide_legend(nrow = 1)) +
  scale_y_continuous(
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1"),
    expand = expansion(add = c(0, 0.01))
  ) +
  coord_cartesian(ylim = c(0, 1.01), clip = "on") +
  labs(x = "", y = "Proportion of Inputs", fill = "Post Synaptic NSC") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = NA, color = NA),
        strip.background = element_rect(color = NA, fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 10, color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 0, size = 12, color = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black", linewidth = 0.5),
        axis.ticks.length.y = unit(2, "mm"),
        axis.line.y = element_line(color = "black", linewidth = 0.5,
                                   arrow = NULL, lineend = "square"),
        axis.line.x.bottom = element_line(color = "black", linewidth = 0.5))
print(pplot)

if (write_plots){
  # Save to pdf
  pdf_filename <- paste0(filename_prop_NSC, ".pdf")
  pdf_path <- file.path(PATH_output, fig_folder_name, pdf_filename)
  ggsave(pdf_path, plot = pplot, width = 26, height = 8, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}
# ------------------------------------------------------------------------------

