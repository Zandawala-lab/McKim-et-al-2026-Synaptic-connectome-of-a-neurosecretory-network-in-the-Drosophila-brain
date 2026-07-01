# -- Figure 7 & 7 S1 -----------------------------------------------------------
# ------------------------------------------------------------------------------
# Code for plotting synaptic outputs of NSCs
# ------------------------------------------------------------------------------

# load packages ----------------------------------------------------------------
library(tidyverse)
library(fafbseg)
library(ggplot2)
library(natverse)
library(ggtext)
library(coconatfly)
library(malecns)
library(neuprintr)

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
fig_folder_name = "Figure_7"
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

# Query caveclient w/python to get NSC output
if (!paste0("NSC_output_filtered_v",v,".csv")  %in%  input_files) {
  reticulate::source_python("./scripts/Figure_7/Figure_7_preparation.py")
}
# Load in file saved from running python script above
# If reticulate not setup or errors, download file from zenodo and read in here
synapses_output =  read_csv(file.path(PATH_input, paste0("tmp/NSC_output_filtered_v", 
                                                         v,".csv")),
                            col_types =  cols(pre_pt_supervoxel_id = col_character(),
                                              pre_pt_root_id = col_character(),
                                              post_pt_supervoxel_id = col_character(),
                                              post_pt_root_id = col_character()))


# Classification file
if (!paste0("classification_v", v, ".csv") %in% list.files(PATH_input)) {
  stop("please go to zenodo and download the classification file provided (see 
       note at top of script) and save it in './input'.")
}

classification <- read_delim(file.path(PATH_input, paste0("classification_v",
                                                          v, ".csv")),
                             delim = ",",
                             escape_double = FALSE,
                             col_types = cols(root_id = col_character(),
                                              flow = col_character()),
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
# ------------------------------------------------------------------------------
# Run python scripts - open .py in R (requires reticulate) or run elsewhere
# otherwise, download files from zenodo and put into: output/neurons_v783/brainmesh/
# ------------------------------------------------------------------------------
# Save brainmesh.obj for 3D plotting -------------------------------------------
# Figure_3_S1_brainmesh.py

# ------------------------------------------------------------------------------
# Figure 7A - Plot output synapses of all NSCs
# ------------------------------------------------------------------------------
synapses_output$pre_x = xyzmatrix(synapses_output$pre_pt_position)[,1]
synapses_output$pre_y = xyzmatrix(synapses_output$pre_pt_position)[,2]
synapses_output$pre_z = xyzmatrix(synapses_output$pre_pt_position)[,3]
synapses_output$post_x = xyzmatrix(synapses_output$post_pt_position)[,1]
synapses_output$post_y = xyzmatrix(synapses_output$post_pt_position)[,2]
synapses_output$post_z = xyzmatrix(synapses_output$post_pt_position)[,3]

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
    synapses_filtered <- synapses_output %>% filter(pre_pt_root_id %in% NSC_ids)
    # Plot the pre-synaptic positions
    plot3d(
      synapses_filtered$pre_x,
      synapses_filtered$pre_y,
      synapses_filtered$pre_z,
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

# Combined density plot of output synapses along the x-axis by NSC type --------
panel_name = "A2"
filename_input2 = paste0(fig_folder_name, panel_name)

x_synapse <- ggplot() +
  labs(x = "x coordinate", y = "# synapses") +
  theme_classic() +
  xlim(208000, 852000) + # based on an approximate of mesh coordinates in codex
  ylim (0, 400)          # fixed the ymax limit manually
for (i in unique(NSC$NSC_name)) {
  # Get the NSC_ids for the current NSC_name
  NSC_ids <- NSC %>% filter(NSC_name == i) %>% pull(NSC_id)
  # Filter the synapses for the current NSC_ids
  synapses_filtered <- synapses_output %>% filter(pre_pt_root_id %in% NSC_ids)
  # Skip cell types that have 1 or less synapses
  if (nrow(synapses_filtered) <= 1) {
    next
  }
  # Plot the pre-synaptic density
  x_synapse <- x_synapse +
    geom_freqpoly(data = synapses_filtered, aes(x = pre_x), linewidth = 1, 
                  bins = 75, color = c(syn_col[i]), na.rm = TRUE)
}

print(x_synapse)


if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename_input2,
                                                             ".pdf"))
  ggsave(pdf_path, plot = x_synapse, width = 6, height = 1.5)
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}


# Combined density plot of output synapses along the y-axis by NSC type --------
panel_name = "A3"
filename_input3 = paste0(fig_folder_name, panel_name)

y_synapse <- ggplot() +
  labs(x = "# synapses", y = "y coordinate") +
  theme_classic() +
  xlim (0, 320) +
  ylim (393600, 81200) # based on an approximate of mesh coordinates in codex

for (i in unique(NSC$NSC_name)) {
  # Get the NSC_ids for the current NSC_name
  NSC_ids <- NSC %>% filter(NSC_name == i) %>% pull(NSC_id)
  # Filter the synapses for the current NSC_ids
  synapses_filtered <- synapses_output %>% filter(pre_pt_root_id %in% NSC_ids)
  # Skip cell types that have 1 or less synapses
  if (nrow(synapses_filtered) <= 1) {
    next
  }
  # Plot the pre-synaptic density
  y_synapse <- y_synapse +
    geom_freqpoly(data = synapses_filtered, aes(y = pre_y), linewidth = 1, 
                  bins = 75, color = c(syn_col[i]), na.rm=TRUE)
}

print(y_synapse)

if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename_input3,
                                                             ".pdf"))
  ggsave(pdf_path, plot = y_synapse, width = 2, height = 3.4)
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Organize connectivity data for all figures
# ------------------------------------------------------------------------------

# Summarize synapses data
NSC_output_sum <- synapses_output %>%
  group_by(pre_pt_root_id, post_pt_root_id) %>%
  summarise(n_synapses = n(),
            .groups = "drop")

# Save this data in tmp
if (write_csv){
  csv_path <- file.path(PATH_input, "tmp", paste0("NSC_output_sum_filtered_v",
                                                  v, ".csv"))
  write.csv(NSC_output_sum, csv_path)
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
colnames(classification_join) <- c("post_pt_root_id", colnames(classification)[-1])
# Filter and keep >= 5 synapses
NSC_output <- left_join(NSC_output_sum[NSC_output_sum$n_synapses >= 5, ], 
                        classification_join, by = "post_pt_root_id")

# Convert NA values to a detectable value and then mutate
NSC_output$nt_type_score[is.na(NSC_output$nt_type_score)] <- -1
NSC_output <- NSC_output %>%
  mutate(nt_type_stringent = ifelse(nt_type_score == -1, "uncertain",
                                    nt_type_stringent))

# Collapse neurotransmitters into 5 categories (GABA, GLUT, ACH, uncertain, other)
NSC_output$nt_type_stringent_corr <- NSC_output$nt_type_stringent
NSC_output <- NSC_output %>%
  mutate(nt_type_stringent_corr = case_when(
    nt_type_stringent_corr %in% c("GABA", "GLUT", "ACH", "uncertain") ~ nt_type_stringent_corr,
    TRUE ~ "other" ))

# Join dfs and rename cols - add info from NSCs
NSC_join <- NSC
NSC_join$NSC_id_old <- NULL
colnames(NSC_join) <- c("name_pre", "pre_pt_root_id", "hemisphere_pre")
NSC_output <- left_join(NSC_output, NSC_join, by = "pre_pt_root_id")
colnames(NSC_join) <- c("name_post", "post_pt_root_id", "hemisphere_post")
NSC_output <- left_join(NSC_output, NSC_join, by = "post_pt_root_id")

# Summarize NSC output data
NSC_output_total <- NSC_output %>%
  group_by(name_pre) %>%
  summarise(n_synapses_total = sum(n_synapses, na.rm = TRUE),
            n_pre_partners_total = length(unique(pre_pt_root_id)))

# add name text to the id in a new column
NSC_output$pre_pt_root_id_name <- paste(NSC_output$pre_pt_root_id,
                                        NSC_output$name_pre, sep = "_")

# ------------------------------------------------------------------------------
# **Save data**
# This is the file that is filtered, organized & processed for all following 
# plots & analyses of output (Figure 7 & Supplementary Files)
# ------------------------------------------------------------------------------
if (write_csv){
  csv_path <- file.path(PATH_input, paste0("NSC_output_classification_v", v, ".csv"))
  write.csv(NSC_output, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Figure 7B - # output neurons to NSCs by super class
# ------------------------------------------------------------------------------

l_NSC_CRZ <- subset(NSC_output, NSC_output$name_pre == "l_NSC_CRZ")
l_NSC_unknown <- subset(NSC_output, NSC_output$name_pre == "l_NSC_unknown")

cell_count_CRZ = l_NSC_CRZ%>%
  group_by(super_class)%>%
  summarise(n_superclass = length(unique(post_pt_root_id)),
            NSC_count = length(unique(pre_pt_root_id)))
cell_count_CRZ$NSC_type <- "l_NSC_CRZ"

cell_count_unknown = l_NSC_unknown%>%
  group_by(super_class)%>%
  summarise(n_superclass = length(unique(post_pt_root_id)),
            NSC_count = length(unique(pre_pt_root_id)))
cell_count_unknown$NSC_type <- "l_NSC_unknown"

NSC_output_superclass <- rbind(cell_count_CRZ, cell_count_unknown)

# Save data to csv
if (write_csv){
  csv_path <- file.path(PATH_output, fig_folder_name, 
                        paste0("Figure_7B_numNeurons_outputNSC_v", v, ".csv"))
  write.csv(NSC_output_superclass, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

#-------------------------------------------------------------------------------

# Summarize NSC output data - total synapses & how many different partners -----
NSC_output_sum_grouped <- NSC_output %>%
  group_by(super_class, name_pre) %>%
  summarise(n_synapses_sum = sum(n_synapses, na.rm = TRUE),
            avrg_synapses = mean(n_synapses, na.rm = TRUE),
            n_pre_partners = length(unique(pre_pt_root_id)),
            n_post_partners = length(unique(post_pt_root_id)),
            .groups = "drop")

# Combine groupings from above^ by NSC
NSC_output_sum_grouped <- left_join(NSC_output_sum_grouped, NSC_output_total, 
                                    by = "name_pre")
# Calculate % of output
NSC_output_sum_grouped$perc_of_output <- NSC_output_sum_grouped$n_synapses_sum / 
                                         NSC_output_sum_grouped$n_synapses_total

NSC_output_sum_grouped$name_pre <- factor(NSC_output_sum_grouped$name_pre,
                                          levels = rev(c("m_NSC_unknown", "l_NSC_unknown", 
                                                         "m_NSC_DILP", "m_NSC_DH44", 
                                                         "m_NSC_DMS", "l_NSC_CRZ",
                                                         "l_NSC_ITP", "l_NSC_DH31", 
                                                         "SEZ_NSC_CAPA", "SEZ_NSC_Hugin")))

NSC_output_sum_grouped$fill <- paste(NSC_output_sum_grouped$super_class, sep = "_")

# Adjust classes for plotting
NSC_output_sum_grouped$super_class_corr <- NSC_output_sum_grouped$super_class
NSC_output_sum_grouped$super_class_corr[NSC_output_sum_grouped$perc_of_output < 0.025] <- "others"
NSC_output_sum_grouped$super_class_corr[is.na(NSC_output_sum_grouped$super_class)] <- "undefined"
NSC_output_sum_grouped$fill_corr <- paste(NSC_output_sum_grouped$super_class_corr)

NSC_output_sum_grouped <- NSC_output_sum_grouped %>%
  group_by(name_pre, fill_corr) %>%
  summarise(perc_of_output = sum(perc_of_output),
            .groups = "drop")

NSC_output_sum_grouped$fill_corr <- factor(NSC_output_sum_grouped$fill_corr,
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
# print(NSC_output_sum_grouped)
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Figure 7C - Proportion output synapses from NSCs by super class
# ------------------------------------------------------------------------------
panel_name = "C"
filename_output_sc = paste0(fig_folder_name, panel_name)

NSC_output_sum_grouped$name_pre <- droplevels(NSC_output_sum_grouped$name_pre)

colored_labels <- setNames(
  paste0("<span style='color:", syn_col[levels(NSC_output_sum_grouped$name_pre)], "'>", 
         custom_labels_map[levels(NSC_output_sum_grouped$name_pre)], "</span>"),
  levels(NSC_output_sum_grouped$name_pre)
)

p <- ggplot(NSC_output_sum_grouped, aes(x = name_pre, y = perc_of_output,
                                        fill = fill_corr)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  facet_grid(scales = "free_x", space = "free_x") +
  scale_fill_manual(values = super_class_colors,
                    guide = guide_legend(nrow = 1)) +
  scale_x_discrete(limits = levels(NSC_output_sum_grouped$name_pre), 
                   labels = colored_labels) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = function(x) ifelse(x %in% c(0, 0.5, 1), x, "")
  ) +
  ylab("Proportion of output synapses") +
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
        axis.text.x = element_markdown(angle = 270, hjust = 0, vjust = 0.5,
                                       colour = "black"),
        axis.text.y = element_text(angle = 0, colour = "black"),
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 8),
        legend.spacing.x = unit(0.1, "cm"))

print(p)

if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename_output_sc,
                                                             ".pdf"))
  ggsave(pdf_path, plot = p, width = 6, height = 10, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 7D & E - Output IDs from NSCs by super class
# ------------------------------------------------------------------------------
# IDs from file put into neuroglancer for visualization

col_to_keep <- c("pre_pt_root_id","cell_type","class","hemibrain_type","name_pre",
               "name_post","hemisphere_post","post_pt_root_id","n_synapses", 
               "nt_type_score", "nt_type", "nt_type_stringent", "nt_type_stringent_corr")

for (i in unique(NSC_output$super_class)) {
  if (!is.na(i)) {
    # Filter the data for the current NSC_ids
    NSC_subset_superclass <- unique(NSC_output[NSC_output$super_class == i,][col_to_keep])
    # Name and save file
    filename <- paste(fig_folder_name, "/Figure_7D_NSC_output_filtered_5_syn_threshold_",
                      i, sep = "")
    outpath <- file.path(PATH_output, paste0(filename,"_v", v,".csv"))
    write.csv(NSC_subset_superclass, outpath)
    cat("Saved CSV to:", normalizePath(outpath), "\n")
  }
}
central_data <- unique(NSC_output[NSC_output$super_class == "central",][col_to_keep])
cat("central - count:", nrow(central_data), "\n")
cat(paste(central_data$post_pt_root_id, collapse = ", "), "\n")

# post_pt_root_ids based on super class are put into
# neuroglancer (https://edit.flywire.ai/) for visualization & screenshots
# neuroglancer links provided as rtf's on zenodo

# ------------------------------------------------------------------------------
# Figure 7E - maleCNS DNg27 neurons
# ------------------------------------------------------------------------------
## Set your Neuprint token: https://natverse.org/neuprintr/#authentication
if (write_mesh) {
  panel_name = "E"
  filename_mCNS = paste0(fig_folder_name, panel_name)
  
  dng27_mcns <- c("11475", "11653")
  
  # choose dataset
  neuprint_dataset <- "male-cns:v0.9"
  
  open3d()
  # Set high resolution for the plot
  par3d(windowRect = c(0, 0, 2560, 1440))
  # Download meshes
  dng27_mcns_meshes <- read_mcns_meshes(dng27_mcns)
  
  # plot meshes
  plot3d(
    dng27_mcns_meshes,
    col = '#c3905f',
    add = TRUE,
    WithNodes = F
  )
  # Add the surface model
  plot3d(
    malecns.surf,
    alpha = 0.1,
    add = TRUE,
    col = "lightgrey"
  )
  plot3d(
    malecnsvnc.surf,
    alpha = 0.1,
    add = TRUE,
    col = "lightgrey"
  )
  # Adjust the view
  view3d(
    fov = 0,
    userMatrix = rotationMatrix(pi, 1, 0, 0) %*%
      rotationMatrix(-pi / 3, 1, 0, 0),
    zoom = 0.8
  )
  
  if (write_plots) {
    # Export as png
    png_filename <- file.path(PATH_output,
                              fig_folder_name,
                              paste0(filename_mCNS, ".png"))
    rgl.snapshot(png_filename)
    cat("Saved PNG to:", normalizePath(png_filename), "\n")
  }
  # Close the 3D plot
  close3d()
}
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Figure 7E & F - CRZ
# ------------------------------------------------------------------------------
panel_name = "E_F"
panel_type = "_CRZ_to_DNg27"
filename_crz = paste0(fig_folder_name, panel_name, panel_type)

crz <- NSC_output[NSC_output$name_pre == "l_NSC_CRZ", ]

crz_summary <- crz %>%
  group_by(pre_pt_root_id, post_pt_root_id, cell_type, side, hemisphere_pre) %>%
  summarise(synapse_sum = sum(n_synapses),
            .groups = "drop")

sorted_crz <- crz_summary %>% filter(!is.na(cell_type)) %>% 
              arrange(desc(pre_pt_root_id), desc(post_pt_root_id), side)

# IDs for E
dng_ids <- unique(sorted_crz$post_pt_root_id)
cat("DNg27 - count:", length(dng_ids), "\n")
cat(paste(dng_ids, collapse = ", "), "\n")

crz_ids <- unique(sorted_crz$pre_pt_root_id)
cat("CRZ - count:", length(crz_ids), "\n")
cat(paste(crz_ids, collapse = ", "), "\n")

if (write_csv){
  csv_path <- file.path(PATH_output, fig_folder_name, 
                        paste0(filename_crz, ".csv"))
  write.csv(sorted_crz, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 7G - Neurotransmitters of outputs
# ------------------------------------------------------------------------------
panel_name = "G_"
filename_nt = paste0(fig_folder_name, panel_name)

# Top individual outputs to NSC subtypes by neurotransmitter
for (i in unique(NSC$NSC_name)) {
  # Get the NSC_ids for the current NSC_name
  NSC_ids <- NSC %>% filter(NSC_name == i) %>% pull(NSC_id)
  # Filter the data for the current NSC_ids
  NSC_subset <- NSC_output %>% filter(pre_pt_root_id %in% NSC_ids)
  if (nrow(NSC_subset) > 0) {
  b <- ggplot(NSC_subset, aes(x = fct_reorder(post_pt_root_id, desc(n_synapses)), 
                              y = n_synapses, fill = nt_type_stringent_corr)) +
    geom_col() +
    labs(x = "", y = "", title = i) +
    scale_fill_manual(values = neurotransmitter_colors) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, color="black"),
          axis.text.y = element_text(size = 8, color="black"),
          axis.line = element_line(color = "black"),
          plot.title = element_text(size = 11, hjust = 0.5, color="black"),
          legend.title = element_blank(),
          text = element_text(size = 12, color="black")) +
    scale_y_continuous(breaks = waiver(), n.breaks = 5, expand = c(0,0))
  print(b)
  if (write_plots){
    pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename_nt, i,
                                                               ".pdf"))
    ggsave(pdf_path, plot = b, width = 12, height = 10, units = "cm")
    cat("Saved PDF to:", normalizePath(pdf_path), "\n")
  }
  }
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Run python scripts - open .py in R (requires reticulate) or run elsewhere
# otherwise, download files from zenodo and put into: output/neurons_v783/
# Note: These were originally run in Figure 3 and reused here
# ------------------------------------------------------------------------------
# Save brainmesh.obj for 3D plotting -------------------------------------------
# Figure_3_S1_brainmesh.py

# Save NSC skeleton objs for 3D plotting ---------------------------------------
# Figure_3_S1_neuronMesh.py

# ------------------------------------------------------------------------------
# Figure 7 S1 - Individual NSC output synapse plots
# ------------------------------------------------------------------------------
if (write_mesh) {
  panel_name = "_S1_"
  filename_output = paste0(fig_folder_name, panel_name)
  
  # Iterate over unique NSC names and plot synapses and meshes
  for (i in unique(NSC$NSC_name)) {
    # Get the NSC_ids for the current NSC_name
    NSC_ids <- NSC %>% filter(NSC_name == i) %>% pull(NSC_id)
    # Filter the synapses for the current NSC_ids
    synapses_filtered <- synapses_output %>% filter(pre_pt_root_id %in% NSC_ids)
    # Get skeletons names for NSCs
    NSC_ids_skel <- paste0(NSC_ids, ".obj")
    # Plot skeletons and output synapses
    open3d()
    # Set high resolution for the plot
    par3d(windowRect = c(0, 0, 3840, 2160))
    for (j in NSC_ids_skel) {
      NSC_skeletons_tmp <- readOBJ(paste0(PATH_output, "/neurons_v783/NSC_mesh/", j))
      plot3d(
        NSC_skeletons_tmp,
        col = "black",
        add = T,
        WithNodes = F
      )
    }
    plot3d(
      synapses_filtered$pre_x,
      synapses_filtered$pre_y,
      synapses_filtered$pre_z,
      col = syn_col[i],
      size = 0.5,
      type = "s",
      add = TRUE
    )
    # Add the surface model
    brainmesh <- readOBJ(paste0(PATH_output, "/neurons_v783/brainmesh/brainmesh.obj"))
    plot3d(brainmesh,
           add = TRUE,
           alpha = 0.1,
           col = "grey")
    # Adjust the view
    view3d(userMatrix = rotationMatrix(90 * pi / 90, 1, 0, 0), zoom = 0.5)
    
    if (write_plots) {
      # Export as png
      png_filename <- file.path(PATH_output,
                                fig_folder_name,
                                paste0(filename_output, i, ".png"))
      rgl.snapshot(png_filename)
      cat("Saved PNG to:", normalizePath(png_filename), "\n")
    }
    # Close the 3D plot
    close3d()
  }
}
