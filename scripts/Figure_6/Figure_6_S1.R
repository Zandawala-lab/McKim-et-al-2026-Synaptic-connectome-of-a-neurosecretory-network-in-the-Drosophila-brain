# -- Figure 6 S1 ---------------------------------------------------------------
# ------------------------------------------------------------------------------
# Code for post-processing and visualization after running influence analysis

# Steps completed prior to this (see InfluenceAnalysis_README.rtf):
# 1. Downloaded data from Harvard dataverse for BANC v888 edgelist and meta to
     # create sqlite file for use in analysis (see python script)
# 2. IDs of interest as source neurons (see python script)
# 3. Python script containing IDs of interest as source neurons and code to run
    # influence pipeline (v0.4) following steps outlined here:
    # https://github.com/DrugowitschLab/ConnectomeInfluenceCalculator
    # Sensory_influence_BANC.py
# 4. Output files (csv) from the influence analysis are used here
    # Example naming scheme: ENS1_5thresh_influence.csv

# Loads influence .csv's generated & plot heatmaps
# ------------------------------------------------------------------------------

# load packages ----------------------------------------------------------------
library(tidyverse)
library(ggtext)
library(gt)

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
fig_folder_name = "Figure_6"
if(!dir.exists(file.path(PATH_output, fig_folder_name))){
  dir.create(file.path(PATH_output, fig_folder_name))
  cat("Created figure folder at:", normalizePath(file.path(PATH_output,
                                                           fig_folder_name)), "\n")
}

# constant based on Bates et al. 2026
c <- 24 

# ------------------------------------------------------------------------------
# Options for saving plots & csv files
# ------------------------------------------------------------------------------
# When running script for the first time, set both to TRUE:
write_plots = TRUE           # TRUE - save/replicate figure plots
                             # FALSE - plots not saved outside of R
write_csv   = TRUE           # TRUE - save processed data associated w/figures
                             # FALSE - data not saved outside of R

# set colors -------------------------------------------------------------------
nsc_colors <- c(
  m_NSC_unknown = "#BFD739",
  l_NSC_unknown = "#009817",
  m_NSC_DILP    = "#BD0023",
  m_NSC_DH44    = "#ffe200",
  m_NSC_DMS     = "#FE7E00",
  l_NSC_CRZ     = "#8100FF",
  l_NSC_ITP     = "#838383",
  l_NSC_DH31    = "#00B4FF",
  SEZ_NSC_CAPA  = "#003BBD",
  SEZ_NSC_Hugin = "#BD00B0"
)

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

scaled_heatmap_palette <- colorRampPalette(c('#2166AC', '#4393C3', '#92C5DE', 
                                             '#F4A582', '#D6604D', '#B2182B'))(100)

# ------------------------------------------------------------------------------
# Load BANC ids for reference
# ------------------------------------------------------------------------------
# provided in input/ (zenodo); can be obtained from BANC meta (Figure 1 S1)
BANC_NSC <- read_delim(file.path(PATH_input, "BANC_NSC_brain_meta_v888.csv"),
                       delim = ",",
                       escape_double = FALSE,
                       col_types = cols(pt_root_id = col_character()),
                       trim_ws = TRUE)
# ------------------------------------------------------------------------------
# Find all CSV files in input folder and subfolders
# ------------------------------------------------------------------------------

csv_files <- list.files(
  path = file.path(PATH_input, "BANC_influence"),
  pattern = "\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

# ------------------------------------------------------------------------------
# Read and organize data by folder
# ------------------------------------------------------------------------------

all_data_by_folder <- list()

for (csv_file in csv_files) {
  
  cat("Found:", csv_file, "\n")
  
  # Extract folder name from file path
  relative_path <- sub(paste0("^", PATH_input), "", csv_file)
  relative_path <- sub("^/", "", relative_path)
  path_parts <- strsplit(relative_path, "/")[[1]]
  
  if (length(path_parts) > 1) {
    folder_name <- path_parts[1]
  } else {
    folder_name <- "root"
  }
  
  # Extract dataset name from filename
  filename <- basename(csv_file)
  dataset_name <- tools::file_path_sans_ext(filename)
  dataset_name <- gsub("_influence_meta", "", dataset_name)
  
  # Read the influence metadata
  noi.influence_meta <- read_csv(csv_file, 
                                 col_types = cols(id = col_character(),
                                                  manc_match = col_character(),
                                                  hemibrain_match = col_character()))
  
  # Add source column
  noi.influence_meta$source_file <- dataset_name
  
  # Initialize list for this folder if it doesn't exist
  if (is.null(all_data_by_folder[[folder_name]])) {
    all_data_by_folder[[folder_name]] <- list()
  }
  
  # Add to list
  all_data_by_folder[[folder_name]][[length(all_data_by_folder[[folder_name]]) + 1]] <- noi.influence_meta
  
}

# ------------------------------------------------------------------------------
# Combine data within each folder group
# ------------------------------------------------------------------------------

combined_data_by_folder <- lapply(all_data_by_folder, function(df_list) {
  bind_rows(df_list)
})
rm(all_data_by_folder)

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Analysis by Endocrine (NSC class)
# ------------------------------------------------------------------------------

# Normalize influence calc & filter by endocrine -------------------------------

infl_nsc <- combined_data_by_folder[["BANC_influence"]] %>%
  group_by(source_file) %>%
  mutate(num_seeds = sum(is_seed)) %>%
  filter(super_class == "visceral_circulatory") %>%  
  filter(!is_seed) %>% #remove self influence
  group_by(source_file, cell_type) %>%    
  summarise(mean = log(sum(`Influence_score_(unsigned)`) / 
                         (n_distinct(id)*unique(num_seeds))) + c,
            .groups = "drop") %>% # change num_seeds to 1 to norm only by targets. If warned with reframe(), grouping went wrong and unique returns > 1 val
  arrange(desc(mean)) %>%
  mutate(
    cell_type = coalesce(cell_type, "unknown"),
    cell_type = gsub("_", " ", cell_type)
  )

# rename categories as unknown
infl_nsc <- infl_nsc %>%
  mutate(cell_type = case_when(
    cell_type == "m-NSC" ~ "m NSC unknown",
    cell_type == "l NSC" ~ "l NSC unknown",
    TRUE ~ cell_type
  ))

review <- infl_nsc %>% filter(is.nan(mean) | is.infinite(mean))
unique(review$source_file)
# none to remove
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 6 S1A - NSC HEATMAP PLOT
# ------------------------------------------------------------------------------
panel_name = "_S1A"
output_filename <- paste0(fig_folder_name, panel_name)

# Fill in missing combinations
data_complete <- infl_nsc %>% ungroup() %>%
  complete(cell_type, source_file, fill = list(mean = NA))

# tidy names
data_complete <- data_complete %>%
  mutate(source_file = str_remove(source_file, "_5thresh_influence")) %>%
  mutate(source_file = str_replace_all(source_file, "_", " "))

data_NSC <- data_complete %>%
  filter(cell_type %in% c("SEZ NSC Hugin", "SEZ NSC CAPA", "l NSC DH31", 
                            "l NSC CRZ","m NSC DMS", "m NSC DH44", "m NSC DILP", 
                            "l NSC unknown", "m NSC unknown")) %>%
  filter(!is.na(source_file)) 

# source_files <- unique(data_complete$source_file) # to view the names
# set the order to show on plot
data_NSC$source_file <- factor(
  data_NSC$source_file,
  levels = rev(c("internal taste sensillum gustatory neuron", "GRN LB4a","GRN LB4b",
                 "GRN amino acids Ir94e","GRN aversive","GRN bitter Gr33a",
                 "GRN heavy metal Ir47a","GRN sugar Gr64f","GRN sugar low salt Gr64f Ir56b",
                 "GRN water ppk28","taste peg gustatory neuron", 
                 "ENS1","ENS2","ENS3","ENS4","ENS5",
                 "ORN Aversive","ORN EggLaying","ORN Food","ORN Pheromonal","ORN Unknown",
                 "internal sensillum tactile neuron", "taste peg tactile neuron", 
                 "bristle neuron","chordotonal organ neuron", 
                 "hygrosensory receptor neuron", "thermosensory receptor neuron"
  )))

# rename
levels(data_NSC$source_file) <- levels(data_NSC$source_file) %>%
  str_replace("^GRN ", "gustatory ") %>%
  str_replace("^ENS(\\d)", "enteric (ens \\1)") %>%
  str_replace("^ORN ", "olfactory ") %>%
  str_replace("EggLaying", "egg-laying") %>%
  str_replace("(olfactory )([A-Z])", function(m) paste0(str_extract(m, "olfactory "), 
                                                        tolower(str_extract(m, "[A-Z]"))))

# NSC name order is based on anatomical location (not including nerve cord)
data_NSC$cell_type <- fct_relevel(
  data_NSC$cell_type,
  "SEZ NSC Hugin", "SEZ NSC CAPA", "l NSC DH31", 
  "l NSC CRZ","m NSC DMS", "m NSC DH44", "m NSC DILP", 
  "l NSC unknown", "m NSC unknown", after = 0)

# Normalize scale 0–1 within neuron_type
data_NSC <- data_NSC %>%
  mutate(mean = (mean - min(mean, na.rm = TRUE)) /
           (max(mean, na.rm = TRUE) - min(mean, na.rm = TRUE)))

# color setup for plots (requires ggtext)
colored_labels <- setNames(
  paste0("<span style='color:", nsc_colors[gsub(" ", "_", levels(data_NSC$cell_type))], "'>", 
         custom_labels_map[gsub(" ", "_", levels(data_NSC$cell_type))], "</span>"),
  levels(data_NSC$cell_type)
)

# Create combined heatmap
p <- ggplot(data_NSC, aes(x = cell_type, y = source_file, fill = mean)) +
  geom_raster() +
  scale_x_discrete(expand = c(0, 0), limits = levels(data_NSC$cell_type), 
                   labels = colored_labels) +
  scale_y_discrete(expand = c(0, 0), limits = levels(data_NSC$source_file)) +
  scale_fill_gradientn(
    colours = scaled_heatmap_palette,
    breaks = c(min(data_NSC$mean, na.rm = TRUE),
               max(data_NSC$mean, na.rm = TRUE)),
    labels = round(c(min(data_NSC$mean, na.rm = TRUE),
                     max(data_NSC$mean, na.rm = TRUE)), 1),
    na.value = "white"
  ) +
  labs(title = "",
       x = "Target",
       y = "Source",
       fill = "Mean Influence") +
  theme_bw() +
  coord_fixed() +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5 
    )
  ) +
  theme(axis.text.x = element_markdown(angle = 270, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(color = "black"),
        panel.grid = element_blank(),
        axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),        
        legend.title.position = "top",
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(1, "cm"),
        legend.position = "right",
        legend.direction = "vertical")

# Draw y-axis lines to define category boundaries
group_boundaries <- c(2, 6, 11, 16)  # levels(data_complete$source_file)

p <- p + geom_hline(yintercept = group_boundaries + 0.5, color = "black", linewidth = 0.5)

print(p)

if(write_plots){
  output_filename_pdf <- file.path(PATH_output, fig_folder_name, 
                                   paste0(output_filename, ".pdf"))
  ggsave(output_filename_pdf, plot = p, width = 12*2, height = 14*1.2, units = "cm")
  cat("Saved PDF to:", normalizePath(output_filename_pdf), "\n")
  
  output_filename_png <- file.path(PATH_output, fig_folder_name, 
                                   paste0(output_filename, ".png"))
  ggsave(output_filename_png, plot = p, width = 12*2, height = 14*1.2, 
         units = "cm", bg = "transparent")
  cat("Saved PNG to:", normalizePath(output_filename_png), "\n")
  
  tif_path <- paste0(tools::file_path_sans_ext(output_filename_pdf), ".tif")
  ggsave(tif_path, plot = p, width = 12*2, height = 14*1.2, units = "cm", 
         bg = "transparent", device = "tiff", dpi = 300)
  cat("Saved TIF to:", normalizePath(tif_path), "\n")
}

# ------------------------------------------------------------------------------
# Analysis by Class
# ------------------------------------------------------------------------------

# if class is missing, use super_class
combined_data_by_folder[["BANC_influence"]] <- combined_data_by_folder[["BANC_influence"]] %>%
  mutate(cell_class = coalesce(cell_class, super_class))

# Add labels -------------------------------------------------------------------

# class labels for endocrine cells are only pars_intercerebralis and pars_lateralis
# "pars_intercerebralis_neurosecretory_cell", "pars_lateralis_neurosecretory_cell"
# "subesophageal_zone_neurosecretory_cell"     
# we want to group them anatomically, this code labels based on medial, lateral and SEZ

banc_nsc_lookup <- BANC_NSC %>%
  select(pt_root_id, cell_type) %>%
  mutate(cell_class_new = case_when(
    startsWith(cell_type, "l") ~ "lateral NSC",
    startsWith(cell_type, "m") ~ "medial NSC",
    startsWith(cell_type, "SEZ") ~ "SEZ NSC"
  ))

combined_data_by_folder[["BANC_influence"]] <- combined_data_by_folder[["BANC_influence"]] %>%
  left_join(banc_nsc_lookup %>% select(pt_root_id, cell_class_new), by = c("id" = "pt_root_id")) %>%
  mutate(cell_class = coalesce(cell_class_new, cell_class)) %>%
  select(-cell_class_new)

# check if any in combined_data not in BANC NSC ids
nsc_check_ids <- combined_data_by_folder[["BANC_influence"]] %>%
  filter(cell_class %in% c("pars_intercerebralis_neurosecretory_cell",
                           "pars_lateralis_neurosecretory_cell",
                           "subesophageal_zone_neurosecretory_cell")) %>%
  distinct(id) %>%
  filter(!id %in% BANC_NSC$pt_root_id)

# check if any in combined_data not in BANC NSC ids
nsc_check <- combined_data_by_folder[["BANC_influence"]] %>%
  filter(cell_class %in% c("pars_intercerebralis_neurosecretory_cell",
                           "pars_lateralis_neurosecretory_cell",
                           "subesophageal_zone_neurosecretory_cell")) %>%
  filter(!id %in% BANC_NSC$pt_root_id)

# remove additional endocrine cells not in McKim et al (also in Fig 3 S4)
combined_data_by_folder[["BANC_influence"]] <- combined_data_by_folder[["BANC_influence"]] %>%
  filter(!cell_type %in% c("CB2748", "DNg28")) %>%
  filter(id != "720575941611801674") # one of the nine ITP neurons is not ITP - remove 

# Normalize influence calc & filter by endocrine -------------------------------

infl_class <- combined_data_by_folder[["BANC_influence"]] %>%
  group_by(source_file) %>%
  mutate(num_seeds = sum(is_seed)) %>%
  filter(!is_seed) %>% #remove self influence
  group_by(source_file, cell_class) %>%
  summarise(mean = log(sum(`Influence_score_(unsigned)`) / 
                         (n_distinct(id)*unique(num_seeds))) + c,
            .groups = "drop") %>% # change num_seeds to 1 to norm only by targets. If warned with reframe(), grouping went wrong and unique returns > 1 val
  arrange(desc(mean)) %>%
  mutate(
    cell_class = coalesce(cell_class, "unknown"),
    cell_class = gsub("_", " ", cell_class)
  )

# need combos of cell class and source file
review <- infl_class %>% filter(is.nan(mean) | is.infinite(mean))
# unique(review$cell_class) # sensory_ascending and trachea

# remove them
infl_class <- infl_class %>% filter(!is.nan(mean) & !is.infinite(mean))

# ------------------------------------------------------------------------------
# Figure 6 S1B - CLASS HEATMAP PLOT
# ------------------------------------------------------------------------------
panel_name = "_S1B"
output_filename <- paste0(fig_folder_name, panel_name)

# Fill in missing combinations
data_complete <- infl_class %>% ungroup() %>%
  complete(cell_class, source_file, fill = list(mean = NA))

# tidy names
data_complete <- data_complete %>%
  mutate(source_file = str_remove(source_file, "_5thresh_influence"),
         source_file = str_replace_all(source_file, "_", " "))

# source_files <- unique(data_complete$source_file) # to view the names

# set the order to show on plot
# source_files <- unique(data_complete$source_file) # to view the names
# set the order to show on plot
data_complete$source_file <- factor(
  data_complete$source_file,
  levels = rev(c("internal taste sensillum gustatory neuron", "GRN LB4a","GRN LB4b",
                 "GRN amino acids Ir94e","GRN aversive","GRN bitter Gr33a",
                 "GRN heavy metal Ir47a","GRN sugar Gr64f","GRN sugar low salt Gr64f Ir56b",
                 "GRN water ppk28","taste peg gustatory neuron", 
                 "ENS1","ENS2","ENS3","ENS4","ENS5",
                 "ORN Aversive","ORN EggLaying","ORN Food","ORN Pheromonal","ORN Unknown",
                 "internal sensillum tactile neuron","taste peg tactile neuron", 
                 "bristle neuron","chordotonal organ neuron", 
                 "hygrosensory receptor neuron","thermosensory receptor neuron"
  )))

# Label order for class - put the NSCs at the end
data_complete$cell_class <- fct_relevel(
  data_complete$cell_class,
  "SEZ NSC", "medial NSC", "lateral NSC",
  after = Inf
)

# rename
levels(data_complete$source_file) <- levels(data_complete$source_file) %>%
  str_replace("^GRN ", "gustatory ") %>%
  str_replace("^ENS(\\d)", "enteric (ens \\1)") %>%
  str_replace("^ORN ", "olfactory ") %>%
  str_replace("EggLaying", "egg-laying") %>%
  str_replace("(olfactory )([A-Z])", function(m) paste0(str_extract(m, "olfactory "), 
                                                        tolower(str_extract(m, "[A-Z]"))))


# Normalize scale 0–1
data_complete <- data_complete %>%
  mutate(mean = (mean - min(mean, na.rm = TRUE)) /
           (max(mean, na.rm = TRUE) - min(mean, na.rm = TRUE)))


# Create combined heatmap
p2 <- ggplot(data_complete, aes(x = cell_class, y = source_file, fill = mean)) +
  geom_raster() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0), limits = levels(data_complete$source_file)) +
  scale_fill_gradientn(
    colours = scaled_heatmap_palette,
    breaks = c(min(data_complete$mean, na.rm = TRUE),
               max(data_complete$mean, na.rm = TRUE)),
    labels = round(c(min(data_complete$mean, na.rm = TRUE),
                     max(data_complete$mean, na.rm = TRUE)), 1),
    na.value = "white"
  ) +
  labs(title = "",
       x = "target",
       y = "source",
       fill = "Mean influence") +
  theme_minimal() +
  coord_fixed() +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5 
    )
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.grid = element_blank(),
        axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.title.position = "bottom",
        legend.key = element_rect(fill = "transparent", color = NA),        
        legend.key.width = unit(1, "null"),
        legend.key.height = unit(0.3, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal")

# Draw y-axis lines to define category boundaries
group_boundaries <- c(2, 6, 11, 16)  # levels(data_complete$source_file)

p2 <- p2 + geom_hline(yintercept = group_boundaries + 0.5, color = "black", linewidth = 0.5)

print(p2)

if(write_plots){
  output_filename_pdf <- file.path(PATH_output, fig_folder_name, paste0(output_filename, ".pdf"))
  ggsave(output_filename_pdf, plot = p2, width = 12*3.5, height = 14*2.7, units = "cm")
  cat("Saved PDF to:", normalizePath(output_filename_pdf), "\n")
  
  output_filename_png <- file.path(PATH_output, fig_folder_name, paste0(output_filename, ".png"))
  ggsave(output_filename_png, plot = p2, width = 12*3.5, height = 14*2.7, units = "cm")
  cat("Saved PNG to:", normalizePath(output_filename_png), "\n")
  
  tif_path <- paste0(tools::file_path_sans_ext(output_filename_pdf), ".tif")
  ggsave(tif_path, plot = p, width = 12*3.5, height = 14*2.7, units = "cm", 
         bg = "transparent", device = "tiff", dpi = 300)
  cat("Saved TIF to:", normalizePath(tif_path), "\n")
}

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 6 S1C - NSC RANK BASED ON THE CLASS ANALYSIS
# ------------------------------------------------------------------------------
panel_name = "_S1C"
csv_filename <- file.path(PATH_output, fig_folder_name, 
                          paste0(fig_folder_name, panel_name, ".csv"))
nsc_ranks <- infl_class %>%
  group_by(source_file) %>%
  mutate(
    rank = rank(desc(mean), ties.method = "min"),
    total = n(),
    percentile = round(rank / total * 100, 1)
  ) %>%
  ungroup() %>%
  filter(cell_class %in% c("SEZ NSC", "lateral NSC", "medial NSC")) %>%
  select(source_file, cell_class, mean, rank, total, percentile) %>%
  arrange(source_file, rank)

if (write_csv) {
  write.csv(nsc_ranks, file = csv_filename, row.names = FALSE)
  cat("Saved CSV to:", normalizePath(csv_filename), "\n")
}

# ------------------------------------------------------------------------------
# Figure 6 S1C - VISUALIZE AS A HEATMAP WITH NUMBERS
# ------------------------------------------------------------------------------

# tidy names
nsc_ranks <- nsc_ranks %>%
  mutate(source_file = str_remove(source_file, "_5thresh_influence")) %>%
  mutate(source_file = str_replace_all(source_file, "_", " "))

nsc_ranks <- nsc_ranks %>%
  mutate(source_file = str_remove(source_file, "_5thresh_influence")) %>%
  mutate(source_file = str_replace_all(source_file, "_", " ")) %>%
  mutate(source_file = str_replace(source_file, "^GRN ", "gustatory ")) %>%
  mutate(source_file = str_replace(source_file, "^ENS(\\d)", "enteric (ens \\1)")) %>%
  mutate(source_file = str_replace(source_file, "^ORN ", "olfactory ")) %>%
  mutate(source_file = str_replace(source_file, "EggLaying", "egg-laying")) %>%
  mutate(source_file = gsub("(olfactory )([A-Z])", "\\1\\L\\2", source_file, perl = TRUE))

source_order <- c("internal taste sensillum gustatory neuron", "gustatory LB4a","gustatory LB4b",
                  "gustatory amino acids Ir94e","gustatory aversive","gustatory bitter Gr33a",
                  "gustatory heavy metal Ir47a","gustatory sugar Gr64f",
                  "gustatory sugar low salt Gr64f Ir56b",
                  "gustatory water ppk28","taste peg gustatory neuron", 
                  "enteric (ens 1)","enteric (ens 2)","enteric (ens 3)",
                  "enteric (ens 4)","enteric (ens 5)",
                  "olfactory aversive","olfactory egg-laying","olfactory food",
                  "olfactory pheromonal","olfactory unknown",
                  "internal sensillum tactile neuron", "taste peg tactile neuron", 
                  "bristle neuron","chordotonal organ neuron", 
                  "hygrosensory receptor neuron","thermosensory receptor neuron")

# set the order for plotting
nsc_class_order <- c("SEZ NSC", "medial NSC", "lateral NSC", "total")

p_ranks <- nsc_ranks %>%
  select(source_file, cell_class, rank, total) %>%
  mutate(source_file = factor(source_file, levels = source_order),
         cell_class = factor(cell_class, levels = nsc_class_order)) %>%
  bind_rows(
    nsc_ranks %>%
      distinct(source_file, total) %>%
      mutate(cell_class = "total", rank = total)
  ) %>%
  mutate(cell_class = factor(cell_class, levels = nsc_class_order)) %>%
  ggplot(aes(x = cell_class, y = source_file, fill = ifelse(cell_class == "total", NA, rank))) +
  geom_raster() +
  geom_text(aes(label = rank), color = "black", size = 3) +
  scale_fill_gradient(low = "black", high = "white", na.value = "white") +
  geom_text(aes(label = rank, color = ifelse(rank < max(rank)/2, "white", "black")), size = 3) +
  scale_color_identity() +
  geom_hline(yintercept = (nrow(distinct(nsc_ranks, source_file)) - c(11,16,21)) + 0.5, 
             color = "red", linewidth = 1) +
  scale_y_discrete(limits = rev(levels(factor(nsc_ranks$source_file, levels = source_order)))) +
  labs(x = "Target", y = "Source", fill = "Rank") +
  theme_minimal() +
  coord_fixed(ratio = 1) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )
print(p_ranks)

if(write_plots){
  base_path <- file.path(PATH_output, fig_folder_name, paste0(fig_folder_name, panel_name))
  
  ggsave(paste0(base_path, ".png"), plot = p_ranks, width = 15, height = 20, units = "cm")
  cat("Saved PNG to:", normalizePath(paste0(base_path, ".png")), "\n")
  
  ggsave(paste0(base_path, ".pdf"), plot = p_ranks, width = 15, height = 20, units = "cm")
  cat("Saved PDF to:", normalizePath(paste0(base_path, ".pdf")), "\n")
  
  ggsave(paste0(base_path, ".tif"), plot = p_ranks, width = 15, height = 20, units = "cm", 
         device = "tiff", dpi = 300)
  cat("Saved TIF to:", normalizePath(paste0(base_path, ".tif")), "\n")
}
