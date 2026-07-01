# -- Figure 6 ------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Code for post-processing and visualization after running influence analysis

# Steps completed prior to this (see InfluenceAnalysis_README.rtf):
# 1. Downloaded data from codex used for connections, classification, etc. to
    # create sqlite file for use in analysis (see python script)
# 2. IDs of interest as source neurons (see python script)
# 3. Python script containing IDs of interest as source neurons and code to run
    # influence pipeline (v0.4) following steps outlined here:
    # https://github.com/DrugowitschLab/ConnectomeInfluenceCalculator
    # scripts/Figure 6/Sensory_influence_FAFB.py
# 4. Output files (csv) from the influence analysis are used here:
    # Example naming scheme: ens1_5thresh_influence.csv

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
# Load data: find all CSV files in input folder and subfolders
# ------------------------------------------------------------------------------

csv_files <- list.files(
  path = file.path(PATH_input, "FAFB_influence"),
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
  noi.influence_meta <- read_csv(csv_file, col_types = cols(id = col_character()))
  
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
# Adding updated FAFB types
# ------------------------------------------------------------------------------

classification <- read_delim(file.path(PATH_input, "classification_v783_02242026.csv"),
                             delim = ",",
                             escape_double = FALSE,
                             col_types = cols(root_id = col_character(), 
                                              flow = col_character()),
                             trim_ws = TRUE)

consolidated_cell_types <- read_delim(file.path(PATH_input,
                                                "consolidated_cell_types_v783_02242026.csv"),
                                      delim = ",",
                                      escape_double = FALSE,
                                      col_types = cols(root_id = col_character(),
                                                       primary_type = col_character(),
                                                       `additional_type(s)` = col_character()),
                                      trim_ws = TRUE)

# rename primary_type to cell_type & additional_types(s) to remove (s) 
consolidated_cell_types <- consolidated_cell_types %>%
  rename(cell_type = primary_type,
         cell_type_additional = `additional_type(s)`)

# join these two by root_id
classification_ct <- left_join(classification, consolidated_cell_types, by = "root_id")
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#  Analysis by Endocrine (NSC class)
# ------------------------------------------------------------------------------

# Normalize influence calc & filter by endocrine -------------------------------

infl_nsc <- combined_data_by_folder[["FAFB_influence"]] %>%
  group_by(source_file) %>%
  mutate(num_seeds = sum(is_seed)) %>%
  filter(super_class == "endocrine") %>%  
  filter(!is_seed) %>% # remove self influence
  group_by(source_file, cell_type) %>%    
  summarise(mean = log(sum(`Influence_score_(unsigned)`) / 
                         (n_distinct(id)*unique(num_seeds))) + c,
            .groups = "drop") %>% # change num_seeds to 1 to norm only by targets. If warned with reframe(), grouping went wrong and unique returns > 1 val
  arrange(desc(mean)) %>%
  mutate(
    cell_type = coalesce(cell_type, "unknown"),
    cell_type = gsub("_", " ", cell_type)
  )
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 6B - NSC HEATMAP PLOT
# ------------------------------------------------------------------------------
panel_name = "B"
output_filename <- paste0(fig_folder_name, panel_name)

# Fill in missing combinations
data_complete <- infl_nsc %>% ungroup() %>%
  complete(cell_type, source_file, fill = list(mean = NA))

# tidy names
data_complete <- data_complete %>%
  mutate(source_file = str_remove(source_file, "_5thresh_influence")) %>%
  mutate(source_file = str_replace_all(source_file, "_", " "))

# source_files <- unique(data_complete$source_file) # to view the names
# set the order to show on plot
data_complete$source_file <- factor(
  data_complete$source_file,
  levels = rev(c("gust aPN sensory grp1", "gust aPN sensory grp2", "gust PN sensory grp2", 
                 "gust PN sensory grp3",
                 "gust SA VTV pro meso meta", "gust bitter", "gust low salt", 
                 "gust sugar water", "gust taste peg",
                 "ens1", "ens2", "ens3", "ens4", "ens5",
                 "ORNs aversive", "ORNs eggLaying", "ORNs food", "ORNs pheromonal", 
                 "ORNs unknown",
                 "mechano aPhN", "mechano AN", "mechano MxLbN", "mechano ON",
                 "hygrosensory", "thermosensory"
  )))

levels(data_complete$source_file) <- levels(data_complete$source_file) %>%
  str_replace("^gust ", "gustatory ") %>%
  str_replace("sugar water", "sugar/water") %>%
  str_replace("^ens(\\d)", "enteric (ens \\1)") %>%
  str_replace("^mechano ", "mechanosensory ") %>%
  str_replace("^ORNs ", "olfactory ") %>%
  str_replace("grp(\\d)", "grp \\1") %>%
  str_replace("eggLaying", "egg-laying")

# Remove ITP before building the plot (infl is low)
data_complete <- data_complete %>% filter(cell_type != "l NSC ITP")

# NSC name order
data_complete$cell_type <- fct_relevel(
  data_complete$cell_type,
  "SEZ NSC Hugin", "SEZ NSC CAPA", "l NSC DH31", "l NSC CRZ",
  "m NSC DMS", "m NSC DH44", "m NSC DILP", "l NSC unknown", "m NSC unknown",
  after = 0)

# Normalize scale 0–1
data_complete <- data_complete %>%
  mutate(mean = (mean - min(mean, na.rm = TRUE)) /
           (max(mean, na.rm = TRUE) - min(mean, na.rm = TRUE))) 

# color setup for plots (requires ggtext)
colored_labels <- setNames(
  paste0("<span style='color:", nsc_colors[gsub(" ", "_",
                                                levels(data_complete$cell_type))], "'>", 
         custom_labels_map[gsub(" ", "_", levels(data_complete$cell_type))], "</span>"),
  levels(data_complete$cell_type)
)

p <- ggplot(data_complete, aes(x = cell_type, y = source_file, fill = mean)) +
  geom_raster() +
  scale_x_discrete(expand = c(0, 0), limits = levels(data_complete$cell_type), 
                   labels = colored_labels) +
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

# ------------------------------------------------------------------------------
# Analysis by Class
# ------------------------------------------------------------------------------

# update with newest classification and cell type info from codex
combined_data_by_folder[["FAFB_influence"]] <- combined_data_by_folder[["FAFB_influence"]] %>%
  select(-class, -super_class) %>%
  left_join(classification_ct %>% select(root_id, class, super_class), 
            by = c("id" = "root_id"))

# if class is missing, use super_class
combined_data_by_folder[["FAFB_influence"]] <- combined_data_by_folder[["FAFB_influence"]] %>%
  mutate(class = coalesce(class, super_class))


# Add labels -------------------------------------------------------------------

# class labels for endocrine cells are only pars_intercerebralis and pars_lateralis
# we want to group them anatomically, this code labels based on medial, lateral and SEZ

combined_data_by_folder[["FAFB_influence"]] <- combined_data_by_folder[["FAFB_influence"]] %>%
  mutate(class = case_when(
    id %in% c("720575940627293201", "720575940643847022","720575940638394446", 
              "720575940618736797", "720575940620829878", "720575940637851069") ~ "SEZ NSC",
    id %in% c("720575940622872870", "720575940650527222", "720575940623011845",
              "720575940635861146", "720575940603765280", "720575940612923390", 
              "720575940628363820", "720575940611254681", "720575940616153371", 
              "720575940631884883", "720575940643539566", "720575940618385718",
              "720575940626619398", "720575940622897639", "720575940625379859",
              "720575940615911380", "720575940623081400", "720575940629974531",
              "720575940615043112", "720575940628199802", "720575940647963513",
              "720575940626375164", "720575940618694827", "720575940632436307",
              "720575940613520351", "720575940604934502", "720575940614623455", 
              "720575940624064295", "720575940608157397", "720575940620628957",
              "720575940638274688", "720575940614107298", "720575940631570380",
              "720575940629938895", "720575940644562286", "720575940628756951", 
              "720575940624606502", "720575940623401703", "720575940610733797", 
              "720575940610467682") ~ "medial NSC",
    id %in% c("720575940636633848", "720575940620332141", "720575940619845552",
              "720575940620503512", "720575940614254534", "720575940622880525",
              "720575940646225006", "720575940642666440", "720575940629808911",
              "720575940629872971", "720575940629924091", "720575940613850262",
              "720575940610895474", "720575940627436035", "720575940609268316",
              "720575940630046506", "720575940631592017", "720575940625721118",
              "720575940625796839", "720575940644303639", "720575940619186240",
              "720575940620926740", "720575940625054901", "720575940622644872",
              "720575940624041203", "720575940628362993", "720575940611186034",
              "720575940626841354", "720575940626878214", "720575940617626235",
              "720575940633416275", "720575940634288735", "720575940625878909",
              "720575940632080429") ~ "lateral NSC",
    id %in% c("720575940619652033", "720575940606239666", "720575940627274578", 
              "720575940645108936", "720575940644676900", "720575940627989609", 
              "720575940634085786", "720575940626767152", "720575940648377209", 
              "720575940612256938") ~ "ens1",
    id %in% c("720575940607152476", "720575940639972816", "720575940627107271",
              "720575940630806267", "720575940613004402", "720575940631462969", 
              "720575940631205580", "720575940605161388", "720575940605682790", 
              "720575940659479937") ~ "ens2",
    id %in% c("720575940610897906", "720575940631070239", "720575940632743649", 
              "720575940630574659", "720575940650760825", "720575940620528990",
              "720575940604782624", "720575940635501560", "720575940638086256",
              "720575940646039092", "720575940611388595", "720575940633895469",
              "720575940634755041", "720575940633129683", "720575940633194029",
              "720575940645524003", "720575940622144648", "720575940623056422",
              "720575940630654539", "720575940625471896", "720575940641265549",
              "720575940637896245", "720575940620983637", "720575940607166731",
              "720575940621207517", "720575940631572141", "720575940625941022",
              "720575940619316603", "720575940633129171") ~ "ens3",
    id %in% c("720575940645493283", "720575940626396803", "720575940625871181", 
              "720575940624910373", "720575940613347474", "720575940620910716", 
              "720575940624905572", "720575940611959011", "720575940633759840",
              "720575940653233569", "720575940639786830", "720575940626311049"
              ) ~ "ens4",
    id %in% c("720575940631963821", "720575940622198772", "720575940633579360",
              "720575940638357813", "720575940631049064", "720575940614641202",
              "720575940630434616", "720575940624724458", "720575940632040749",
              "720575940624502269", "720575940625785731", "720575940621106465",
              "720575940623951591", "720575940606163842", "720575940622575754",
              "720575940604355872", "720575940617034713", "720575940629887695"
              ) ~ "ens5",
    TRUE ~ class  # keep existing class for everything else
  ))

# Normalize influence calc & filter by endocrine -------------------------------

infl_class <- combined_data_by_folder[["FAFB_influence"]] %>%
  group_by(source_file) %>%
  mutate(num_seeds = sum(is_seed)) %>%
  filter(!is_seed) %>% #remove self influence
  group_by(source_file, class) %>%
  summarise(mean = log(sum(`Influence_score_(unsigned)`) / 
                         (n_distinct(id)*unique(num_seeds))) + c,
            .groups = "drop") %>% # change num_seeds to 1 to norm only by targets. If warned with reframe(), grouping went wrong and unique returns > 1 val
  arrange(desc(mean)) %>%
  mutate(
    class = coalesce(class, "unknown"),
    class = gsub("_", " ", class)
  )
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure C - CLASS HEATMAP PLOT
# ------------------------------------------------------------------------------
panel_name = "C"
output_filename <- paste0(fig_folder_name, panel_name)

# Fill in missing combinations
data_complete <- infl_class %>% ungroup() %>%
  complete(class, source_file, fill = list(mean = NA))

# tidy names
data_complete <- data_complete %>%
  mutate(source_file = str_remove(source_file, "_5thresh_influence")) %>%
  mutate(source_file = str_replace_all(source_file, "_", " "))

# source_files <- unique(data_complete$source_file) # to view the names

data_complete$source_file <- factor(
  data_complete$source_file,
  levels = rev(c("gust aPN sensory grp1", "gust aPN sensory grp2", "gust PN sensory grp2", 
                 "gust PN sensory grp3",
                 "gust SA VTV pro meso meta", "gust bitter", "gust low salt",
                 "gust sugar water", "gust taste peg", 
                 "ens1", "ens2", "ens3", "ens4", "ens5",
                 "ORNs aversive", "ORNs eggLaying", "ORNs food", "ORNs pheromonal", 
                 "ORNs unknown",
                 "mechano aPhN", "mechano AN", "mechano MxLbN", "mechano ON",
                 "hygrosensory", "thermosensory"
  )))

levels(data_complete$source_file) <- levels(data_complete$source_file) %>%
  str_replace("^gust ", "gustatory ") %>%
  str_replace("sugar water", "sugar/water") %>%
  str_replace("^ens(\\d)", "enteric (ens \\1)") %>%
  str_replace("^mechano ", "mechanosensory ") %>%
  str_replace("^ORNs ", "olfactory ") %>%
  str_replace("grp(\\d)", "grp \\1") %>%
  str_replace("eggLaying", "egg-laying")

# Label order for class - put the NSCs at the end
data_complete$class <- fct_relevel(
  data_complete$class,
  "SEZ NSC", "medial NSC", "lateral NSC",
  after = Inf
)

# Normalize scale 0–1
data_complete <- data_complete %>%
  mutate(mean = (mean - min(mean, na.rm = TRUE)) /
           (max(mean, na.rm = TRUE) - min(mean, na.rm = TRUE)))


# Create combined heatmap
p2 <- ggplot(data_complete, aes(x = class, y = source_file, fill = mean)) +
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
       x = "Target",
       y = "Source",
       fill = "Mean Influence") +
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
        axis.ticks = element_line(color = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.background = element_blank(),
        legend.title.position = "bottom",
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.3, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal")

# Draw y-axis lines to define category boundaries
group_boundaries <- c(2, 6, 11, 16)  # levels(data_complete$source_file)

p2 <- p2 + geom_hline(yintercept = group_boundaries + 0.5, color = "black", linewidth = 0.5)

print(p2)


if(write_plots){
  output_filename_pdf <- file.path(PATH_output, fig_folder_name, paste0(output_filename, ".pdf"))
  ggsave(output_filename_pdf, plot = p2, width = 12*2, height = 14*1.2, units = "cm")
  cat("Saved PDF to:", normalizePath(output_filename_pdf), "\n")
  
  output_filename_png <- file.path(PATH_output, fig_folder_name, paste0(output_filename, ".png"))
  ggsave(output_filename_png, plot = p2, width = 12*2, height = 14*1.2, units = "cm")
  cat("Saved PNG to:", normalizePath(output_filename_png), "\n")
  
  tif_path <- paste0(tools::file_path_sans_ext(output_filename_pdf), ".tif")
  ggsave(tif_path, plot = p, width = 12*2, height = 14*1.2, units = "cm", 
         bg = "transparent", device = "tiff", dpi = 300)
  cat("Saved TIF to:", normalizePath(tif_path), "\n")
}

# ------------------------------------------------------------------------------
# Figure 6D - NSC RANK BASED ON THE CLASS ANALYSIS
# ------------------------------------------------------------------------------
panel_name = "D"
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
  filter(class %in% c("SEZ NSC", "lateral NSC", "medial NSC")) %>%
  select(source_file, class, mean, rank, total, percentile) %>%
  arrange(source_file, rank)

if (write_csv) {
  write.csv(nsc_ranks, file = csv_filename, row.names = FALSE)
  cat("Saved CSV to:", normalizePath(csv_filename), "\n")
}

# ------------------------------------------------------------------------------
# Figure D - VISUALIZE AS A HEATMAP WITH NUMBERS
# ------------------------------------------------------------------------------

# tidy names
nsc_ranks <- nsc_ranks %>%
  mutate(source_file = str_remove(source_file, "_5thresh_influence")) %>%
  mutate(source_file = str_replace_all(source_file, "_", " ")) %>%
  mutate(source_file = source_file %>%
           str_replace("^gust ", "gustatory ") %>%
           str_replace("sugar water", "sugar/water") %>%
           str_replace("^ens(\\d)", "enteric (ens \\1)") %>%
           str_replace("^mechano ", "mechanosensory ") %>%
           str_replace("^ORNs ", "olfactory ") %>%
           str_replace("grp(\\d)", "grp \\1") %>%
           str_replace("eggLaying", "egg-laying"))

class_order <- c("SEZ NSC", "medial NSC", "lateral NSC", "total")

source_order <- c("gustatory aPN sensory grp 1", "gustatory aPN sensory grp 2", 
                  "gustatory PN sensory grp 2", "gustatory PN sensory grp 3",
                  "gustatory SA VTV pro meso meta", "gustatory bitter", 
                  "gustatory low salt", "gustatory sugar/water", "gustatory taste peg",
                  "enteric (ens 1)", "enteric (ens 2)", "enteric (ens 3)", 
                  "enteric (ens 4)", "enteric (ens 5)",
                  "olfactory aversive", "olfactory egg-laying", "olfactory food", 
                  "olfactory pheromonal", "olfactory unknown",
                  "mechanosensory aPhN", "mechanosensory AN", "mechanosensory MxLbN", 
                  "mechanosensory ON",
                  "hygrosensory", "thermosensory")

p_ranks <- nsc_ranks %>%
  select(source_file, class, rank, total) %>%
  mutate(source_file = factor(source_file, levels = source_order)) %>%
  bind_rows(
    nsc_ranks %>%
      distinct(source_file, total) %>%
      mutate(class = "total", rank = total)
  ) %>%
  mutate(class = factor(class, levels = class_order)) %>%
  ggplot(aes(x = class, y = source_file, fill = ifelse(class == "total", NA, rank))) +
  geom_raster() +
  geom_text(aes(label = rank), color = "black", size = 3) +
  scale_fill_steps(low = "black", high = "white", na.value = "white",
                   limits = c(0, 40), breaks = c(0, 10, 20, 30, 40)) +
  geom_text(aes(label = rank, color = ifelse(rank < max(rank)/2, "white", "black")), size = 3) +
  scale_color_identity() +
  geom_hline(yintercept = (nrow(distinct(nsc_ranks, source_file)) - c(9, 14, 19, 23)) + 0.5,
             color = "red", linewidth = 1) +
  scale_y_discrete(limits = rev(levels(factor(nsc_ranks$source_file, 
                                              levels = source_order)))) +
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
