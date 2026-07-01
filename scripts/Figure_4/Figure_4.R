# --- Figure 4 -----------------------------------------------------------------
# ------------------------------------------------------------------------------
# Code for direct, indirect and sensory/nonsensory interneuron inputs to NSC
# interneuron_input_filtered needs to generated separately by downloading data 
# from Python or download from zenodo: 10.5281/zenodo.20795211
# ------------------------------------------------------------------------------


# load packages ----------------------------------------------------------------
library(natverse)
library(tidyverse)
library(fafbseg)

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
fig_folder_name = "Figure_4"
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

sensory_colors <- c("unknown_sensory" = "black",
                    "enteric"         = "#a85060",
                    "mechanosensory"  = "#edc2c8",
                    "gustatory"       = "#ca4a5a")

# load data --------------------------------------------------------------------
# Filtered and processed data file generated from Figure_3.R and saved to input/
NSC_input <- read_delim(file.path(PATH_input, paste0("NSC_input_classification_v",
                                                     v, ".csv")),
                        delim = ",",
                        escape_double = FALSE,
                        col_types = cols(pre_pt_root_id = col_character(),
                                         post_pt_root_id = col_character()),
                        trim_ws = TRUE)

# Provided with paper and code (zenodo)
NSC <- read_delim(file.path(PATH_input, paste0("NSC_v", v, ".csv")),
                  delim = ",",
                  escape_double = FALSE,
                  col_types = cols(NSC_id = col_character()),
                  trim_ws = TRUE)

# Obtained from coconatfly in Figure_4_S1.R (provided in zenodo)
enteric_groups <- read_csv(file.path(PATH_input, "enteric_groups.csv"),
                          col_types = cols(id = col_character(), 
                                           supervoxel_id = col_character()))

# Prepare data for caveclient to get input IDs to interneurons(NSC_inputs) -----
# remove sensory super_class from direct inputs to NSCs
interneurons_to_NSC <- NSC_input %>% filter(!grepl("sensory", super_class,
                                                   ignore.case = TRUE))

unique_interneurons <- interneurons_to_NSC %>% select(pre_pt_root_id) %>% distinct()
# only done first time through - this file goes to caveclient (python script)
if (write_csv){
  csv_path <- file.path(PATH_input, paste0("unique_interneurons_v", v, ".csv"))
  write.csv(unique_interneurons, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# Load connectivity data -------------------------------------------------------
# This uses python script to identify connections with the IDs from above
# if file already output, skip this and move onto next to pull in the data file
if (!paste0("interneuron_input_filtered_v",v,".csv") %in% input_files) {
  reticulate::source_python("./scripts/Figure_4/Figure_4_preparation.py")
}

# Pull in the data file we saved above with the new info on connections
# input (pre) to interneurons (post)
# If reticulate not setup or errors, download file from zenodo and read in here
interneurons_curated <- read_csv(file.path(PATH_input, 
                                           paste0("tmp/interneuron_input_filtered_v", 
                                                  v, ".csv")),
                                 col_types = cols(
                                   pre_pt_supervoxel_id = col_character(),
                                   pre_pt_root_id = col_character(),
                                   post_pt_supervoxel_id = col_character(),
                                   post_pt_root_id = col_character()
                                 ))


# Check for classification file ------------------------------------------------
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

# adjust classification so enteric is a class
classification <- classification %>%
  mutate(class = ifelse(class == "unknown_sensory" & root_id %in% enteric_groups$id,
                        "enteric", class))

# Data Org & processing --------------------------------------------------------

# Summarize synapses data ------------------------------------------------------
interneurons_input_sum <- interneurons_curated %>%
  group_by(pre_pt_root_id, post_pt_root_id) %>%
  summarise(n_synapses = n(),
            .groups = "drop")

# Save to csv
if (write_csv){
  csv_path <- file.path(PATH_input, "tmp", paste0("interneurons_input_sum_v",
                                                  v, ".csv"))
  write.csv(interneurons_input_sum, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# Join data for further analysis & filter >= 5 synapses
classification_join <- classification
colnames(classification_join) <- c("pre_pt_root_id", colnames(classification)[-1])
interneurons_input <- left_join(interneurons_input_sum[interneurons_input_sum$n_synapses >= 5, ], 
                                classification_join, by = "pre_pt_root_id")
# Save to csv
if (write_csv){
  csv_path <- file.path(PATH_output, fig_folder_name, 
                        paste0("NSC_interneurons_allinput_threshold5_v", v, ".csv"))
  write.csv(interneurons_input, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# Subset classification and adjust col names to include 'int' for interneuron info (these are the post_ids)
classification_int <- select(classification, root_id, super_class, class,
                             cell_type, hemibrain_type)
colnames(classification_int) <- c("post_pt_root_id", "int_super_class","int_class", 
                                  "int_cell_type", "int_hemibrain_type")
interneurons_input <- left_join(interneurons_input, classification_int,
                                by = "post_pt_root_id")

# Remove interneurons (int_) that do not have super_class assigned because these are orphan fragments
interneurons_input_filtered <- interneurons_input %>% drop_na(int_super_class)
# Remove inputs to interneurons which comes from orphan fragments (marked by NA under super_class)
interneurons_input_filtered <- interneurons_input_filtered %>% drop_na(super_class)

# add NSC info from NSC_input
NSC_input_subset <- NSC_input %>% select("pre_pt_root_id","name_post","post_pt_root_id_name")

colnames(NSC_input_subset)[1] <- "post_pt_root_id" # change to match with other df

# join input -> INTs with INTs -> NSCs
ind_sensory_NSC <- left_join(interneurons_input_filtered, NSC_input_subset, 
                             by = "post_pt_root_id", relationship = "many-to-many")

# add enteric label to NSC input
NSC_input <- NSC_input %>%
  mutate(class = ifelse(class == "unknown_sensory" & pre_pt_root_id %in% 
                          enteric_groups$id, "enteric", class))

# For Figure 5 (next), prep for upstream neuron info ---------------------------
al_interneurons_input_filtered <- interneurons_input_filtered[grep("^AL", interneurons_input_filtered$class), ] 
al_interneurons_input_ids <- al_interneurons_input_filtered %>% 
                             select(pre_pt_root_id) %>% 
                             distinct()
# dataframe
if (write_csv){
  csv_path <- file.path(PATH_input, paste0("al_interneurons_input_filtered_v",
                                           v, ".csv"))
  write.csv(al_interneurons_input_filtered, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# ids for python
if (write_csv){
  csv_path <- file.path(PATH_input, paste0("al_interneurons_input_filtered_ids_v",
                                           v, ".csv"))
  write.csv(al_interneurons_input_ids, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Sensory direct inputs to NSCs doughnut plot - Fig 4 Panel A - Direct sensory 
# ------------------------------------------------------------------------------
panel_name = "A_sensorydirect_to_"
filename_sensory_direct = paste0(fig_folder_name, panel_name)

sensory_direct_to_NSC <- NSC_input %>% filter(grepl("sensory", super_class,
                                                    ignore.case = TRUE))
  
for (i in unique(sensory_direct_to_NSC$name_post)) {
  # Get the data for this NSC type (name_post)
  input_df <- unique(sensory_direct_to_NSC[sensory_direct_to_NSC$name_post == i,])
  # count of ids
  neuron_count <- length(unique(input_df$pre_pt_root_id))
  class_grp <- input_df %>%
    group_by(class, pre_pt_root_id) %>%
    summarise(n_count = n(), .groups = "drop")
  
  class_count <- class_grp %>%
    group_by(class) %>% 
    summarise(n_count = length(unique(pre_pt_root_id)), .groups = "drop")

plot <- ggplot(class_count, aes(x=2,y=n_count, fill = class)) + 
    geom_bar(stat = "identity", width =1) +
    scale_fill_manual(values = sensory_colors) +
    coord_polar(theta = "y") +
    theme_void() + 
    xlim(0.5,2.5) +
    theme(legend.position = "none") +
  annotate("text", x = 0.5, y = 0, label = neuron_count, size = 12, hjust = 0.5,
           vjust = 0.5)
  print(plot)
  
  filepath <- file.path(PATH_output, fig_folder_name, paste0(filename_sensory_direct, 
                                                             i, ".pdf"))
  cat("Saving:", filepath, "\n")
  ggsave(filepath, plot = plot, width = 5, height = 5, units = "cm")
  
  }
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Get sensory and nonsensory inputs separately to INTs
sensory_interneurons_to_NSC <- ind_sensory_NSC %>% filter(grepl("sensory", 
                                                                super_class, 
                                                                ignore.case = TRUE)) 


nonsensory_interneurons_to_NSC <-  anti_join(ind_sensory_NSC, sensory_interneurons_to_NSC, 
                                             by = "post_pt_root_id")
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Sensory INTs to NSC - Fig 4 Panel A - Sensory interneuron layer 
# ------------------------------------------------------------------------------
panel_name = "A_sensoryINTs_to_"
filename_sensory_int = paste0(fig_folder_name, panel_name)

for (i in unique(sensory_interneurons_to_NSC$name_post)) {
  # Get the data for this NSC type (name_post)
  input_df <- unique(sensory_interneurons_to_NSC[sensory_interneurons_to_NSC$name_post == i,])
  super_class_grp <- input_df %>%
    group_by(int_super_class, post_pt_root_id) %>%
    summarise(n_count = n(), .groups = "drop")
  # count of ids for the plot
  neuron_count <- length(unique(input_df$post_pt_root_id))
  
  super_class_count <- super_class_grp %>%
    group_by(int_super_class) %>% 
    summarise(n_count = length(unique(post_pt_root_id)), .groups = "drop")
  
  plot <- ggplot(super_class_count, aes(x=2,y=n_count, fill = int_super_class)) + 
    geom_bar(stat = "identity", width =1) +
    scale_fill_manual(values = super_class_colors) +
    coord_polar(theta = "y") +
    theme_void() + 
    xlim(0.5,2.5) +
    theme(legend.position = "none") +
    annotate("text", x = 0.5, y = 0, label = neuron_count, size = 12, hjust = 0.5,
             vjust = 0.5)
  print(plot)
  
  filepath <- file.path(PATH_output, fig_folder_name, paste0(filename_sensory_int, 
                                                             i, ".pdf"))
  cat("Saving:", filepath, "\n")
  ggsave(filepath, plot = plot, width = 5, height = 5, units = "cm")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Sensory to INTs to NSC- Fig 4 Panel A - Indirect sensory layer
# ------------------------------------------------------------------------------
panel_name = "A_SensorytoINTs_to_"
filename_sensory_to_int = paste0(fig_folder_name, panel_name)

for (i in unique(sensory_interneurons_to_NSC$name_post)) {
  # Get the data for this NSC type (name_post)
  input_df <- unique(sensory_interneurons_to_NSC[sensory_interneurons_to_NSC$name_post == i,])
  class_grp <- input_df %>%
    group_by(class, pre_pt_root_id) %>%
    summarise(n_count = n(), .groups = "drop")
  # count of ids for the plot
  neuron_count <- length(unique(input_df$pre_pt_root_id))
  
  class_count <- class_grp %>%
    group_by(class) %>% 
    summarise(n_count = length(unique(pre_pt_root_id)), .groups = "drop")
  
  plot <- ggplot(class_count, aes(x=2,y=n_count, fill = class)) + 
    geom_bar(stat = "identity", width =1) +
    scale_fill_manual(values = sensory_colors) +
    coord_polar(theta = "y") +
    theme_void() + 
    xlim(0.5,2.5) +
    theme(legend.position = "none") +
    annotate("text", x = 0.5, y = 0, label = neuron_count, size = 12, hjust = 0.5,
             vjust = 0.5)
  print(plot)
  
  filepath <- file.path(PATH_output, fig_folder_name, paste0(filename_sensory_to_int, 
                                                             i, ".pdf"))
  cat("Saving:", filepath, "\n")
  ggsave(filepath, plot = plot, width = 5, height = 5, units = "cm")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Nonsensory INTs to NSC -  - Fig 4 Panel A - Nonsensory interneuron layer
# ------------------------------------------------------------------------------
panel_name = "A_NonsensoryINTs_to_"
filename_nonsensory = paste0(fig_folder_name, panel_name)

for (i in unique(nonsensory_interneurons_to_NSC$name_post)) {
  # Get the data for this NSC type (name_post)
  input_df <- unique(nonsensory_interneurons_to_NSC[nonsensory_interneurons_to_NSC$name_post == i,])
  super_class_grp <- input_df %>%
    group_by(int_super_class, post_pt_root_id) %>%
    summarise(n_count = n(), .groups = "drop")
  # count of ids for the plot
  neuron_count <- length(unique(input_df$post_pt_root_id))
  
  super_class_count <- super_class_grp %>%
    group_by(int_super_class) %>% 
    summarise(n_count = length(unique(post_pt_root_id)), .groups = "drop")
  
  plot <- ggplot(super_class_count, aes(x=2,y=n_count, fill = int_super_class)) + 
    geom_bar(stat = "identity", width =1) +
    scale_fill_manual(values = super_class_colors) +
    coord_polar(theta = "y") +
    theme_void() + 
    xlim(0.5,2.5) +
    theme(legend.position = "none") +
    annotate("text", x = 0.5, y = 0, label = neuron_count, size = 12, hjust = 0.5,
             vjust = 0.5)
  print(plot)
  
  filepath <- file.path(PATH_output, fig_folder_name, paste0(filename_nonsensory, 
                                                             i, ".pdf"))
  cat("Saving:", filepath, "\n")
  ggsave(filepath, plot = plot, width = 5, height = 5, units = "cm")
  }
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Panel 4B: top #s for codex brain plots
# Sensory -> NSC by enteric, mechano, gust
# ------------------------------------------------------------------------------
panel_name = "B_Directsensory_to_NSC_from_"
filename_direct = paste0(fig_folder_name, panel_name)

for (i in unique(sensory_direct_to_NSC$class)) {
  # Get the data for this class
  input_df <- unique(sensory_direct_to_NSC[sensory_direct_to_NSC$class == i,])
  # get the ids
  ids <- as.data.frame(unique(input_df$pre_pt_root_id))
  colnames(ids) <- c("pre_pt_root_id")
  full_df <- filter(interneurons_input_filtered, pre_pt_root_id %in% ids$pre_pt_root_id)
  # now group by pre_id
  top_tmp = full_df %>%
    group_by(pre_pt_root_id)%>%
    summarise(total_syn = sum(n_synapses), .groups = "drop")
  # go through and order
  top10_tmp = top_tmp[order(-top_tmp$total_syn),]$pre_pt_root_id[1:10]
  top_input = full_df[full_df$pre_pt_root_id %in% top10_tmp,]
  # Create a data frame with pre_id and their respective order
  order_df = data.frame(pre_pt_root_id = top10_tmp, order = 1:10)
  # Merge the order information
  top_input = merge(top_input, order_df, by = "pre_pt_root_id")
  # Save to csv
  if (write_csv){
    csv_path <- file.path(PATH_output, fig_folder_name, paste0(filename_direct, 
                                                               i, "_top_10",".csv"))
    write.csv(top_input, csv_path)
    cat("Saved CSV to:", normalizePath(csv_path), "\n")
  }
  
  cat(paste0(unique(order_df$pre_pt_root_id), collapse = ","), "\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Panel 4C: top #s for codex brain plots
# Sensory -> Sensory INT -> NSC by enteric, mechano, gust
# ------------------------------------------------------------------------------
panel_name = "C_Sensoryinput_tosensoryINT_to_NSC_from_"
filename_sensoryINT = paste0(fig_folder_name, panel_name)

for (i in unique(sensory_interneurons_to_NSC$class)) {
  # Get the data for this class
  input_df <- unique(sensory_interneurons_to_NSC[sensory_interneurons_to_NSC$class == i,])
  # get the ids
  ids <- as.data.frame(unique(input_df$pre_pt_root_id))
  colnames(ids) <- c("pre_pt_root_id")
  full_df <- filter(interneurons_input_filtered, pre_pt_root_id %in% ids$pre_pt_root_id)
  # now group by pre_id
  top_tmp = full_df %>%
    group_by(pre_pt_root_id)%>%
    summarise(total_syn = sum(n_synapses), .groups = "drop")
  # go through and order
  top10_tmp = top_tmp[order(-top_tmp$total_syn),]$pre_pt_root_id[1:10]
  top_input = full_df[full_df$pre_pt_root_id %in% top10_tmp,]
  # Create a data frame with pre_id and their respective order
  order_df = data.frame(pre_pt_root_id = top10_tmp, order = 1:10)
  # Merge the order information
  top_input = merge(top_input, order_df, by = "pre_pt_root_id")
  # Save to csv
  if (write_csv){
    csv_path <- file.path(PATH_output, fig_folder_name, paste0(filename_sensoryINT, 
                                                               i, "_top_10",".csv"))
    write.csv(top_input, csv_path)
    cat("Saved CSV to:", normalizePath(csv_path), "\n")
  }
  
  cat(paste0(unique(order_df$pre_pt_root_id), collapse = ","), "\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Panel 4D: indirect inputs to NSC by subclass
# ------------------------------------------------------------------------------
panel_name = "D_sensory_grp_table"
filename_ind = paste0(fig_folder_name, panel_name)

sensory_grp <- sensory_interneurons_to_NSC %>%
  group_by(class, sub_class) %>%
  summarise(n_neurons = n_distinct(pre_pt_root_id), .groups = "drop") %>%
  arrange(class, desc(n_neurons))

# Save to csv
if (write_csv){
  csv_path <- file.path(PATH_output, fig_folder_name, paste0(filename_ind, ".csv"))
  write.csv(sensory_grp, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}
# ------------------------------------------------------------------------------
