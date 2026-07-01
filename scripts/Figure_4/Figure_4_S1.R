# -- Figure 4 S1 ---------------------------------------------------------------
# ------------------------------------------------------------------------------
# Code for A & B
# ------------------------------------------------------------------------------

# load packages ----------------------------------------------------------------
library(tidyverse)
library(coconatfly)
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

# load data --------------------------------------------------------------------
# download flywire dataset - uncomment if next line doesn't run
# fafbseg::download_flywire_release_data()

# obtain meta data for all unknown sensory
unknown_sensory <- cf_meta(cf_ids(flywire='/cell_class:unknown_sensory'))

# get the enterics only
enteric_groups <- unknown_sensory %>%
  filter(grepl("^ENS[1-5]$", type))

panel_name = "_S1A_enteric_groups"
filename_ent = paste0(fig_folder_name, panel_name)

# save csv
if (write_csv){
  csv_path <- file.path(PATH_output, fig_folder_name, paste0(filename_ent, ".csv"))
  write.csv(enteric_groups, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# ------------------------------------------------------------------------------
# Panel A - IDs for neuroglancer
# ------------------------------------------------------------------------------
# Loop through ENS groups 1-5, print id list and count for each
for (g in paste0("ENS", 1:5)) {
  group_data <- enteric_groups %>% filter(type == g)
  
  cat("\n", g, "- count:", nrow(group_data), "\n")
  cat(paste(group_data$id, collapse = ", "), "\n")
}
# IDs were then copy and pasted into neuroglancer for visualization
# See associated figure .rtf files with links
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Panel B - Cosine similarity for ENS groups
# ------------------------------------------------------------------------------
panel_name = "_S1B_ENS_cosine_plot.pdf"
filename_ent_sim = paste0(fig_folder_name, panel_name)

# Cosine similarity plot
outpath <- file.path(PATH_output, fig_folder_name, filename_ent_sim)
pdf(outpath, width = 10, height = 10)
cf_cosine_plot(
  cf_ids('/type:ENS.+', datasets = c("flywire")), 
  threshold = 1, 
  interactive = FALSE, 
  partners = c("outputs", "inputs"), 
  method = "ward.D2"
)
dev.off()

cat("Saved to:", outpath, "\n")

# ------------------------------------------------------------------------------

