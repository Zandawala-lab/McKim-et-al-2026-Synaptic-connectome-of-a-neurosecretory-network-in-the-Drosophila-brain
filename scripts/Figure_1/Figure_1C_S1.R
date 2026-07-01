# --- Figure 1C and S1 ---------------------------------------------------------
# ------------------------------------------------------------------------------
# 1C: Code for NSC clustering analysis in FAFB v783 & 
# S1: Plotting NSC neuron reconstructions in BANC v888 and maleCNS v0.9
# ------------------------------------------------------------------------------

# load packages ----------------------------------------------------------------
library(tidyverse)
library(coconatfly)
library(bancr)
library(fafbseg)
library(natverse)
library(malecns)
library(neuprintr)
library(dplyr)
library(rgl)

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

# Read version file (provided) & get # to use in output filenames - v783
v = read_delim(file.path(PATH_input,"version.csv"),
               col_types  =  cols(version  =  col_character()),delim  =  ";")
v = v$version[1]

# Folder name for saving figs + check if exists
fig_folder_name = "Figure_1"
if(!dir.exists(file.path(PATH_output, fig_folder_name))){
  dir.create(file.path(PATH_output, fig_folder_name))
  cat("Created figure folder at:", normalizePath(file.path(PATH_output, fig_folder_name)), "\n")
}

# ------------------------------------------------------------------------------
# Options for saving plots & csv files
# ------------------------------------------------------------------------------

# When running script for the first time, first two to TRUE:
write_plots = TRUE           # TRUE - save/replicate figure plots
                             # FALSE - plots not saved outside of R
write_csv = TRUE             # TRUE - save processed data associated w/figures
                             # FALSE - data not saved outside of R
write_mesh = FALSE           # FALSE - do not generate mesh plots by default (takes awhile to render)

# set colors -------------------------------------------------------------------
m_NSC_unknown  = "#BFD739"
l_NSC_unknown  = "#009817"
m_NSC_DILP     = "#BD0023"
m_NSC_DH44     = "#ffe200"
m_NSC_DMS      = "#FE7E00"
l_NSC_CRZ      = "#8100FF"
l_NSC_ITP      = "#838383"
l_NSC_DH31     = "#00B4FF"
SEZ_NSC_CAPA   = "#003BBD"
SEZ_NSC_Hugin  = "#BD00B0"
l_NSC_unknown2 = "black" # maleCNS and BANC

syn_col <- c("m_NSC_unknown" = m_NSC_unknown,"l_NSC_unknown" = l_NSC_unknown,
             "m_NSC_DILP" = m_NSC_DILP,"m_NSC_DH44" = m_NSC_DH44,
             "m_NSC_DMS" = m_NSC_DMS,"l_NSC_CRZ" = l_NSC_CRZ,
             "l_NSC_ITP" = l_NSC_ITP,"l_NSC_DH31" = l_NSC_DH31,
             "SEZ_NSC_CAPA" = SEZ_NSC_CAPA,"SEZ_NSC_Hugin" = SEZ_NSC_Hugin,
             "l_NSC_unknown2" = l_NSC_unknown2)

# load data --------------------------------------------------------------------
# Download from zenodo: 10.5281/zenodo.20795211
# Provided w/paper: McKim et al. NSC connectome: https://doi.org/10.7554/eLife.102684
# Supplementary Table 3 - has labels and IDs across datasets

NSC <- read_delim(file.path(PATH_input, paste0("NSC_v", v, ".csv")),
                  delim = ",",
                  escape_double = FALSE,
                  col_types = cols(NSC_id = col_character()),
                  trim_ws = TRUE)

# ------------------------------------------------------------------------------
# Figure 1C - FAFB - Clustering of endocrine cells
# ------------------------------------------------------------------------------
panel_name = "C"
filename_clustering = paste0(fig_folder_name, panel_name)

# may need to set the version
fafbseg::flywire_connectome_data_version(set = as.numeric(v)) #783

# Data prep
colnames(NSC)[2] <- 'id'

# Get metadata
endocrine <- cf_meta(cf_ids('/super_class:endocrine', datasets = c("flywire")))

# Join dfs
endocrine = left_join(endocrine,NSC,by = "id")
endocrine$hemisphere_name_id = paste(
  endocrine$hemisphere,
  endocrine$NSC_name,endocrine$id,sep="_")

# Plot - for review 
par(mar = c(10, 10, 4, 2) + 0.1)
endocrine %>%
  with(cf_cosine_plot(key, threshold = 1, labRow = hemisphere_name_id, interactive = FALSE, partners = "inputs", method = "ward.D2"))

# Plot - save to pdf
if (write_plots) {
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename_clustering, ".pdf"))
  pdf(pdf_path, width = 20, height = 16)
  par(mar = c(10, 10, 4, 2) + 0.1)
  endocrine %>%
    with(
      cf_cosine_plot(
        key,
        threshold = 1,
        labRow = hemisphere_name_id,
        interactive = FALSE,
        partners = "inputs",
        method = "ward.D2"
      )
    )
  dev.off()
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}

# ------------------------------------------------------------------------------
#  Figure_1 S1A - BANC - Plot neuron meshes of endocrine cells
# ------------------------------------------------------------------------------

# Download BANC meta data v888 -------------------------------------------------
# all BANC meta for v888
banc_meta <- banc_codex_annotations(source = c("gcs"))

# pt_root_id comes in as character
#class(banc_meta$pt_root_id)
banc_meta <- banc_meta %>%
  mutate(pt_root_id = as.character(pt_root_id))

# subset all endocrine cells in brain and VNC
BANC_endocrine <- subset(
  banc_meta,
  super_class %in% c(
    "visceral_circulatory"))

# subset only the brain endocrine cells
BANC_brain_NSC <- subset(
  BANC_endocrine,
  cell_class %in% c(
    "subesophageal_zone_neurosecretory_cell",
    "pars_intercerebralis_neurosecretory_cell",
    "pars_lateralis_neurosecretory_cell"))

# remove additional endocrine cells not in current paper (McKim et al)
BANC_brain_NSC <- BANC_brain_NSC %>%
  filter(!cell_type %in% c("CB2748", "DNg28") & !is.na(cell_type))

# one of the nine ITP neurons is not ITP - remove 
BANC_brain_NSC <- BANC_brain_NSC %>%
  filter(pt_root_id != "720575941611801674")

# Table 1 - counts per NSC class
BANC_endocrine_count <- BANC_brain_NSC %>%
  group_by(cell_type) %>%
  summarise(count = length(unique(pt_root_id)))

# save to csv
if (write_csv) {
  csv_path <- file.path(PATH_output, fig_folder_name, 
                        "Table_1_BANC_endocrine_count.csv")
  write.csv(BANC_endocrine_count, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# save to csv (used in Figure 3 S3 & 7 S3; Supp Table 3)
if (write_csv) {
  csv_path <- file.path(PATH_input, "BANC_NSC_brain_meta_v888.csv")
  write.csv(BANC_brain_NSC, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}
# ------------------------------------------------------------------------------

if (write_mesh) {
  # Neuron plots ---------------------------------------------------------------
  panel_name = "_S1A_"
  panel_type = "BANC_meshes_all_endocrine"
  filename_banc = paste0(fig_folder_name, panel_name, panel_type)
  
  # Iterate over unique endocrine cell types, download meshes and plot all meshes
  open3d()
  par3d(windowRect = c(0, 0, 3840, 2160))
  view3d(userMatrix = rotationMatrix(pi, 0, 1, 1) %*% rotationMatrix(pi, 0, 0, 1))
  
  for (i in unique(BANC_brain_NSC$cell_type)) {
    # Get the bodyids for the current cell type
    endocrine_ids <- BANC_brain_NSC %>% filter(cell_type == i) %>% pull(root_id)
    # Download meshes
    endocrine_meshes <- banc_read_neuron_meshes(endocrine_ids)
    # Plot meshes
    plot3d(
      endocrine_meshes,
      col = syn_col[i],
      add = TRUE,
      WithNodes = F
    )
  }
  # Add the surface model and zoom in
  plot3d(banc_brain_neuropil.surf,
         col = "grey",
         alpha = 0.1) #lightgrey
  par3d(zoom = 0.53)
  
  if (write_plots) {
    # Export as png
    png_filename <- file.path(PATH_output,
                              fig_folder_name,
                              paste0(filename_banc, ".png"))
    rgl.snapshot(png_filename)
    cat("Saved PNG to:", normalizePath(png_filename), "\n")
  }
  # Close the 3D plot
  close3d()
}
# ------------------------------------------------------------------------------
# Figure 1 S1B - maleCNS - Plot neuron meshes of endocrine cells
# ------------------------------------------------------------------------------

# Download maleCNS meta data v0.9 ----------------------------------------------

## Set your Neuprint token: https://natverse.org/neuprintr/#authentication

# choose dataset
neuprint_dataset <- "male-cns:v0.9"

# get meta data for endocrine cells and plot -----------------------------------
endocrine_meta <- mcns_neuprint_meta("/superclass:cb_endocrine")

# add a new column with McKim et al nomenclature
endocrine_meta <- endocrine_meta %>%
  mutate(
    MckimType = case_when(
      flywireType == "IPC" ~ "m_NSC_DILP",
      flywireType == "DMS" ~ "m_NSC_DMS",
      flywireType == "DH44" ~ "m_NSC_DH44",
      flywireType == "mNSC_unknown" ~ "m_NSC_unknown",
      flywireType == "CAPA" ~ "SEZ_NSC_CAPA",
      flywireType == "Hugin-RG" ~ "SEZ_NSC_Hugin",
      flywireType == "ITP" ~ "l_NSC_ITP",
      flywireType == "CRZ" ~ "l_NSC_CRZ",
      flywireType == "DNES1" ~ "l_NSC_unknown",
      flywireType == "DNES2" ~ "l_NSC_DH31",
      flywireType == "DNES3" ~ "l_NSC_DH31",
      is.na(flywireType) ~ "l_NSC_unknown2",
      TRUE ~ flywireType))

# one of l_NSC_unknown2 is l_NSC_CRZ - reclassify
endocrine_meta <- endocrine_meta %>%
  mutate(
    MckimType = ifelse(
      bodyid == 35813 & MckimType == "l_NSC_unknown2",
      "l_NSC_CRZ",
      MckimType))

maleCNS_brain_NSC <- endocrine_meta %>%
  filter(!(bodyid == 134996 & MckimType == "l_NSC_unknown2"))

# Table 1 - counts per NSC class
endocrine_count <- maleCNS_brain_NSC %>%
  group_by(MckimType) %>%
  summarise(count = length(unique(bodyid)))

# save to csv
if (write_csv) {
  csv_path <- file.path(PATH_output, fig_folder_name, 
                     "Table_1_maleCNS_endocrine_count.csv")
  write.csv(endocrine_count, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# save to csv (used in Figure 3 S3 & 7 S3; Supp Table 3)
if (write_csv) {
  csv_path <- file.path(PATH_input, "maleCNS_brain_NSC_meta_v0.9.csv")
  write.csv(maleCNS_brain_NSC, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# ------------------------------------------------------------------------------

if (write_mesh) {
  # Neuron plots ---------------------------------------------------------------
  panel_name = "_S1B_"
  panel_type = "maleCNS_meshes_all_endocrine"
  filename_mcns = paste0(fig_folder_name, panel_name, panel_type)
  
  # Iterate over unique endocrine cell types, download meshes and plot all meshes
  open3d()
  # Set high resolution for the plot
  par3d(windowRect = c(0, 0, 3840, 2160))
  
  for (i in unique(endocrine_meta$MckimType)) {
    # Get the bodyids for the current cell type
    endocrine_ids <- endocrine_meta %>% filter(MckimType == i) %>% pull(bodyid)
    # Download meshes
    endocrine_meshes <- read_mcns_meshes(endocrine_ids)
    # plot meshes
    plot3d(
      endocrine_meshes,
      col = syn_col[i],
      add = TRUE,
      WithNodes = F
    )
  }
  
  # Add the surface model
  plot3d(
    malecns.surf,
    alpha = 0.1,
    add = TRUE,
    col = "lightgrey"
  )
  # Adjust the view
  view3d(userMatrix = rotationMatrix(90 * pi / 90, 1, 0, 0), zoom = 0.5)
  
  if (write_plots) {
    # Export as png
    png_filename <- file.path(PATH_output,
                              fig_folder_name,
                              paste0(filename_mcns, ".png"))
    rgl.snapshot(png_filename)
    cat("Saved PNG to:", normalizePath(png_filename), "\n")
  }
  # Close the 3D plot
  close3d()
}
