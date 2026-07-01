# -- Figure 3 S1 ---------------------------------------------------------------
# ------------------------------------------------------------------------------
# Code for plotting input synapses of all NSCs with skeletons 
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
fig_folder_name = "Figure_3"
if(!dir.exists(file.path(PATH_output, fig_folder_name))){
  dir.create(file.path(PATH_output, fig_folder_name))
  cat("Created figure folder at:", normalizePath(file.path(PATH_output, fig_folder_name)), "\n")
}

# ------------------------------------------------------------------------------
# Options for saving plots & csv files
# ------------------------------------------------------------------------------
# When running script for the first time, set both to TRUE:
write_plots = TRUE           # TRUE - save/replicate figure plots
                             # FALSE - plots not saved outside of R

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
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Run python scripts - open .py in R (requires reticulate) or run elsewhere
# otherwise, download files from zenodo and put into: output/neurons_v783/
# ------------------------------------------------------------------------------
# Save brainmesh.obj for 3D plotting -------------------------------------------
# Figure_3_S1_brainmesh.py

# Save NSC skeleton objs for 3D plotting ---------------------------------------
# Figure_3_S1_neuronMesh.py

# ------------------------------------------------------------------------------
# Figure 3 S1 - Plot input synapses of all NSCs
# ------------------------------------------------------------------------------
panel_name = "_S1_"
filename_input = paste0(fig_folder_name, panel_name)

synapses_input$pre_x = xyzmatrix(synapses_input$pre_pt_position)[,1]
synapses_input$pre_y = xyzmatrix(synapses_input$pre_pt_position)[,2]
synapses_input$pre_z = xyzmatrix(synapses_input$pre_pt_position)[,3]
synapses_input$post_x = xyzmatrix(synapses_input$post_pt_position)[,1]
synapses_input$post_y = xyzmatrix(synapses_input$post_pt_position)[,2]
synapses_input$post_z = xyzmatrix(synapses_input$post_pt_position)[,3]

# Iterate over unique NSC names and plot synapses and meshes
for (i in unique(NSC$NSC_name)) {
  # Get the NSC_ids for the current NSC_name
  NSC_ids <- NSC %>% filter(NSC_name == i) %>% pull(NSC_id)
  # Filter the synapses for the current NSC_ids
  synapses_filtered <- synapses_input %>% filter(post_pt_root_id %in% NSC_ids)
  # Get skeletons names for NSCs
  NSC_ids_skel <- paste0(NSC_ids, ".obj")
  # Plot skeletons and input synapses
  open3d()
  # Set high resolution for the plot
  par3d(windowRect = c(0, 0, 3840, 2160))
  for (j in NSC_ids_skel) {
    NSC_skeletons_tmp <- readOBJ(paste0(PATH_output, "/neurons_v783/NSC_mesh/",j))
    plot3d(NSC_skeletons_tmp, col = "black", add= T, WithNodes = F)
  } 
  plot3d(synapses_filtered$post_x, synapses_filtered$post_y, synapses_filtered$post_z,
         col = syn_col[i], size = 0.5, type = "s", add = TRUE)
  # Add the surface model
  brainmesh <- readOBJ(paste0(PATH_output, "/neurons_v783/brainmesh/brainmesh.obj"))
  plot3d(brainmesh, add = TRUE, alpha = 0.1, col = "grey")
  # Adjust the view
  view3d(userMatrix = rotationMatrix(90 * pi / 90, 1, 0, 0), zoom = 0.5)  
  if(write_plots) {
    # Export as png
    png_filename <- file.path(PATH_output, fig_folder_name, paste0(filename_input, i, ".png"))
    rgl.snapshot(png_filename)
    cat("Saved PNG to:", normalizePath(png_filename), "\n")
  }
  # Close the 3D plot
  close3d()
}


