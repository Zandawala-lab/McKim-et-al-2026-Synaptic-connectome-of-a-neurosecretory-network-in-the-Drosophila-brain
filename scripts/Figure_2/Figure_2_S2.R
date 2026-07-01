# -- Figure 2 S2 ---------------------------------------------------------------
# ------------------------------------------------------------------------------
# Code for plotting:
# B,C,D inputs to DH44 and DMS NSC in maleCNS &
# H,I,J mean gray value for DCVs across datasets
# ------------------------------------------------------------------------------

# load packages ----------------------------------------------------------------
library(tidyverse)
library(fafbseg)
library(natverse)
library(malecns)
library(neuprintr)
library(dplyr)
library(ggplot2)
library(multcompView)
library(tidyplots)
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

# Folder name for saving figs + check if exists
fig_folder_name = "Figure_2"
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
write_mesh = FALSE           # FALSE - do not generate mesh plots by default (takes awhile to render)


# set colors -------------------------------------------------------------------
m_NSC_DH44 = "#ffe200"
m_NSC_DMS = "#FE7E00"
DMS = "darkgrey"

plot_colors <- c("m_NSC_DH44" = m_NSC_DH44,"m_NSC_DMS" = m_NSC_DMS,"DMS" = DMS)

# load and filter data ---------------------------------------------------------

# ------------------
# maleCNS plots (B-D)
# ------------------
# this file was generated in Figure_1C_S1 script (S1B part)
maleCNS_brain_NSC <- read_delim(file.path(PATH_input, "maleCNS_brain_NSC_meta_v0.9.csv"),
                                delim = ",",
                                escape_double = FALSE,
                                col_types = cols(bodyid = col_character()),
                                trim_ws = TRUE)

maleCNS_DH44 <- maleCNS_brain_NSC %>%
  filter(MckimType == "m_NSC_DH44")

maleCNS_DMS <- maleCNS_brain_NSC %>%
  filter(MckimType == "m_NSC_DMS")

# choose dataset
neuprint_dataset <- "male-cns:v0.9"

# download DH44 inputs
maleCNS_DH44_input <- mcns_connection_table(maleCNS_DH44$bodyid, 
                                            partners='inputs', 
                                            moredetails = TRUE, 
                                            threshold = 1,
                                            roi = NULL,
                                            by.roi = FALSE)

# filter
maleCNS_DH44_5thresh_input <- maleCNS_DH44_input %>%
  filter(weight >= 5)

length(unique(maleCNS_DH44_5thresh_input))

# download DMS inputs
maleCNS_DMS_input <- mcns_connection_table(maleCNS_DMS$bodyid, 
                                           partners='inputs', 
                                           moredetails = TRUE, 
                                           threshold = 1,
                                           roi = NULL,
                                           by.roi = FALSE)

# filter
maleCNS_DMS_5thresh_input <- maleCNS_DMS_input %>%
  filter(weight >= 5)

length(unique(maleCNS_DMS_5thresh_input))

# ------------------
# DCV plots (H-J)
# ------------------
FAFB_DCV = read_delim(file.path(PATH_input,"Figure2_S2_FAFB_DCV.csv"),
                      col_types = cols(id = col_character(), 
                                       Cell_type = col_character()), delim = ",",)

BANC_DCV <- read_delim(file.path(PATH_input, "Figure2_S2_BANC_DCV.csv"),
                       col_types = cols(id = col_character(), 
                                        Cell_type = col_character()), delim = ",")

maleCNS_DCV <- read_delim(file.path(PATH_input, "Figure2_S2_maleCNS_DCV.csv"),
                          col_types = cols(id = col_character(), 
                                           Cell_type = col_character()), delim = ",")

# ------------------------------------------------------------------------------
#                                     maleCNS
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Figure 2 S2B - maleCNS - DH44 & inputs
# ------------------------------------------------------------------------------
if (write_mesh) {
  panel_name = "_S2B_"
  panel_type = "maleCNS_DH44_5thresh_input"
  filename_mcns_dh44 = paste0(fig_folder_name, panel_name, panel_type)
  
  # download meshes and plot all meshes together
  open3d()
  # Set high resolution for the plot
  par3d(windowRect = c(0, 0, 3840, 2160))
  # plot DH44 neurons
  for (i in unique(maleCNS_DH44$bodyid)) {
    # Download meshes
    DH44_meshes <- read_mcns_meshes(i)
    # plot meshes
    plot3d(
      DH44_meshes,
      col = "#ffe200",
      add = TRUE,
      WithNodes = F
    )
  }
  # plot DH44 inputs
  for (i in unique(maleCNS_DH44_5thresh_input$partner)) {
    # Download meshes
    input_meshes <- read_mcns_meshes(i)
    # plot meshes
    plot3d(
      input_meshes,
      col = "magenta",
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
                              paste0(filename_mcns_dh44, ".png"))
    rgl.snapshot(png_filename)
    cat("Saved PNG to:", normalizePath(png_filename), "\n")
  }
  # Close the 3D plot
  close3d()
  # ----------------------------------------------------------------------------
  
  # ----------------------------------------------------------------------------
  # Figure 2 S2C - maleCNS - DMS NSC & inputs
  # ----------------------------------------------------------------------------
  panel_name = "_S2C_"
  panel_type = "maleCNS_DMS_NSC_5thresh_input"
  filename_mcns_dms_nsc = paste0(fig_folder_name, panel_name, panel_type)
  
  # download meshes and plot all meshes together
  open3d()
  # Set high resolution for the plot
  par3d(windowRect = c(0, 0, 3840, 2160))
  # plot DH44 neurons
  for (i in unique(maleCNS_DMS$bodyid)) {
    # Download meshes
    DMS_meshes <- read_mcns_meshes(i)
    # plot meshes
    plot3d(
      DMS_meshes,
      col = "#FE7E00",
      add = TRUE,
      WithNodes = F
    )
  }
  # plot DH44 inputs
  for (i in unique(maleCNS_DMS_5thresh_input$partner)) {
    # Download meshes
    input_meshes <- read_mcns_meshes(i)
    # plot meshes
    plot3d(
      input_meshes,
      col = "magenta",
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
                              paste0(filename_mcns_dms_nsc, ".png"))
    rgl.snapshot(png_filename)
    cat("Saved PNG to:", normalizePath(png_filename), "\n")
  }
  # Close the 3D plot
  close3d()
  # ----------------------------------------------------------------------------
  
  # ----------------------------------------------------------------------------
  # Figure 2 S2D - maleCNS - DMS DNp32 & inputs
  # ----------------------------------------------------------------------------
  panel_name = "_S2D_"
  panel_type = "maleCNS_DMS_DNp32_5thresh_input"
  filename_mcns_dms_dnp = paste0(fig_folder_name, panel_name, panel_type)
  
  DMS_DNp32ids <- c("524968", "556286")
  
  # download meshes and plot all meshes together
  open3d()
  # Set high resolution for the plot
  par3d(windowRect = c(0, 0, 3840, 2160))
  # plot DH44 neurons
  for (i in unique(DMS_DNp32ids)) {
    # Download meshes
    DMS_meshes <- read_mcns_meshes(i)
    # plot meshes
    plot3d(
      DMS_meshes,
      col = "#c3905f",
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
                              paste0(filename_mcns_dms_dnp, ".png"))
    rgl.snapshot(png_filename)
    cat("Saved PNG to:", normalizePath(png_filename), "\n")
  }
  # Close the 3D plot
  close3d()
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#                     Figure 2 S2 HIJ - DCVs across datasets
# ------------------------------------------------------------------------------
panel_name = "_S2HIJ"
filename_dcv = paste0(fig_folder_name, panel_name)

# Add dataset labels
FAFB_DCV$dataset    <- "FAFB"
maleCNS_DCV$dataset <- "maleCNS"
BANC_DCV$dataset    <- "BANC"

# Combine datasets
combined_DCV <- bind_rows(FAFB_DCV, maleCNS_DCV, BANC_DCV)

# Reset factor levels with new labels
combined_DCV$Cell_type <- factor(combined_DCV$Cell_type,
                                 levels = rev(c("m_NSC_DH44", "m_NSC_DMS", "DMS")))

combined_DCV$dataset <- factor(combined_DCV$dataset,
                               levels = c("FAFB", "maleCNS", "BANC"))

# Create each panel by filtering from combined_DCV
p1 <- combined_DCV %>%
  filter(dataset == "FAFB") %>%
  tidyplot(x = Cell_type, y = Mean, color = Cell_type) %>%
  add_boxplot() %>%
  add_data_points_beeswarm() %>%
  add_test_pvalue(method = "tukey_hsd", step.increase = 0.25) %>%
  adjust_y_axis_title("Mean gray value") %>%
  adjust_y_axis(limits = c(20, 170), breaks = seq(25, 150, by = 25)) %>%
  adjust_colors(new_colors = plot_colors) %>%
  adjust_title("FAFB") %>%
  adjust_legend_title("Cell type") %>%
  adjust_legend_position("none") %>%
  remove_x_axis_labels() %>%
  remove_x_axis_title() %>%
  remove_x_axis_ticks() %>%
  remove_caption() +
  scale_x_discrete(labels = scales::label_parse()) +
  theme(axis.text.x = element_blank(),
        legend.text = element_text(size = 8),
        plot.margin = margin(t = 20, r = 5, b = 5, l = 5, unit = "mm"))

p2 <- combined_DCV %>%
  filter(dataset == "maleCNS") %>%
  tidyplot(x = Cell_type, y = Mean, color = Cell_type) %>%
  add_boxplot() %>%
  add_data_points_beeswarm() %>%
  add_test_pvalue(method = "tukey_hsd", step.increase = 0.25) %>%
  adjust_y_axis_title("Mean gray value") %>%
  adjust_y_axis(limits = c(20, 170), breaks = seq(25, 150, by = 25)) %>%
  adjust_colors(new_colors = plot_colors) %>%
  adjust_title("maleCNS") %>%
  adjust_legend_title("Cell type") %>%
  adjust_legend_position("none") %>%
  remove_x_axis_labels() %>%
  remove_x_axis_title() %>%
  remove_x_axis_ticks() %>%
  remove_caption() +
  scale_x_discrete(labels = scales::label_parse()) +
  theme(axis.text.x = element_blank(),
        legend.text = element_text(size = 8),
        plot.margin = margin(t = 20, r = 5, b = 5, l = 5, unit = "mm"))

p3 <- combined_DCV %>%
  filter(dataset == "BANC") %>%
  tidyplot(x = Cell_type, y = Mean, color = Cell_type) %>%
  add_boxplot() %>%
  add_data_points_beeswarm() %>%
  add_test_pvalue(method = "tukey_hsd", step.increase = 0.25) %>%
  adjust_y_axis_title("Mean gray value") %>%
  adjust_y_axis(limits = c(20, 170), breaks = seq(25, 150, by = 25)) %>%
  adjust_colors(new_colors = plot_colors) %>%
  adjust_title("BANC") %>%
  adjust_legend_title("Cell type") %>%
  adjust_legend_position("right") %>%
  remove_x_axis_labels() %>%
  remove_x_axis_title() %>%
  remove_x_axis_ticks() %>%
  remove_caption() +
  scale_x_discrete(labels = scales::label_parse()) +
  theme(axis.text.x = element_blank(),
        legend.text = element_text(size = 8),
        plot.margin = margin(t = 20, r = 5, b = 5, l = 5, unit = "mm"))

# Combine with patchwork
plot_combined <- p1 | p2 | p3

print(plot_combined)

if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename_dcv, ".pdf"))
  ggsave(pdf_path, plot = plot_combined, width = 24, height = 15, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}
