# -- Figure 7 - S2 -------------------------------------------------------------
# ------------------------------------------------------------------------------
# Code for plotting output of NSCs at low threshold (>= 2 synapses)
# ------------------------------------------------------------------------------

# load packages ----------------------------------------------------------------
library(tidyverse)
library(fafbseg)
library(ggplot2)
library(natverse)
library(ggtext)

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
                        "undefined"          = "#999999", # was white
                        "optic"              = "white",   # was gray
                        "descending"         = "#c3905f",  
                        "visual_centrifugal" = "#f46d43",
                        "central"            = "#3384b8", 
                        "visual_projection"  = "#c8d684", 
                        "sensory"            = "#b73545", 
                        "endocrine"          = "#a8c8c0", # was grey
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
synapses_output =  read_csv(file.path(PATH_input, paste0("tmp/NSC_output_filtered_v",v,".csv")),
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
# Organize connectivity data for all figures
# ------------------------------------------------------------------------------

# Summarize synapses data
NSC_output_sum <- synapses_output %>%
  group_by(pre_pt_root_id, post_pt_root_id) %>%
  summarise(n_synapses = n(),
            .groups = "drop")

# Large classification table with neurotransmitters and other annotations
classification <- left_join(classification, neurotransmitters, by = "root_id")
# Reclassify neurotransmitters based on stringent cutoff
classification$nt_type_stringent <- classification$nt_type
classification <- classification %>%
  mutate(nt_type_stringent = ifelse(nt_type_score <= 0.62, "uncertain",
                                    nt_type_stringent))

# Join data for further analysis
classification_join <- classification
colnames(classification_join) <- c("post_pt_root_id", colnames(classification)[-1])
# Filter and keep >= 2 synapses
NSC_output <- left_join(NSC_output_sum[NSC_output_sum$n_synapses >= 2, ], 
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
# plots & analyses of output
# ------------------------------------------------------------------------------
if (write_csv){
  csv_path <- file.path(PATH_input, paste0("NSC_output_classification_2thresh_v",
                                           v, ".csv"))
  write.csv(NSC_output, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}
# ------------------------------------------------------------------------------

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
                                                          "endocrine",
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
# Figure 7 S2A - Proportion output synapses from NSCs by super class
# ------------------------------------------------------------------------------
panel_name = "_S2A"
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
  ggsave(pdf_path, plot = p, width = 12, height = 10, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 7 S2B - NSC #'s & super class
# ------------------------------------------------------------------------------
cell_count = NSC_output%>%
  group_by(super_class)%>%
  summarize(n_NSC = length(unique(pre_pt_root_id)),
            n_superclass = length(unique(post_pt_root_id))) %>% arrange(desc(n_NSC))

# Save data to csv
if (write_csv){
  csv_path <- file.path(PATH_output, fig_folder_name, 
                        paste0("Figure_7_S2B_numNeurons_outputNSC_v", v, ".csv"))
  write.csv(cell_count, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# ------------------------------------------------------------------------------
# Figure 7 S2C  - Output IDs from NSCs by super class
# ------------------------------------------------------------------------------
# IDs from csv were put into neuroglancer for visualization

col_to_keep <- c("pre_pt_root_id","cell_type","class","hemibrain_type","name_pre",
               "name_post","hemisphere_post","post_pt_root_id","n_synapses", 
               "nt_type_score", "nt_type", "nt_type_stringent", "nt_type_stringent_corr")

ids_to_print <- c("central", "endocrine", "descending")  # NA/undefined handled separately

for (i in unique(NSC_output$super_class)) {
  if (is.na(i)) {
    NSC_subset_superclass <- unique(NSC_output[is.na(NSC_output$super_class), col_to_keep])
    filename_part <- "NA"
  } else {
    NSC_subset_superclass <- unique(NSC_output[NSC_output$super_class == i, col_to_keep])
    filename_part <- i
  }
  
  filename_part <- gsub(" ", "_", filename_part)
  
  # Name and save file
  filename <- paste(fig_folder_name, 
                    "/Figure_7_S2C_NSC_output_filtered_2_syn_threshold_",
                    i, sep = "")
  outpath <- file.path(PATH_output, paste0(filename,"_v", v,".csv"))
  write.csv(NSC_subset_superclass, outpath)
  cat("Saved CSV to:", normalizePath(outpath), "\n")
  
  # Print unique IDs only for the 4 categories of interest
  if (is.na(i) || filename_part %in% ids_to_print) {
    unique_ids <- unique(NSC_subset_superclass$post_pt_root_id)
    unique_ids <- unique_ids[!is.na(unique_ids)]
    cat(filename_part, "- count:", length(unique_ids), "\n")
    cat(paste(unique_ids, collapse = ", "), "\n")
  }
}

# post_pt_root_ids based on super class are put into
# neuroglancer (https://edit.flywire.ai/) for visualization & screenshots
# neuroglancer links provided as rtf's on zenodo
#-------------------------------------------------------------------------------
