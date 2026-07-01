# -- Figure 3 - S3 -------------------------------------------------------------
# ------------------------------------------------------------------------------
# Code for visualizing input to NSCs in BANC v888 & maleCNS v0.9
# ------------------------------------------------------------------------------

# load packages ----------------------------------------------------------------
library(tidyverse)
library(bancr)
library(ggplot2)
library(malecns)
library(neuprintr)
library(rgl)
library(ggtext)

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
fig_folder_name = "Figure_3"
if(!dir.exists(file.path(PATH_output, fig_folder_name))){
  dir.create(file.path(PATH_output, fig_folder_name))
  cat("Created figure folder at:", normalizePath(file.path(PATH_output,
                                                           fig_folder_name)), "\n")
}
# ------------------------------------------------------------------------------
# Options for saving plots & csv files
# ------------------------------------------------------------------------------

# When running script for the first time, use provided settings:
write_plots = TRUE           # TRUE - save/replicate figure plots
                             # FALSE - plots not saved outside of R
write_csv = TRUE             # TRUE - save processed data associated w/figures
                             # FALSE - data not saved outside of R
load_data = TRUE             # TRUE - use data generated from Figure 1_S1
                             # FALSE - start from scratch in this script

# set colors: ------------------------------------------------------------------
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

super_class_colors <- c("others"                  = "#5f50a1", 
                        "undefined"               = "#999999", 
                        "optic"                   = "white",
                        "descending"              = "#c3905f",  
                        "visual_centrifugal"      = "#f46d43",
                        "central_brain_intrinsic" = "#3384b8", 
                        "visual_projection"       = "#c8d684", 
                        "sensory"                 = "#b73545",
                        "sensory_ascending"       = "#d6808b",
                        "ascending"               = "#a9d5a3",
                        "visceral_circulatory"    = "#a8c8c0", # BANC
                        "endocrine"               = "#a8c8c0", # maleCNS
                        "motor"                   = "black")

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

# all BANC meta for v888
banc_meta <- banc_codex_annotations(source = c("gcs"))

# pt_root_id comes in as character
#class(banc_meta$pt_root_id)
banc_meta <- banc_meta %>%
  mutate(pt_root_id = as.character(pt_root_id))

# Save csv (used in Figure 7 S3 & Supp Table 3)
if (write_csv){
  csv_path <- file.path(PATH_input, paste0("BANC_meta_v888.csv"))
  write.csv(banc_meta, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# Filtered NSC data for BANC and maleCNS (from Figure 1C_S1)
if (load_data){
  BANC_brain_NSC <- read_delim(file.path(PATH_input, "BANC_NSC_brain_meta_v888.csv"),
                               delim = ",",
                               escape_double = FALSE,
                               col_types = cols(pt_root_id = col_character()),
                               trim_ws = TRUE)
  
  maleCNS_brain_NSC <- read_delim(file.path(PATH_input, "maleCNS_brain_NSC_meta_v0.9.csv"),
                                  delim = ",",
                                  escape_double = FALSE,
                                  col_types = cols(bodyid = col_character()),
                                  trim_ws = TRUE)
  
  
} else {
  
  # BANC --------------------------------------
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
  
  # save to csv (used in Figure 3 S3 & 7 S3; Supp Table 3)
  if (write_csv) {
    csv_path <- file.path(PATH_input, "BANC_NSC_brain_meta_v888.csv")
    write.csv(BANC_brain_NSC, csv_path)
    cat("Saved CSV to:", normalizePath(csv_path), "\n")
  }
  
  # maleCNS --------------------------------------
  
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
  
  # save to csv (used in Figure 3 S3 & 7 S3; Supp Table 3)
  if (write_csv) {
    csv_path <- file.path(PATH_input, "maleCNS_brain_NSC_meta_v0.9.csv")
    write.csv(maleCNS_brain_NSC, csv_path)
    cat("Saved CSV to:", normalizePath(csv_path), "\n")
  }
}

# Connectivity data for BANC & maleCNS -----------------------------------------

# download connectivity data for BANC 
el <- banc_edgelist(version = "v2", source = "gcs")

# edgelist contains autapses - filter out
el_filtered <- el %>%
  filter(pre != post)

# Save to csv (used in Figure 7 S3)
if (write_csv) {
  csv_path <- file.path(PATH_input, "BANC_edgelist_no_autapses_v888.csv")
  write.csv(el_filtered, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# choose dataset
neuprint_dataset <- "male-cns:v0.9"

# download connectivity data for maleCNS 
maleCNS_NSC_input <- mcns_connection_table(maleCNS_brain_NSC$bodyid, 
                                              partners='inputs', 
                                              moredetails = TRUE, 
                                              threshold = 1,
                                              roi = NULL,
                                              by.roi = FALSE)

maleCNS_NSC_input <- maleCNS_NSC_input %>%
  mutate(bodyid = as.character(bodyid)) %>%
  left_join(
    maleCNS_brain_NSC %>% select(bodyid, MckimType),
    by = "bodyid"
  )

# remove autapses if any
maleCNS_NSC_input <- maleCNS_NSC_input %>%
  filter(bodyid != partner)
# no autapses

# check naming
unique(maleCNS_NSC_input$superclass)

# adjust to match plotting code for colors, etc.
maleCNS_NSC_input <- maleCNS_NSC_input %>%
  mutate(superclass = case_when(
    superclass == "cb_intrinsic"       ~ "central_brain_intrinsic",
    superclass == "cb_sensory"         ~ "sensory",
    superclass == "ascending_neuron"   ~ "ascending",
    superclass == "descending_neuron"  ~ "descending",
    superclass == "cb_motor"           ~ "motor",
    superclass == "sensory_descending" ~ "sensory",
    TRUE                               ~ superclass
  ))

if (write_csv) {
  csv_path <- file.path(PATH_input, "maleCNS_NSC_input_v0.9.csv")
  write.csv(maleCNS_NSC_input, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Analyze input connectivity data for NSC in BANC
# ------------------------------------------------------------------------------
BANC_brain_NSC <- BANC_brain_NSC %>%
  select(pt_root_id, cell_class, cell_type) %>%
  rename(
    NSC_root_id = pt_root_id,
    NSC_cell_class = cell_class,
    NSC_cell_type  = cell_type)

BANC_brain_NSC$NSC_root_id_name <- paste(BANC_brain_NSC$NSC_root_id, 
                                         BANC_brain_NSC$NSC_cell_type, sep = "_")

# add input partners of NSC
BANC_NSC_input <- BANC_brain_NSC %>% 
  left_join(el_filtered, by = c("NSC_root_id" = "post"))

# add meta data for synaptic partners
BANC_NSC_input <- BANC_NSC_input %>%
  left_join(banc_meta, by = c("pre" = "pt_root_id"))

# filter synapses
BANC_NSC_input <- BANC_NSC_input %>%
  filter(count >= 5)

if (write_csv) {
  csv_path <- file.path(PATH_input, "BANC_NSC_input_v888.csv")
  write.csv(BANC_NSC_input, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#                                   BANC
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Panel A - BANC # of inputs
# ------------------------------------------------------------------------------
panel_name = '_S3A'
panel_type = "_BANC_NSC_input_summary_v888"
filename_banc = paste0(fig_folder_name, panel_name, panel_type)

na_rows <- BANC_NSC_input %>%
  filter(is.na(super_class))
# 34 cells with no super_class

# change the superclass to undefined for these
BANC_NSC_input <- BANC_NSC_input %>%
  mutate(super_class = ifelse(is.na(super_class), "undefined", super_class))

# table for plot, with:
# super_class, count per super_class, count per NSC, total synapses
connectivity_count <- BANC_NSC_input %>%
  group_by(super_class) %>%
  summarise(n_superclass = length(unique(pre)),
            n_NSC = length(unique(NSC_root_id)),
            n_synapses_total = sum(count, na.rm = TRUE))

if (write_csv) {
  csv_path <- file.path(PATH_output, fig_folder_name, paste0(filename_banc, ".csv"))
  write.csv(connectivity_count, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}


# Prep for plot 
# Categorize n_synapses_total into 4 bins
connectivity_count <- connectivity_count %>%
  mutate(synapse_bin = case_when(
    n_synapses_total >= 10000 ~ ">10000",
    n_synapses_total >= 1000  ~ "1000-10000",
    n_synapses_total >= 100   ~ "100-1000",
    TRUE                      ~ "<100"
  )) %>%
  # Set the order for the legend/mapping
  mutate(synapse_bin = factor(synapse_bin, 
                              levels = c(">10000", "1000-10000", "100-1000",
                                         "<100"))) %>%
  arrange(desc(n_superclass)) %>%
  mutate(y_pos = rev(row_number()))

# Define the manual thicknesses for each category
bin_widths <- c(">10000" = 5, "1000-10000" = 3, "100-1000" = 1.5, "<100" = 0.5)

p <- ggplot(connectivity_count) +
  # Main Arrows - mapped to the categorical 'synapse_bin'
  geom_segment(aes(x = 2.5, xend = 4.5, y = y_pos, yend = y_pos, 
                   linewidth = synapse_bin, color = super_class),
               lineend = "butt", linejoin = "mitre",
               arrow = arrow(length = unit(0.3, "cm"), type = "closed", angle = 20)) + 
  
  # Text Labels
  geom_text(aes(x = 2.4, y = y_pos, label = paste0(super_class, " (", n_superclass, ")")), 
            hjust = 1, size = 4.5) +
  geom_text(aes(x = 4.7, y = y_pos, label = n_NSC), 
            hjust = 0, size = 4.5) +
  
  # Headers
  annotate("text", x = 2.4, y = max(connectivity_count$y_pos) + 0.8, 
           label = "Super class\n(# neurons)", fontface = "bold", size = 4.5, hjust = 1) +
  annotate("text", x = 4.7, y = max(connectivity_count$y_pos) + 0.8, 
           label = "NSC\n(# neurons)", fontface = "bold", size = 4.5, hjust = 0) +
  
  # --- Categorial legend ---
  # Manually place four segments at y=0, space out across the x-axis
  annotate("text", x = 1.2, y = 0, label = "synapses", size = 4, fontface = "italic") +
  
  # Bin 1: >10000
  geom_segment(aes(x = 1.8, xend = 2.1, y = 0, yend = 0),
               linewidth = bin_widths[">10000"], color = "black") +
  annotate("text", x = 2.2, y = 0, label = ">10000", hjust = 0, size = 3.5) +
  
  # Bin 2: 1000-10000
  geom_segment(aes(x = 3.1, xend = 3.4, y = 0, yend = 0),
               linewidth = bin_widths["1000-10000"], color = "black") +
  annotate("text", x = 3.5, y = 0, label = "1000-10000", hjust = 0, size = 3.5) +
  
  # Bin 3: 100-1000
  geom_segment(aes(x = 4.7, xend = 5.0, y = 0, yend = 0),
               linewidth = bin_widths["100-1000"], color = "black") +
  annotate("text", x = 5.1, y = 0, label = "100-1000", hjust = 0, size = 3.5) +
  
  # Bin 4: <100
  geom_segment(aes(x = 6.0, xend = 6.3, y = 0, yend = 0),
               linewidth = bin_widths["<100"], color = "black") +
  annotate("text", x = 6.4, y = 0, label = "<100", hjust = 0, size = 3.5) +
  
  # Scales
  scale_color_manual(values = super_class_colors) +
  scale_linewidth_manual(values = bin_widths) +
  
  # Theme
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(30, 30, 30, 30)) +
  xlim(0, 7.5) + 
  ylim(-0.5, max(connectivity_count$y_pos) + 1.5)

print(p)

if (write_plots){
  # Save to pdf
  pdf_filename <- paste0(fig_folder_name, panel_name, ".pdf")
  pdf_path <- file.path(PATH_output, fig_folder_name, pdf_filename)
  ggsave(pdf_path, plot = p, width = 30, height = 14, units = "cm")
  cat("Saved combined PDF to:", normalizePath(pdf_path), "\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Panel B - BANC prop input synapses
# ------------------------------------------------------------------------------
panel_name = '_S3B'
panel_type = "_BANC_NSC_prop_input"
filename_banc_prop = paste0(fig_folder_name, panel_name, panel_type)

# Table for plot, with:
# super_class, NSC_cell_type, synapses sum, avg synapses, pre & post partners (unique)
NSC_input_sum_grouped <- BANC_NSC_input %>%
  group_by(super_class, NSC_cell_type) %>%
  summarise(n_synapses_sum = sum(count, na.rm = TRUE),
            avrg_synapses = mean(count, na.rm = TRUE),
            n_pre_partners = length(unique(pre)),
            n_post_partners = length(unique(NSC_root_id)),
            .groups = "drop")

# Table for plot, with:
# NSC_cell_type, # synapses total, # pre partners (unique)
NSC_input_total <- BANC_NSC_input %>%
  group_by(NSC_cell_type) %>%
  summarise(n_synapses_total = sum(count, na.rm = TRUE),
            n_pre_partners_total = length(unique(pre)))

# Join tables by NSC_cell_type
# Calc %age of input
NSC_input_sum_grouped <- left_join(NSC_input_sum_grouped, NSC_input_total, 
                                   by = "NSC_cell_type")
NSC_input_sum_grouped$perc_of_input <- NSC_input_sum_grouped$n_synapses_sum / 
                                       NSC_input_sum_grouped$n_synapses_total

# Adjust naming of unknown m and l neurons
NSC_input_sum_grouped$NSC_cell_type <- recode(NSC_input_sum_grouped$NSC_cell_type,
                                              "l_NSC" = "l_NSC_unknown",
                                              "m-NSC" = "m_NSC_unknown")

# no ITP "l_NSC_ITP"
NSC_input_sum_grouped$NSC_cell_type <- factor(NSC_input_sum_grouped$NSC_cell_type,
                                              levels = rev(c("m_NSC_unknown", "l_NSC_unknown", 
                                                             "m_NSC_DILP","m_NSC_DH44", 
                                                             "m_NSC_DMS", "l_NSC_CRZ",
                                                             "l_NSC_DH31", "SEZ_NSC_CAPA", 
                                                             "SEZ_NSC_Hugin")))

# Adjust classes for plotting - threshold for others & undefined
# Use fill_corr for plot colors 
NSC_input_sum_grouped$super_class_corr <- NSC_input_sum_grouped$super_class
NSC_input_sum_grouped$super_class_corr[NSC_input_sum_grouped$perc_of_input < 0.025] <- "others"
NSC_input_sum_grouped$super_class_corr[is.na(NSC_input_sum_grouped$super_class)] <- "undefined"
NSC_input_sum_grouped$fill_corr <- paste(NSC_input_sum_grouped$super_class_corr)

# simplify table for plot, with:
# NSC_cell_type, fill_corr (determines plot colors), % of input
NSC_input_sum_grouped <- NSC_input_sum_grouped %>%
  group_by(NSC_cell_type, fill_corr) %>%
  summarise(perc_of_input = sum(perc_of_input),
            .groups = "drop")

# set fill_corr as factor for plot
NSC_input_sum_grouped$fill_corr <- factor(NSC_input_sum_grouped$fill_corr,
                                          levels = rev(c("others",
                                                         "undefined",
                                                         "optic",
                                                         "descending",
                                                         "visual_centrifugal",
                                                         "central_brain_intrinsic",
                                                         "visual_projection",
                                                         "sensory",
                                                         "sensory_ascending",
                                                         "ascending",
                                                         "visceral_circulatory",
                                                         "motor")))

# Diagnostic print
# print(NSC_input_sum_grouped)

if (write_csv) {
  csv_path <- file.path(PATH_output, fig_folder_name,
                        paste0(filename_banc_prop, ".csv"))
  write.csv(NSC_input_sum_grouped, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# Add text to appear with superscript on x-axis
colored_labels <- setNames(
  paste0("<span style='color:", syn_col[levels(NSC_input_sum_grouped$NSC_cell_type)], "'>", 
         custom_labels_map[levels(NSC_input_sum_grouped$NSC_cell_type)], "</span>"),
  levels(NSC_input_sum_grouped$NSC_cell_type)
)

# no ITP data
NSC_input_sum_grouped$NSC_cell_type <- droplevels(NSC_input_sum_grouped$NSC_cell_type)

# Plot proportions
p <- ggplot(NSC_input_sum_grouped, aes(x = NSC_cell_type, y = perc_of_input,
                                       fill = fill_corr)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  facet_grid(scales = "free_x", space = "free_x") +
  scale_fill_manual(values = super_class_colors,
                    guide = guide_legend(nrow = 1)) +
  scale_x_discrete(limits = levels(NSC_input_sum_grouped$NSC_cell_type), 
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
        axis.text.x = element_markdown(angle = 270, hjust = 0, vjust = 0.5,
                                       colour = "black"),
        axis.text.y = element_text(angle = 0, colour = "black"),
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 8),
        legend.spacing.x = unit(0.1, "cm"))
print(p)

# Save
if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(fig_folder_name,
                                                             panel_name, ".pdf"))
  ggsave(pdf_path, plot = p, width = 12, height = 10, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#                                   maleCNS
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Panel C - maleCNS # of inputs
# ------------------------------------------------------------------------------
panel_name = '_S3C'
panel_type = "_maleCNS_NSC_input_summary_v0.9"
filename_mcns = paste0(fig_folder_name, panel_name, panel_type)

# threshold of 5 synapses
maleCNS_NSC_input <- maleCNS_NSC_input %>%
  filter(weight >= 5)

na_rows <- maleCNS_NSC_input %>%
  filter(is.na(superclass))
# 62 cells with no superclass # ** still the same

# change NA to undefined for superclass
maleCNS_NSC_input <- maleCNS_NSC_input %>%
  mutate(superclass = replace_na(superclass, "undefined"))

# table for plot, with:
# super_class, count per super_class, count per NSC, total synapses
connectivity_count = maleCNS_NSC_input%>%
  group_by(superclass)%>%
  summarise(n_superclass = length(unique(partner)),
            n_NSC = length(unique(bodyid)),
            n_synapses_total = sum(weight, na.rm = TRUE))

# Save
if (write_csv) {
  csv_path <- file.path(PATH_output, fig_folder_name, paste0(filename_mcns, ".csv"))
  write.csv(connectivity_count, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# Prep for plot 
# Categorize n_synapses_total into 4 bins
connectivity_count <- connectivity_count %>%
  mutate(synapse_bin = case_when(
    n_synapses_total >= 10000 ~ ">10000",
    n_synapses_total >= 1000  ~ "1000-10000",
    n_synapses_total >= 100   ~ "100-1000",
    TRUE                      ~ "<100"
  )) %>%
  # Set the order for the legend/mapping
  mutate(synapse_bin = factor(synapse_bin, 
                              levels = c(">10000", "1000-10000", "100-1000",
                                         "<100"))) %>%
  arrange(desc(n_superclass)) %>%
  mutate(y_pos = rev(row_number()))

# Define the manual thicknesses for each category
bin_widths <- c(">10000" = 5, "1000-10000" = 3, "100-1000" = 1.5, "<100" = 0.5)

p <- ggplot(connectivity_count) +
  geom_segment(aes(x = 2.5, xend = 4.5, y = y_pos, yend = y_pos, 
                   linewidth = synapse_bin, color = superclass),
               lineend = "butt", linejoin = "mitre",
               arrow = arrow(length = unit(0.3, "cm"), type = "closed", angle = 20)) + 
  
  geom_text(aes(x = 2.4, y = y_pos, label = paste0(superclass, " (", n_superclass, ")")), 
            hjust = 1, size = 4.5) +
  geom_text(aes(x = 4.7, y = y_pos, label = n_NSC), 
            hjust = 0, size = 4.5) +
  
  annotate("text", x = 2.4, y = max(connectivity_count$y_pos) + 0.8, 
           label = "Super class\n(# neurons)", fontface = "bold", size = 4.5, hjust = 1) +
  annotate("text", x = 4.7, y = max(connectivity_count$y_pos) + 0.8, 
           label = "NSC\n(# neurons)", fontface = "bold", size = 4.5, hjust = 0) +
  
  annotate("text", x = 1.2, y = 0, label = "synapses", size = 4, fontface = "italic") +
  
  geom_segment(aes(x = 1.8, xend = 2.1, y = 0, yend = 0),
               linewidth = bin_widths[">10000"], color = "black") +
  annotate("text", x = 2.2, y = 0, label = ">10000", hjust = 0, size = 3.5) +
  
  geom_segment(aes(x = 3.1, xend = 3.4, y = 0, yend = 0),
               linewidth = bin_widths["1000-10000"], color = "black") +
  annotate("text", x = 3.5, y = 0, label = "1000-10000", hjust = 0, size = 3.5) +
  
  geom_segment(aes(x = 4.7, xend = 5.0, y = 0, yend = 0),
               linewidth = bin_widths["100-1000"], color = "black") +
  annotate("text", x = 5.1, y = 0, label = "100-1000", hjust = 0, size = 3.5) +
  
  geom_segment(aes(x = 6.0, xend = 6.3, y = 0, yend = 0),
               linewidth = bin_widths["<100"], color = "black") +
  annotate("text", x = 6.4, y = 0, label = "<100", hjust = 0, size = 3.5) +
  
  scale_color_manual(values = super_class_colors) +
  scale_linewidth_manual(values = bin_widths) +
  
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(30, 30, 30, 30)) +
  xlim(0, 7.5) + 
  ylim(-0.5, max(connectivity_count$y_pos) + 1.5)

print(p)

if (write_plots){
  # Save to pdf
  pdf_filename <- paste0(fig_folder_name, panel_name, ".pdf")
  pdf_path <- file.path(PATH_output, fig_folder_name, pdf_filename)
  ggsave(pdf_path, plot = p, width = 30, height = 14, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}


# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Panel D - maleCNS prop input synapses
# ------------------------------------------------------------------------------
panel_name = '_S3D'
panel_type = "_maleCNS_NSC_prop_input"
filename_mcns_prop = paste0(fig_folder_name, panel_name, panel_type)

# Table for plot, with:
# super_class, NSC_cell_type, synapses sum, avg synapses, pre & post partners (unique)
NSC_input_sum_grouped <- maleCNS_NSC_input %>%
  group_by(superclass, MckimType) %>%
  summarise(n_synapses_sum = sum(weight, na.rm = TRUE),
            avrg_synapses = mean(weight, na.rm = TRUE),
            n_pre_partners = length(unique(partner)),
            n_post_partners = length(unique(bodyid)),
            .groups = "drop")

# Table for plot, with:
# NSC_cell_type, # synapses total, # pre partners (unique)
NSC_input_total <- maleCNS_NSC_input %>%
  group_by(MckimType) %>%
  summarise(n_synapses_total = sum(weight, na.rm = TRUE),
            n_pre_partners_total = length(unique(partner)))

# Join tables by NSC_cell_type
# Calc %age of input
NSC_input_sum_grouped <- left_join(NSC_input_sum_grouped, NSC_input_total, by = "MckimType")
NSC_input_sum_grouped$perc_of_input <- NSC_input_sum_grouped$n_synapses_sum / 
                                       NSC_input_sum_grouped$n_synapses_total

# NSC_cell_type as factor for plot order
NSC_input_sum_grouped$MckimType <- factor(NSC_input_sum_grouped$MckimType,
                                              levels = rev(c("m_NSC_unknown", "l_NSC_unknown",
                                                             "m_NSC_DILP","m_NSC_DH44", 
                                                             "m_NSC_DMS", "l_NSC_CRZ",
                                                             "l_NSC_ITP", "l_NSC_DH31", 
                                                             "SEZ_NSC_CAPA", "SEZ_NSC_Hugin")))

# Check before plotting
# unique(NSC_input_sum_grouped$superclass)

# Adjust classes for plotting - threshold for others & undefined
# Use fill_corr for plot colors 
NSC_input_sum_grouped$super_class_corr <- NSC_input_sum_grouped$superclass
NSC_input_sum_grouped$super_class_corr[NSC_input_sum_grouped$perc_of_input < 0.025] <- "others"
NSC_input_sum_grouped$super_class_corr[is.na(NSC_input_sum_grouped$superclass)] <- "undefined"
NSC_input_sum_grouped$fill_corr <- paste(NSC_input_sum_grouped$super_class_corr)

# simplify table for plot, with:
# NSC_cell_type, fill_corr (determines plot colors), % of input
NSC_input_sum_grouped <- NSC_input_sum_grouped %>%
  group_by(MckimType, fill_corr) %>%
  summarise(perc_of_input = sum(perc_of_input),
            .groups = "drop")

# set fill_corr as factor for plot
NSC_input_sum_grouped$fill_corr <- factor(NSC_input_sum_grouped$fill_corr,
                                          levels = rev(c("others",
                                                         "undefined",
                                                         "optic",
                                                         "descending",
                                                         "visual_centrifugal",
                                                         "central_brain_intrinsic",
                                                         "visual_projection",
                                                         "sensory",
                                                         "sensory_ascending",
                                                         "ascending",
                                                         "visceral_circulatory",
                                                         "endocrine",
                                                         "motor")))

# Diagnostic print
# print(NSC_input_sum_grouped)

if (write_csv) {
  csv_path <- file.path(PATH_output, fig_folder_name, paste0(filename_mcns_prop, ".csv"))
  write.csv(NSC_input_sum_grouped, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# Add text to appear with superscript on x-axis
colored_labels <- setNames(
  paste0("<span style='color:", syn_col[levels(NSC_input_sum_grouped$MckimType)], "'>", 
         custom_labels_map[levels(NSC_input_sum_grouped$MckimType)], "</span>"),
  levels(NSC_input_sum_grouped$MckimType)
)


# Plot proportions
p <- ggplot(NSC_input_sum_grouped, aes(x = MckimType, y = perc_of_input, fill = fill_corr)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  facet_grid(scales = "free_x", space = "free_x") +
  scale_fill_manual(values = super_class_colors,
                    guide = guide_legend(nrow = 1)) +
  scale_x_discrete(limits = levels(NSC_input_sum_grouped$MckimType), 
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
        axis.text.x = element_markdown(angle = 270, hjust = 0, vjust = 0.5,
                                       colour = "black"),
        axis.text.y = element_text(angle = 0, colour = "black"),
        axis.ticks.x = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 8),
        legend.spacing.x = unit(0.1, "cm"))
print(p)

# Save
if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(fig_folder_name,
                                                             panel_name, ".pdf"))
  ggsave(pdf_path, plot = p, width = 22, height = 14, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}
# ------------------------------------------------------------------------------