# -- Figure 7 - S3 -------------------------------------------------------------
# ------------------------------------------------------------------------------
# Code for plotting synaptic outputs of NSCs in BANC v888 and maleCNS v0.9
# ------------------------------------------------------------------------------

# load packages ----------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(malecns)
library(neuprintr)
library(bancr)

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

super_class_colors <- c("others"                               = "#5f50a1", 
                        "undefined"                            = "#999999",
                        "optic"                                = "white",
                        "descending"                           = "#c3905f",  
                        "visual_centrifugal"                   = "#f46d43",
                        "central"                              = "#3384b8",
                        "central_brain_intrinsic"              = "#3384b8",
                        "visual_projection"                    = "#c8d684", 
                        "sensory"                              = "#b73545", 
                        "ascending"                            = "#a9d5a3",
                        "visceral_circulatory"                 = "#a8c8c0",
                        "endocrine"                            = "#a8c8c0",
                        "cb_endocrine"                         = "#a8c8c0",
                        "motor"                                = "black")

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

# Assumes Figure_1C_S1 (S1 section) and Figure_3_S3.R code have been run to
# download files
# if not, can go back and start there to obtain as needed (or download from zenodo)

# load NSC data - BANC (Fig 1 S1)
BANC_brain_NSC <- read_delim(file.path(PATH_input, "BANC_NSC_brain_meta_v888.csv"),
                             delim = ",",
                             escape_double = FALSE,
                             col_types = cols(pt_root_id = col_character()),
                             trim_ws = TRUE)

# load BANC meta (Fig 3 S3)
banc_meta <- read.csv(
  file.path(PATH_input, "BANC_meta_v888.csv"),
  colClasses = c(pt_root_id  = "character"),
  stringsAsFactors = FALSE,
  strip.white = TRUE)

# load connectivity data (BANC) (Fig 3 S3)
el_filtered <- read_delim(
  file.path(PATH_input, "BANC_edgelist_no_autapses_v888.csv"),
  delim = ",",
  escape_double = FALSE,
  col_types = cols(
    pre = col_character(), 
    post = col_character()),
  trim_ws = TRUE)

# load NSC data - maleCNS  (Fig 1 S1)
maleCNS_brain_NSC <- read_delim(file.path(PATH_input, "maleCNS_brain_NSC_meta_v0.9.csv"),
                                delim = ",",
                                escape_double = FALSE,
                                col_types = cols(bodyid = col_character()),
                                trim_ws = TRUE)
# ------------------------------------------------------------------------------
# Download connectivity data for maleCNS 
# ------------------------------------------------------------------------------
## Set your Neuprint token: https://natverse.org/neuprintr/#authentication
# choose dataset
neuprint_dataset <- "male-cns:v0.9"

# file also provided above instead 
maleCNS_NSC_output <- mcns_connection_table(maleCNS_brain_NSC$bodyid, 
                                           partners='outputs', 
                                           moredetails = TRUE, 
                                           threshold = 1,
                                           roi = NULL,
                                           by.roi = FALSE)

maleCNS_NSC_output <- maleCNS_NSC_output %>%
  mutate(bodyid = as.character(bodyid)) %>%
  left_join(
    maleCNS_brain_NSC %>% select(bodyid, MckimType),
    by = "bodyid"
  )

# remove autapses if any
maleCNS_NSC_output <- maleCNS_NSC_output %>%
  filter(bodyid != partner)
# no autapses

# check naming
unique(maleCNS_NSC_output$superclass)

# adjust to match plotting code for colors, etc.
maleCNS_NSC_output <- maleCNS_NSC_output %>%
  mutate(superclass = case_when(
    superclass == "cb_intrinsic"       ~ "central_brain_intrinsic",
    superclass == "cb_sensory"         ~ "sensory",
    superclass == "ascending_neuron"   ~ "ascending",
    superclass == "descending_neuron"  ~ "descending",
    superclass == "cb_motor"           ~ "motor",
    superclass == "sensory_descending" ~ "sensory",
    TRUE                               ~ superclass
  ))


# save to csv
if (write_csv) {
  csv_path <- file.path(PATH_input, "maleCNS_NSC_output_v0.9.csv")
  write.csv(maleCNS_NSC_output, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Analyze output connectivity data for NSC in BANC
# ------------------------------------------------------------------------------
BANC_brain_NSC <- BANC_brain_NSC %>%
  select(pt_root_id, cell_class, cell_type) %>%
  rename(
    NSC_root_id = pt_root_id,
    NSC_cell_class = cell_class,
    NSC_cell_type  = cell_type)

BANC_brain_NSC$NSC_root_id_name <- paste(BANC_brain_NSC$NSC_root_id, 
                                         BANC_brain_NSC$NSC_cell_type, sep = "_")

# add output partners of NSC
BANC_NSC_output <- BANC_brain_NSC %>% 
  left_join(el_filtered, by = c("NSC_root_id" = "pre"))

# add meta data for synaptic partners
BANC_NSC_output <- BANC_NSC_output %>%
  left_join(banc_meta, by = c("post" = "pt_root_id"))

# filter synapses
BANC_NSC_output <- BANC_NSC_output %>%
  filter(count >= 5)

if (write_csv) {
  csv_path <- file.path(PATH_input, "BANC_NSC_output_v888.csv")
  write.csv(BANC_NSC_output, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#                                   BANC
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Panel A - BANC outputs
# ------------------------------------------------------------------------------
panel_name = '_S3A'
panel_type = "_BANC_NSC_output_summary_v888"
filename_banc = paste0(fig_folder_name, panel_name, panel_type)

na_rows <- BANC_NSC_output %>%
  filter(is.na(super_class))
# all entries have a super class

# table for plot, with:
# super_class, count per super_class, count per NSC, total synapses
connectivity_count = BANC_NSC_output %>%
  group_by(super_class) %>%
  summarise(n_superclass = length(unique(post)),
            n_NSC = length(unique(NSC_root_id)),
            n_synapses_total = sum(count, na.rm = TRUE))

# save as csv
if (write_csv){
  csv_path <- file.path(PATH_output, fig_folder_name, 
                        paste0(fig_folder_name, "_S3_BANC_NSC_output_summary_v888.csv"))
  write.csv(connectivity_count, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# Sort by n_superclass descending - needed for plot
connectivity_count <- connectivity_count %>%
  arrange(desc(n_superclass)) %>%
  mutate(y_pos = row_number())

# Reverse y_pos so the highest n_superclass is at the top - needed for plot
connectivity_count$y_pos <- rev(connectivity_count$y_pos)

# Define scale limit slightly above max
max_limit <- 25

# Plot
p <- ggplot(connectivity_count) +
  geom_segment(aes(x = 2.5, xend = 4.5, y = y_pos, yend = y_pos, 
                   linewidth = n_synapses_total, color = super_class),
               lineend = "butt", 
               linejoin = "mitre",
               arrow = arrow(length = unit(0.35, "cm"), 
                             type = "closed", 
                             angle = 20)) + 
  
  geom_text(aes(x = 2.4, y = y_pos, label = n_NSC), 
            hjust = 1, size = 4.5) +
  
  geom_text(aes(x = 4.7, y = y_pos, 
                label = paste0(super_class, " (", n_superclass, ")")), 
            hjust = 0, size = 4.5) +
  
  annotate("text", x = 2.4, y = max(connectivity_count$y_pos) + 0.8, 
           label = "NSC\n(# neurons)", fontface = "bold", size = 4.5, hjust = 1) +
  annotate("text", x = 4.7, y = max(connectivity_count$y_pos) + 0.8, 
           label = "Super class\n(# neurons)", fontface = "bold", size = 4.5,
           hjust = 0) +
  
  annotate("text", x = 1.6, y = 0, label = "synapses", size = 4,
           fontface = "italic") +
  
  # High Benchmark (20)
  geom_segment(aes(x = 2.4, xend = 2.7, y = 0, yend = 0, linewidth = 20), 
               color = "black", lineend = "butt", data = data.frame(x=1)) +
  annotate("text", x = 2.8, y = 0, label = "20", hjust = 0, size = 3.5) +
  
  # Low Benchmark (5)
  geom_segment(aes(x = 3.6, xend = 3.9, y = 0, yend = 0, linewidth = 5), 
               color = "black", lineend = "butt", data = data.frame(x=1)) +
  annotate("text", x = 4.0, y = 0, label = "5", hjust = 0, size = 3.5) +
  
  scale_color_manual(values = super_class_colors) +
  scale_linewidth_continuous(range = c(0.8, 7), limits = c(0, max_limit)) +
  
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(legend.position = "none", 
        plot.margin = margin(30, 40, 30, 40)) +
  xlim(0, 8) +
  ylim(-0.5, max(connectivity_count$y_pos) + 1.5)

print(p)

if (write_plots){
  # Save to pdf
  pdf_filename <- paste0(fig_folder_name, panel_name, ".pdf")
  pdf_path <- file.path(PATH_output, fig_folder_name, pdf_filename)
  ggsave(pdf_path, plot = p, width = 32, height = 12, units = "cm")
  cat("Saved combined PDF to:", normalizePath(pdf_path), "\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Panel B - BANC prop output synapses
# ------------------------------------------------------------------------------
panel_name = '_S3B'
panel_type = "_BANC_NSC_prop_input"
filename_banc_prop = paste0(fig_folder_name, panel_name, panel_type)

# Table for plot, with:
# super_class, NSC_cell_type, synapses sum, avg synapses, pre & post partners (unique)
NSC_output_sum_grouped <- BANC_NSC_output %>%
  group_by(super_class, NSC_cell_type) %>%
  summarise(n_synapses_sum = sum(count, na.rm = TRUE),
            avrg_synapses = mean(count, na.rm = TRUE),
            n_pre_partners = length(unique(NSC_root_id)),
            n_post_partners = length(unique(post)),
            .groups = "drop")

# Table for plot, with:
# NSC_cell_type, # synapses total, # pre partners (unique)
NSC_output_total <- BANC_NSC_output %>%
  group_by(NSC_cell_type) %>%
  summarise(n_synapses_total = sum(count, na.rm = TRUE),
            n_pre_partners_total = length(unique(NSC_root_id)))

# Join tables by NSC_cell_type
# Calc %age of input
NSC_output_sum_grouped <- left_join(NSC_output_sum_grouped, NSC_output_total, 
                                    by = "NSC_cell_type")
NSC_output_sum_grouped$perc_of_output <- NSC_output_sum_grouped$n_synapses_sum / 
                                         NSC_output_sum_grouped$n_synapses_total

# NSC_cell_type as factor for plot order
NSC_output_sum_grouped$NSC_cell_type <- factor(NSC_output_sum_grouped$NSC_cell_type,
                                              levels = rev(c("m_NSC_unknown", "l_NSC_unknown", 
                                                             "m_NSC_DILP","m_NSC_DH44", 
                                                             "m_NSC_DMS", "l_NSC_CRZ",
                                                             "l_NSC_ITP", "l_NSC_DH31", 
                                                             "SEZ_NSC_CAPA", "SEZ_NSC_Hugin")))

# check category
# unique(NSC_output_sum_grouped$super_class)

# Adjust classes for plotting - threshold for others & undefined
# Use fill_corr for plot colors 
NSC_output_sum_grouped$super_class_corr <- NSC_output_sum_grouped$super_class
NSC_output_sum_grouped$super_class_corr[NSC_output_sum_grouped$perc_of_output < 0.025] <- "others"
NSC_output_sum_grouped$super_class_corr[is.na(NSC_output_sum_grouped$super_class)] <- "undefined"
NSC_output_sum_grouped$fill_corr <- paste(NSC_output_sum_grouped$super_class_corr)

# simplify table for plot, with:
# NSC_cell_type, fill_corr (determines plot colors), % of output
NSC_output_sum_grouped <- NSC_output_sum_grouped %>%
  group_by(NSC_cell_type, fill_corr) %>%
  summarise(perc_of_output = sum(perc_of_output),
            .groups = "drop")

# set fill_corr as factor for plot
NSC_output_sum_grouped$fill_corr <- factor(NSC_output_sum_grouped$fill_corr,
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
# print(NSC_output_sum_grouped)

if (write_csv) {
  csv_path <- file.path(PATH_output, fig_folder_name, paste0(filename_banc_prop,
                                                             ".csv"))
  write.csv(NSC_output_sum_grouped, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}


# Add text to appear with superscript on x-axis
colored_labels <- setNames(
  paste0("<span style='color:", syn_col[levels(NSC_output_sum_grouped$NSC_cell_type)], "'>", 
         custom_labels_map[levels(NSC_output_sum_grouped$NSC_cell_type)], "</span>"),
  levels(NSC_output_sum_grouped$NSC_cell_type)
)

# drop unused levels (in this case, only have CRZ that appears)
NSC_output_sum_grouped$NSC_cell_type <- droplevels(NSC_output_sum_grouped$NSC_cell_type)

# Plot proportions
p <- ggplot(NSC_output_sum_grouped, aes(x = NSC_cell_type, y = perc_of_output,
                                        fill = fill_corr)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  facet_grid(scales = "free_x", space = "free_x") +
  scale_fill_manual(values = super_class_colors,
                    guide = guide_legend(nrow = 1)) +
  scale_x_discrete(limits = levels(NSC_output_sum_grouped$NSC_cell_type), 
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

# Save
if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(fig_folder_name,
                                                             panel_name, ".pdf"))
  ggsave(pdf_path, plot = p, width = 10, height = 13, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
#                                   maleCNS
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Panel C - maleCNS outputs
# ------------------------------------------------------------------------------
panel_name = '_S3C'
panel_type = "_maleCNS_NSC_output_summary_v0.9"
filename_mcns = paste0(fig_folder_name, panel_name, panel_type)

# threshold of 5 synapses
maleCNS_NSC_output <- maleCNS_NSC_output %>%
  filter(weight >= 5)

na_rows <- maleCNS_NSC_output %>%
  filter(is.na(superclass))
# 0 cells with no superclass

# table for plot, with:
# super_class, count per super_class, count per NSC, total synapses
connectivity_count = maleCNS_NSC_output%>%
  group_by(superclass)%>%
  summarise(n_superclass = length(unique(partner)),
            n_NSC = length(unique(bodyid)),
            n_synapses_total = sum(weight, na.rm = TRUE))

# Save
if (write_csv) {
  csv_path <- file.path(PATH_output, fig_folder_name, paste0(filename_mcns,
                                                             ".csv"))
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

# Define scale limit slightly above max
max_limit <- 317

# Plot
p <- ggplot(connectivity_count) +
  # Main arrows
  geom_segment(aes(x = 2.5, xend = 4.5, y = y_pos, yend = y_pos, 
                   linewidth = n_synapses_total, color = superclass),
               lineend = "butt", 
               linejoin = "mitre",
               arrow = arrow(length = unit(0.35, "cm"), 
                             type = "closed", 
                             angle = 20)) + 
  
  # Left: NSC (# neurons)
  # Positioned at 2.4 (left of arrow start)
  geom_text(aes(x = 2.4, y = y_pos, label = n_NSC), 
            hjust = 1, size = 4.5) +
  
  # Right: Super class (# neurons)
  # Positioned at 4.7 (right of arrow head)
  geom_text(aes(x = 4.7, y = y_pos, 
                label = paste0(superclass, " (", n_superclass, ")")), 
            hjust = 0, size = 4.5) +
  
  # Swapped Headers
  annotate("text", x = 2.4, y = max(connectivity_count$y_pos) + 0.8, 
           label = "NSC\n(# neurons)", fontface = "bold", size = 4.5, hjust = 1) +
  annotate("text", x = 4.7, y = max(connectivity_count$y_pos) + 0.8, 
           label = "Super class\n(# neurons)", fontface = "bold", size = 4.5,
           hjust = 0) +
  
  # --- Continuous Legend ---
  annotate("text", x = 1.6, y = 0, label = "synapses", size = 4, fontface = "italic") +
  
  # High benchmark (300)
  geom_segment(aes(x = 2.4, xend = 2.7, y = 0, yend = 0, linewidth = 300), 
               color = "black", lineend = "butt", data = data.frame(x=1)) +
  annotate("text", x = 2.8, y = 0, label = "300", hjust = 0, size = 3.5) +
  
  # Low benchmark (50)
  geom_segment(aes(x = 3.6, xend = 3.9, y = 0, yend = 0, linewidth = 50), 
               color = "black", lineend = "butt", data = data.frame(x=1)) +
  annotate("text", x = 4.0, y = 0, label = "50", hjust = 0, size = 3.5) +
  
  # Scales
  scale_color_manual(values = super_class_colors) +
  scale_linewidth_continuous(range = c(0.8, 7), limits = c(0, max_limit)) +
  
  # Formatting
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(legend.position = "none", 
        plot.margin = margin(30, 40, 30, 40)) +
  xlim(0, 8) +
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
# Panel D - maleCNS prop output synapses
# ------------------------------------------------------------------------------
panel_name = '_S3D'
panel_type = "_maleCNS_NSC_prop_input"
filename_mcns_prop = paste0(fig_folder_name, panel_name, panel_type)


# Table for plot, with:
# super_class, NSC_cell_type, synapses sum, avg synapses, pre & post partners (unique)
NSC_output_sum_grouped <- maleCNS_NSC_output %>%
  group_by(superclass, MckimType) %>%
  summarise(n_synapses_sum = sum(weight, na.rm = TRUE),
            avrg_synapses = mean(weight, na.rm = TRUE),
            n_pre_partners = length(unique(bodyid)),
            n_post_partners = length(unique(partner)),
            .groups = "drop")

# Table for plot, with:
# NSC_cell_type, # synapses total, # pre partners (unique)
NSC_output_total <- maleCNS_NSC_output %>%
  group_by(MckimType) %>%
  summarise(n_synapses_total = sum(weight, na.rm = TRUE),
            n_pre_partners_total = length(unique(bodyid)))

# Join tables by NSC_cell_type
# Calc %age of input
NSC_output_sum_grouped <- left_join(NSC_output_sum_grouped, NSC_output_total, 
                                    by = "MckimType")
NSC_output_sum_grouped$perc_of_output <- NSC_output_sum_grouped$n_synapses_sum / 
                                         NSC_output_sum_grouped$n_synapses_total

# NSC_cell_type as factor for plot order
NSC_output_sum_grouped$MckimType <- factor(NSC_output_sum_grouped$MckimType,
                                          levels = rev(c("m_NSC_unknown", "l_NSC_unknown",
                                                         "m_NSC_DILP","m_NSC_DH44", 
                                                         "m_NSC_DMS", "l_NSC_CRZ",
                                                         "l_NSC_ITP", "l_NSC_DH31", 
                                                         "SEZ_NSC_CAPA", "SEZ_NSC_Hugin")))

# Check before plotting
# unique(NSC_output_sum_grouped$superclass)

# Adjust classes for plotting - threshold for others & undefined
# Use fill_corr for plot colors 
NSC_output_sum_grouped$super_class_corr <- NSC_output_sum_grouped$superclass
NSC_output_sum_grouped$super_class_corr[NSC_output_sum_grouped$perc_of_output < 0.025] <- "others"
NSC_output_sum_grouped$super_class_corr[is.na(NSC_output_sum_grouped$superclass)] <- "undefined"
NSC_output_sum_grouped$fill_corr <- paste(NSC_output_sum_grouped$super_class_corr)

# simplify table for plot, with:
# NSC_cell_type, fill_corr (determines plot colors), % of output
NSC_output_sum_grouped <- NSC_output_sum_grouped %>%
  group_by(MckimType, fill_corr) %>%
  summarise(perc_of_output = sum(perc_of_output),
            .groups = "drop")

# set fill_corr as factor for plot
NSC_output_sum_grouped$fill_corr <- factor(NSC_output_sum_grouped$fill_corr,
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
                                                         "cb_endocrine",
                                                         "motor")))

# Diagnostic print
#print(NSC_output_sum_grouped)


if (write_csv) {
  csv_path <- file.path(PATH_output, fig_folder_name, paste0(filename_mcns_prop,
                                                             ".csv"))
  write.csv(NSC_output_sum_grouped, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}


# Add text to appear with superscript on x-axis
colored_labels <- setNames(
  paste0("<span style='color:", syn_col[levels(NSC_output_sum_grouped$MckimType)], "'>", 
         custom_labels_map[levels(NSC_output_sum_grouped$MckimType)], "</span>"),
  levels(NSC_output_sum_grouped$MckimType)
)

# drop unused levels
NSC_output_sum_grouped$MckimType <- droplevels(NSC_output_sum_grouped$MckimType)

# Plot proportions
p <- ggplot(NSC_output_sum_grouped, aes(x = MckimType, y = perc_of_output,
                                        fill = fill_corr)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  facet_grid(scales = "free_x", space = "free_x") +
  scale_fill_manual(values = super_class_colors,
                    guide = guide_legend(nrow = 1)) +
  scale_x_discrete(limits = levels(NSC_output_sum_grouped$MckimType), 
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

# Save
if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(fig_folder_name,
                                                             panel_name, ".pdf"))
  ggsave(pdf_path, plot = p, width = 22, height = 14, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}

# ------------------------------------------------------------------------------