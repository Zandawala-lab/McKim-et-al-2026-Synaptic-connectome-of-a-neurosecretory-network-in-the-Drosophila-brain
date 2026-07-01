# -- Figure 3 - S4 -------------------------------------------------------------
# ------------------------------------------------------------------------------
# A & B - top input and proportion input to NSC by neurotransmitter
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Note: Assumes Figure_3.R was run to generate NSC_input file imported below
# Can be run if file downloaded from zenodo and accessible in input directory:
# NSC_input_classification_v783.csv
# ------------------------------------------------------------------------------

# load packages ----------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(ggtext)
library(patchwork)
library(cowplot)

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
  cat("Created figure folder at:", normalizePath(file.path(PATH_output,
                                                           fig_folder_name)), "\n")
}

# ------------------------------------------------------------------------------
# Options for saving plots & csv files
# ------------------------------------------------------------------------------
# When running script for the first time, set to TRUE:
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

neurotransmitter_colors <- c("ACH"       = "#B9529F",
                             "GLUT"      = "#51B848",
                             "GABA"      = "orange", 
                             "uncertain" = "black",
                             "other"     = "lightgrey")

# load data --------------------------------------------------------------------

# Filtered and processed data file generated from Figure_3_and_S5.R 
# and saved to input/
NSC_input <- read_delim(file.path(PATH_input, paste0("NSC_input_classification_v",
                                                     v, ".csv")),
                        delim = ",",
                        escape_double = FALSE,
                        col_types = cols(pre_pt_root_id = col_character(),
                                         post_pt_root_id = col_character()),
                        trim_ws = TRUE)

NSC <- read_delim(file.path(PATH_input, paste0("NSC_v", v, ".csv")),
                  delim = ",",
                  escape_double = FALSE,
                  col_types = cols(NSC_id = col_character()),
                  trim_ws = TRUE)

# ------------------------------------------------------------------------------
# Figure 3 S4A - Individual input to NSC types by neurotransmitter
# Individual plots, followed by combined plots
# ------------------------------------------------------------------------------
panel_name = "_S4A_"
filename_input_NSC = paste0(fig_folder_name, panel_name)

plot_list <- list()  # collect plots for combined figure

for (i in unique(NSC$NSC_name)) {
  # Get the NSC_ids for the current NSC_name
  NSC_ids <- NSC %>% filter(NSC_name == i) %>% pull(NSC_id)
  # Filter the data for the current NSC_ids
  NSC_subset <- NSC_input %>% filter(post_pt_root_id %in% NSC_ids)
  
  # Aggregate the data to get the sum of n_synapses for each pre_pt_root_id
  NSC_aggregated <- NSC_subset %>%
    group_by(pre_pt_root_id, nt_type_stringent_corr) %>%
    summarise(n_synapses = sum(n_synapses), .groups = 'drop')
  
  # Select top 20 pre_pt_root_id based on n_synapses
  top20_subset <- NSC_aggregated %>%
    group_by(pre_pt_root_id) %>%
    summarise(n_synapses = sum(n_synapses)) %>%
    top_n(20, n_synapses)
  
  NSC_subset_top20 <- NSC_aggregated %>%
    filter(pre_pt_root_id %in% top20_subset$pre_pt_root_id) %>%
    mutate(
      nt_type_stringent_corr = factor(nt_type_stringent_corr, 
                                      levels = names(neurotransmitter_colors))
    )
  
  # Skip if no data (otherwise warning)
  if (nrow(NSC_subset_top20) == 0) {
    message("Skipping ", i, " — no data after filtering")
    next
  }
  
  # Shorten ids for combined plotting
  # All start with the same 9 digits: 720575940 (last 9 are unique)
  NSC_subset_top20 <- NSC_subset_top20 %>%
    mutate(pre_pt_root_id_short = substr(as.character(pre_pt_root_id), 10, 18))
  
  # --- Individual plot (full axis labels) -------------------------------------
  b_individual <- ggplot(NSC_subset_top20, aes(x = fct_reorder(pre_pt_root_id,
                                                               desc(n_synapses)), 
                                               y = n_synapses, 
                                               fill = nt_type_stringent_corr)) +
    geom_col() +
    labs(x = "", y = "", title = i) +
    scale_fill_manual(values = neurotransmitter_colors) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, color="black"),
          axis.text.y = element_text(size = 8, color="black"),
          axis.line = element_line(color = "black"),
          plot.title = element_text(size = 11, hjust = 0.5, color="black"),
          legend.title = element_blank(),
          text = element_text(size = 12, color="black")) +
    scale_y_continuous(breaks = waiver(), n.breaks = 5, expand = c(0,0))
  
  if (write_plots){
    pdf_filename <- paste0(filename_input_NSC, "top20_individual_inputs_NT_",
                           i, ".pdf")
    pdf_path <- file.path(PATH_output, fig_folder_name, pdf_filename)
    ggsave(pdf_path, plot = b_individual, width = 6, height = 3)
    cat("Saved PDF to:", normalizePath(pdf_path), "\n")
  }
  
  # --- Combined-plot version (shortened ids) ----------------------------------
  b_combined <- ggplot(NSC_subset_top20, aes(x = fct_reorder(pre_pt_root_id_short, 
                                                             desc(n_synapses)), 
                                             y = n_synapses, 
                                             fill = nt_type_stringent_corr)) +
    geom_col() +
    labs(x = "", y = "", title = i) +
    scale_fill_manual(values = neurotransmitter_colors) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5,
                                     color = "black"),
          axis.text.y = element_text(size = 8, color="black"),
          axis.line = element_line(color = "black"),
          plot.title = element_text(size = 11, hjust = 0.5, color="black"),
          legend.position = "none",
          text = element_text(size = 12, color="black")) +
    scale_y_continuous(breaks = waiver(), n.breaks = 5, expand = c(0,0))
  
  plot_list[[i]] <- b_combined  # store for combined figure
}

# ------------------------------------------------------------------------------
# *Combined plot* - assembled from plot_list above
# ------------------------------------------------------------------------------

# Legend adjustments
all_nt_types <- unique(unlist(lapply(plot_list, 
                                     function(p) levels(droplevels(p$data$nt_type_stringent_corr)))))
has_na <- any(unlist(lapply(plot_list, function(p) any(is.na(p$data$nt_type_stringent_corr)))))

legend_categories <- all_nt_types
legend_colors <- neurotransmitter_colors[all_nt_types]

if (has_na) {
  legend_categories <- c(legend_categories, "NA")
  legend_colors <- c(legend_colors, "NA" = "grey")
}

legend_plot <- ggplot(data.frame(nt_type_stringent_corr = factor(legend_categories,
                                                                 levels = legend_categories)), 
                      aes(x = 1, y = 1, fill = nt_type_stringent_corr)) +
  geom_col() +
  scale_fill_manual(values = legend_colors, name = "Neurotransmitter", drop = TRUE) +
  theme(legend.position = "bottom")

shared_legend <- cowplot::get_legend(legend_plot)

# Set the order
desired_order <- rev(names(syn_col))
desired_order <- desired_order[desired_order %in% names(plot_list)]
plot_list <- plot_list[desired_order]

# Combined plot
combined_plot <- wrap_plots(plot_list, ncol = 2)
final_plot <- combined_plot / shared_legend + plot_layout(heights = c(20, 1))
print(final_plot)

if (write_plots){
  # Save to pdf
  pdf_filename <- paste0(filename_input_NSC, "top20_individual_inputs_combined_NT",
                         ".pdf")
  pdf_path <- file.path(PATH_output, fig_folder_name, pdf_filename)
  ggsave(pdf_path, plot = final_plot, width = 12, height = 15)
  cat("Saved combined PDF to:", normalizePath(pdf_path), "\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 3 S4B - Proportion of input to NSC types by neurotransmitter
# ------------------------------------------------------------------------------

# Organize and process data ----------------------------------------------------
# Summarize NSC input data - by NT
NSC_input_sum_NT <- NSC_input %>%
  group_by(nt_type_stringent_corr, name_post) %>%
  summarise(n_synapses_sum = sum(n_synapses, na.rm = TRUE),
            avrg_synapses = mean(n_synapses, na.rm = TRUE),
            n_pre_partners = length(unique(pre_pt_root_id)),
            n_post_partners = length(unique(post_pt_root_id)),
            .groups = "drop")

# Summarize NSC input data - total synapses & how many different pre neurons to NSC
NSC_input_total <- NSC_input %>%
  group_by(name_post) %>%
  summarise(n_synapses_total = sum(n_synapses, na.rm = TRUE),
            n_pre_partners_total = length(unique(pre_pt_root_id)))

# Combine groupings from above^ by NSC (post)
NSC_input_sum_NT <- left_join(NSC_input_sum_NT, NSC_input_total,
                              by = "name_post")
# Calculate % of input
NSC_input_sum_NT$perc_of_input <- NSC_input_sum_NT$n_synapses_sum / 
                                  NSC_input_sum_NT$n_synapses_total

NSC_input_sum_NT$name_post <- factor(NSC_input_sum_NT$name_post,
                                     levels = rev(c("m_NSC_unknown", "l_NSC_unknown",
                                                    "m_NSC_DILP", "m_NSC_DH44",
                                                    "m_NSC_DMS", "l_NSC_CRZ",
                                                    "l_NSC_ITP", "l_NSC_DH31",
                                                    "SEZ_NSC_CAPA", "SEZ_NSC_Hugin")))

NSC_input_sum_NT$fill <- paste(NSC_input_sum_NT$nt_type_stringent_corr, sep = "_")

NSC_input_sum_NT <- NSC_input_sum_NT %>%
  group_by(name_post, fill) %>%
  summarise(perc_of_input = sum(perc_of_input),
            .groups = "drop")

NSC_input_sum_NT$fill <- factor(NSC_input_sum_NT$fill,
                                levels = rev(c("uncertain",
                                               "other",
                                               "GABA",
                                               "GLUT",
                                               "ACH")))


# Diagnostic print
# print(NSC_input_sum_NT)
# ------------------------------------------------------------------------------

# Proportion plot --------------------------------------------------------------
panel_name = "_S4B"
filename_input_NT = paste0(fig_folder_name, panel_name)

if (write_csv) {
  csv_path <- file.path(PATH_output, fig_folder_name, paste0(filename_input_NT, 
                                                            "_v", v, ".csv"))
  write.csv(NSC_input_sum_NT, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# no ITP data
NSC_input_sum_NT$name_post <- droplevels(NSC_input_sum_NT$name_post)

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

colored_labels <- setNames(
  paste0("<span style='color:", syn_col[levels(NSC_input_sum_NT$name_post)], "'>", 
         custom_labels_map[levels(NSC_input_sum_NT$name_post)], "</span>"),
  levels(NSC_input_sum_NT$name_post)
)

q <- ggplot(NSC_input_sum_NT, aes(x = name_post, y = perc_of_input, fill = fill)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  facet_grid(scales = "free_x", space = "free_x") +
  scale_fill_manual(values = neurotransmitter_colors,
                    guide = guide_legend(nrow = 1)) +
  scale_x_discrete(limits = levels(NSC_input_sum_NT$name_post), 
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
print(q)

if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename_input_NT, ".pdf"))
  ggsave(pdf_path, plot = q, width = 12, height = 10, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}
# ------------------------------------------------------------------------------