# --- Figure 5 -----------------------------------------------------------------
# ------------------------------------------------------------------------------
# Code for full pathway from ORN->AL*->INTs->NSCs
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
fig_folder_name = "Figure_5"
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

neurotransmitter_colors <- c("ACH"       = "#B9529F",
                             "GLUT"      = "#51B848",
                             "GABA"      = "orange", 
                             "uncertain" = "black",
                             "other"     = "lightgrey")

behav_signif <- c("Aversive"  = "#ED2024",
                  "Egg-laying" = "#3953A4", 
                  "Food"       = "#6ABD45", 
                  "Pheromonal" = "#8150A0",
                  "Unknown"    = "black")

# load data:--------------------------------------------------------------------
# Filtered and processed data file generated from Figure_3.R 
# saved to input/ & data is: INTs -> NSCs
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

# This is the filtered datafile: AL* -> INTs:
al_interneurons_input <- read_csv(file.path(PATH_input, 
                                            paste0("al_interneurons_input_filtered_v", 
                                                   v, ".csv")),
                      col_types = cols(
                        pre_pt_root_id = col_character(),
                        post_pt_root_id = col_character()
                      ))

# Query caveclient w/python to get ORN input to AL*
if (!paste0("input_to_interneurons_filtered_v",v,".csv")  %in%  input_files) {
  reticulate::source_python("./scripts/Figure_5/Figure_5_preparation.py", 
                            envir = NULL)
}
al_interneurons_input_curated <- read_csv(file.path(PATH_input,
                                                    paste0("tmp/input_to_interneurons_filtered_v",
                                                           v, ".csv")),
                               col_types = cols(
                                 pre_pt_supervoxel_id = col_character(),
                                 pre_pt_root_id = col_character(),
                                 post_pt_supervoxel_id = col_character(),
                                 post_pt_root_id = col_character()
                               ))

# Check for classification file
if (!paste0("classification_v", v, ".csv") %in% list.files(PATH_input)) {
  stop("please go to zenodo and download the classification file provided (see 
       note at top of script) and save it in './input'.")
}

classification <- read_delim(file.path(PATH_input, paste0("classification_v", v, ".csv")),
                             delim = ",",
                             escape_double = FALSE,
                             col_types = cols(root_id = col_character(), 
                                              flow = col_character()),
                             trim_ws = TRUE)

# Update classification
classification[classification$hemibrain_type %in% c("HBeyelet"),]$super_class <- "sensory"


# Pull in behavioral classification file (provided with files on zenodo)
load(file.path(PATH_input, file='glomeruli and behaviors.rda'))
colnames(glom_and_behav) <- c("PN_types", "Behav_Signif", "Lit") 


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

# Reclassify neurotransmitters based on stringent cutoff
neurotransmitters$nt_type_stringent <- neurotransmitters$nt_type
neurotransmitters <- neurotransmitters %>%
  mutate(nt_type_stringent = ifelse(nt_type_score <= 0.62, "uncertain", 
                                    nt_type_stringent))

# Collapse neurotransmitters into 5 categories (GABA, GLUT, ACH, uncertain, other)
neurotransmitters$nt_type_stringent_corr <- neurotransmitters$nt_type_stringent
neurotransmitters <- neurotransmitters %>%
  mutate(nt_type_stringent_corr = case_when(
    nt_type_stringent_corr %in% c("GABA", "GLUT", "ACH", "uncertain") ~ nt_type_stringent_corr,
    TRUE ~ "other" ))
  
# Data processing and org ------------------------------------------------------

# Summarize synapses data - ORN -> AL*
input_to_al_neurons_sum <- al_interneurons_input_curated %>%
  group_by(pre_pt_root_id, post_pt_root_id) %>%
  summarise(n_synapses = n(),
            .groups = "drop")

if (write_csv){
  csv_path <- file.path(PATH_input, paste0("tmp/input_to_al_neurons_sum_v",
                                           v, ".csv"))
  write.csv(input_to_al_neurons_sum, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# Join data for further analysis
classification_join <- classification
colnames(classification_join) <- c("pre_pt_root_id", colnames(classification)[-1])
input_to_al_neurons <- left_join(input_to_al_neurons_sum[input_to_al_neurons_sum$n_synapses >= 5, ], 
                                 classification_join, by = "pre_pt_root_id")

if (write_csv){
  csv_path <- file.path(PATH_output, fig_folder_name,
                        paste0("input_to_al_neurons_allinput_threshold5_v",
                               v, ".csv"))
  write.csv(input_to_al_neurons, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# Subset classification and adjust col names to include 'al' for neuron info (post_pt_root_id)
classification_al <- select(classification, root_id, super_class, class,
                            sub_class, cell_type, hemibrain_type)
colnames(classification_al) <- c("post_pt_root_id", "al_super_class","al_class", 
                                 "al_subclass", "al_cell_type", "al_hemibrain_type")
input_to_al_neurons <- left_join(input_to_al_neurons, classification_al, 
                                 by = "post_pt_root_id")

# Check how many rows with NA for pre
rows_na_input <- input_to_al_neurons[is.na(input_to_al_neurons$super_class),] #24
# No na rows for post
rows_na_al <- input_to_al_neurons[is.na(input_to_al_neurons$al_super_class),]

# Remove inputs that are orphan fragments (marked by NA under super_class)
input_to_al_neurons_filtered <- input_to_al_neurons %>% drop_na(super_class)

# Assign values to hemibrain_type from cell_type or class if they are NA
input_to_al_neurons_filtered$hemibrain_type <- ifelse(is.na(input_to_al_neurons_filtered$hemibrain_type),
                                                      ifelse(is.na(input_to_al_neurons_filtered$cell_type),
                                                             input_to_al_neurons_filtered$class,
                                                             input_to_al_neurons_filtered$cell_type),
                                                      input_to_al_neurons_filtered$hemibrain_type)
# Check how many - 7
rows_na_hemi <- input_to_al_neurons_filtered[is.na(input_to_al_neurons_filtered$hemibrain_type),]

# Filter out rows where hemibrain_type is still NA
input_to_al_neurons_filtered <- input_to_al_neurons_filtered[!is.na(input_to_al_neurons_filtered$hemibrain_type), ]  

# Subset based on super_class sensory
sensory_input_to_al_neurons_filtered <- input_to_al_neurons_filtered %>% 
                                       filter(grepl("sensory", super_class, 
                                                    ignore.case = TRUE)) 
# 435 pre / 40 post (unique)


# Filter based on olf first - 505 rows
olf_input_to_al_neurons <- input_to_al_neurons_filtered[grep("olfactory", 
                                                             input_to_al_neurons_filtered$class), ] 
# Add the behav and lit info
olf_names <- as.data.frame(str_split_fixed(olf_input_to_al_neurons$hemibrain_type, '_',2))
colnames(olf_names) <- c("ORN", "PN_types")
# Merge PN types col to the orig df
olf_input_to_al_neurons$PN_types <- olf_names$PN_types
# Join olf df with collated info from glom_and_behav
olf_input_to_al_neurons <-left_join(olf_input_to_al_neurons, glom_and_behav, 
                                    by="PN_types")

# ------------------------------------------------------------------------------
# Combine data to get full pathway:  ORN -> AL* -> INTs -> NSCs
# ------------------------------------------------------------------------------
unique_olf <- unique(olf_input_to_al_neurons$pre_pt_root_id) #321
unique_al <- unique(olf_input_to_al_neurons$post_pt_root_id) #14

# now take the ids from unique_al and find them in: AL* -> INTs
matching_rows <- al_interneurons_input %>%
  filter(pre_pt_root_id %in% unique_al)

# 15 rows
unique_al2 <- unique(matching_rows$pre_pt_root_id) #14
unique_int <- unique(matching_rows$post_pt_root_id) #7

# now take the ids from unique_int and find them in: INTs -> NSC
# 10 rows
matching_ints<- NSC_input %>%
  filter(pre_pt_root_id %in% unique_int)

unique_int2 <- unique(matching_ints$pre_pt_root_id) #7
unique_NSC <- unique(matching_ints$post_pt_root_id) #8

# organize int_to_NSC col names (INTs -> NSC) ----------------------------------
# remove empty
int_to_NSC <- matching_ints %>%
  select(-`...1`, -name_pre, -hemisphere_pre, -nerve)

stay_as_is <- c("post_pt_root_id", "name_post", "hemisphere_post", 
                "post_pt_root_id_name")

pre_cols <- setdiff(colnames(int_to_NSC), stay_as_is)

int_to_NSC <- int_to_NSC %>%
  rename_with(~ paste0("int_", .x), .cols = all_of(pre_cols))

# Change post to nsc
colnames(int_to_NSC) <- gsub("post", "nsc", colnames(int_to_NSC))

# final change- pre (INT) becomes post so that it can be matched in next step 
# to dataframe with AL* (pre) -> INT (post)
int_to_NSC <- int_to_NSC %>%
  rename(post_pt_root_id = int_pre_pt_root_id,
         nsc_id = nsc_pt_root_id)

connections <- full_join(matching_rows, 
                         int_to_NSC %>% select(-int_super_class, -int_class,
                                               -int_cell_type, -int_hemibrain_type), 
                         by = "post_pt_root_id", 
                         relationship = "many-to-many")

# organize AL* -> INTs ---------------------------------------------------------
connect <- connections  # duplicate dataset to preserve original

# remove empty
connect <- connections %>%
  select(-`...1`, -nerve)

exclude <- c("pre_pt_root_id", "post_pt_root_id")

al_cols <- colnames(connect)[!grepl("int_|nsc", colnames(connect)) & 
                               !colnames(connect) %in% exclude]

connect <- connect %>%
  rename_with(~ paste0("al_", .x), .cols = all_of(al_cols))

# final change - post (INT) becomes int_id and pre (AL) becomes post so can
# be matched in next step to dataframe with ORNs (pre) -> AL* (post)
connect <- connect %>%
  rename(int_id = post_pt_root_id,
         post_pt_root_id = pre_pt_root_id)

# let's sort to pretty it up
al_cols <- colnames(connect)[grepl("^al_", colnames(connect))]
int_cols <- colnames(connect)[grepl("^int_", colnames(connect))]
nsc_cols <- colnames(connect)[grepl("nsc", colnames(connect), ignore.case = TRUE)]

desired_order <- c(setdiff(colnames(connect), c(al_cols, int_cols, nsc_cols)), 
                   al_cols, int_cols, nsc_cols)

# Rearrange the dataframe columns
connect_org <- connect[, desired_order]

# now join again
# now need a file that has unique 
# olf_input_to_al_neurons #321 ORNs -> 14 ALs

all_connects <- full_join(olf_input_to_al_neurons, 
                          connect_org %>% select(-al_super_class, -al_class,
                                                 -al_cell_type, -al_hemibrain_type), 
                          by = "post_pt_root_id", relationship = "many-to-many")

# organize ORN -> AL* ----------------------------------------------------------
all_data <- all_connects # duplicate dataset to preserve original

# remove empty
all_data <- all_connects %>%
  select(-sub_class, -hemilineage)

# everything without al/int/nsc label gets orn_ prefix,
# except pre_pt_root_id and post_pt_root_id
exclude <- c("pre_pt_root_id", "post_pt_root_id")

orn_cols <- colnames(all_data)[!grepl("al_|int_|nsc", colnames(all_data)) & 
                                 !colnames(all_data) %in% exclude]

all_data <- all_data %>%
  rename_with(~ paste0("orn_", .x), .cols = all_of(orn_cols))

# final change - post (AL) becomes al_id and pre (ORN) becomes orn_id
all_data <- all_data %>%
  rename(al_id = post_pt_root_id,
         orn_id = pre_pt_root_id)

# organize final order - left to right: ORN -> AL* -> INT -> NSC) --------------

# let's sort to pretty it up
al_cols <- colnames(all_data)[grepl("^al_", colnames(all_data))]
# put id first
ordered_al <- c(al_cols[which(al_cols == "al_id")], setdiff(al_cols, "al_id"))
int_cols <- colnames(all_data)[grepl("^int_", colnames(all_data))]
# put id first
ordered_int <- c(int_cols[which(int_cols == "int_id")], setdiff(int_cols, "int_id"))
orn_cols <- colnames(all_data)[grepl("^orn_", colnames(all_data))]
# put id first
ordered_orn <- c(orn_cols[which(orn_cols == "orn_id")], setdiff(orn_cols, "orn_id"))
nsc_cols <- colnames(all_data)[grepl("nsc", colnames(all_data), ignore.case = TRUE)]
# put id first
ordered_nsc <- c(nsc_cols[which(nsc_cols == "nsc_id")], setdiff(nsc_cols, "nsc_id"))

desired_order <- c(ordered_orn, ordered_al, ordered_int, ordered_nsc)

# Rearrange the dataframe columns
all_data <- all_data[, desired_order]

# Save as csv
if (write_csv){
  csv_path <- file.path(PATH_output, fig_folder_name, 
                        paste0(fig_folder_name, "_all_data_connections_v", v, ".csv"))
  write.csv(all_data, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# Run unique and check #s again
# ------------------------------------------------------------------------------
count_olf <- unique(all_data$orn_id) #321
count_al <- unique(all_data$al_id) #14
count_int <- unique(all_data$int_id) #7
count_nsc <- unique(all_data$nsc_id) #8


# add NT information to all data -----------------------------------------------
# This was used starting from Panel I,J,K

nt_lookup <- neurotransmitters %>%
  select(root_id, nt_type_stringent_corr)

all_data <- all_data %>%
  left_join(nt_lookup %>% rename(orn_id = root_id, 
                                 orn_nt_type_stringent_corr = nt_type_stringent_corr),
            by = "orn_id") %>%
  left_join(nt_lookup %>% rename(al_id = root_id,
                                 al_nt_type_stringent_corr = nt_type_stringent_corr),
            by = "al_id") %>%
  left_join(nt_lookup %>% rename(nsc_id = root_id,
                                 nsc_nt_type_stringent_corr = nt_type_stringent_corr),
            by = "nsc_id")

# ------------------------------------------------------------------------------
# Figure 5B - Numbers for plot
# ------------------------------------------------------------------------------
# ORN behav count
olf_behav <- all_data %>%
  group_by(orn_Behav_Signif) %>%
  summarise(n_count = length(unique(orn_id)))

# AL count
al_count <- all_data %>%
  group_by(al_class, al_subclass) %>%
  summarise(unique_count = length(unique(al_id)),
            .groups = "drop")

# INT count
int_class_sub_hemi_grp <- all_data %>%
  group_by(int_class, int_sub_class, int_hemibrain_type) %>%
  summarise(n_count = length(unique(int_id)),
            .groups = "drop")

# NSC count
nsc_grp <- all_data %>%
  group_by(name_nsc) %>%
  summarise(n_count = length(unique(nsc_id)))

olf_behav <- olf_behav %>% arrange(desc(n_count))
al_count <- al_count %>% arrange(desc(unique_count))
int_class_sub_hemi_grp <- int_class_sub_hemi_grp %>% arrange(desc(n_count))
nsc_grp <- nsc_grp %>% arrange(desc(n_count))

save_dir <- file.path(PATH_output, fig_folder_name)
sink(file.path(save_dir, "Figure_5B.txt"), split = TRUE)

cat("=== ORN behavior significance counts ===\n")
print(olf_behav)

cat("\n=== AL class/subclass counts ===\n")
print(al_count)

cat("\n=== INT class/sub_class/hemibrain_type counts ===\n")
print(int_class_sub_hemi_grp)

cat("\n=== NSC counts ===\n")
print(nsc_grp)
all_counts <- list(
  ORN = olf_behav,
  AL = al_count,
  INT = int_class_sub_hemi_grp,
  NSC = nsc_grp
)

for (name in names(all_counts)) {
  cat("\n===", name, "===\n")
  print(all_counts[[name]])
}

sink() 

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 5C - ORN by behavior and synapses
# ------------------------------------------------------------------------------
panel_name = "C"
filename_orn_behav = paste0(fig_folder_name, panel_name)

# ORN by hemibraintype (PN types) & behav signif
olf_behav_grp <- all_data %>%
  group_by(orn_PN_types, orn_Behav_Signif) %>%
  summarise(n_count = length(unique(orn_id)),
            .groups = "drop") %>%
  arrange(desc(n_count))

# Save to csv
if (write_csv){
  csv_path <- file.path(PATH_output, fig_folder_name, paste0(filename_orn_behav, ".csv"))
  write.csv(olf_behav_grp, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}


# Update the factor levels of Behav_Signif based on the sorted order
olf_behav_grp$orn_PN_types <- factor(olf_behav_grp$orn_PN_types, 
                                     levels = unique(olf_behav_grp$orn_PN_types[order(-olf_behav_grp$n_count)]))
olf_behav_grp$orn_Behav_Signif <- factor(olf_behav_grp$orn_Behav_Signif, 
                                         levels = unique(olf_behav_grp$orn_Behav_Signif[order(-olf_behav_grp$n_count)]))

# Bar chart - count of # of neurons
g <- ggplot(olf_behav_grp, aes(x = fct_reorder(orn_PN_types, desc(n_count)),
                               y = n_count, fill = orn_Behav_Signif)) +
  geom_col() +
  labs(x = "ORN type", y = "# of neurons") +
  scale_fill_manual(values = behav_signif) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1,
                                   color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.title.x = element_text(size = 14, color="black"),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(size = 16, hjust = 0.5, color="black"),
        legend.text = element_text(size = 10, color="black"),   
        legend.title = element_blank(),
        text = element_text(size = 12, color="black")) + 
  scale_y_continuous(breaks = seq(0, 60, by = 10), limits = c(0, 60),
                     expand = c(0,0))

print(g)

if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename_orn_behav,
                                                             ".pdf"))
  ggsave(pdf_path, plot = g, width = 16, height = 8, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 5D - # of synapses by behav
# ------------------------------------------------------------------------------
panel_name = "D"
filename_orn_syn = paste0(fig_folder_name, panel_name)

# Only unique combinations across ids
df_unique <- all_data %>%
  distinct(orn_id, orn_n_synapses, al_id, .keep_all = TRUE)

result <- df_unique %>%
  group_by(orn_Behav_Signif) %>%
  summarise(synapse_sum = sum(orn_n_synapses))

# Sort in descending order
olf_synapse_grp <- result %>%
  arrange(desc(synapse_sum))

# Save to csv
if (write_csv){
  csv_path <- file.path(PATH_output, fig_folder_name, paste0(filename_orn_syn,
                                                             ".csv"))
  write.csv(olf_synapse_grp, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# Update the factor levels of Behav_Signif based on the sorted order
olf_synapse_grp$orn_Behav_Signif <- factor(olf_synapse_grp$orn_Behav_Signif,
                                           levels = unique(olf_synapse_grp$orn_Behav_Signif[order(-olf_synapse_grp$synapse_sum)]))

# Create the plot
b <- ggplot(olf_synapse_grp, aes(x = fct_reorder(orn_Behav_Signif,
                                                 desc(synapse_sum)),
                                 y = synapse_sum, fill = orn_Behav_Signif)) +
  geom_col() +
  labs(y = "# of synapses") +
  scale_fill_manual(values = behav_signif) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        axis.line = element_line(color = "black"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        legend.title = element_blank(), 
        legend.position = "none",
        text = element_text(size = 12, color = "black")) +
  scale_y_continuous(breaks = seq(0, 7500, by = 1500), limits = c(0, 7500),
                     expand = c(0,0))

print(b)

if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename_orn_syn,
                                                             ".pdf"))
  ggsave(pdf_path, plot = b, width = 6, height = 8, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 5E - Top 10 ORN type (IDs for codex)
# ------------------------------------------------------------------------------
panel_name = "E_ORN_top10"
filename_orn_top = paste0(fig_folder_name, panel_name)

input_df <- unique(all_data[all_data$orn_class == 'olfactory',])
top_synapes_orn <- input_df %>%
  group_by(orn_hemibrain_type) %>%
  summarise(sum_synapses = sum(orn_n_synapses))
top10_tmp = top_synapes_orn[order(-top_synapes_orn$sum_synapses),]$orn_hemibrain_type[1:10]
top_input = all_data[all_data$orn_hemibrain_type %in% top10_tmp,]
# Create a data frame with id and their respective order
order_df = data.frame(orn_hemibrain_type = top10_tmp, order = 1:10)
# Merge the order information
top_input = merge(top_input, order_df, by = "orn_hemibrain_type")

# Save to csv
if (write_csv){
  csv_path <- file.path(PATH_output, fig_folder_name, paste0(filename_orn_top,
                                                             ".csv"))
  write.csv(top_input, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# IDs for neuroglancer
for (g in unique(top_input$orn_Behav_Signif)) {
  group_data <- top_input %>% filter(orn_Behav_Signif == g)
  unique_ids <- unique(group_data$orn_id)
  
  cat("\n", g, "- count:", length(unique_ids), "\n")
  cat(paste(unique_ids, collapse = ", "), "\n")
}
# neuroglancer (https://edit.flywire.ai/) for visualization & screenshots
# neuroglancer links provided as rtf's on zenodo
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 5F - AL by glomeruli
# ------------------------------------------------------------------------------
panel_name = "F"
panel_type = "_al_class_sub_hemi_grp"
filename_al = paste0(fig_folder_name, panel_name)
filename_al_csv = paste0(fig_folder_name, panel_name, panel_type)

# al by class, subclass and hemibrain type
al_class_sub_hemi_grp <- all_data %>%
  group_by(al_subclass, al_hemibrain_type) %>%
  summarise(n_count = length(unique(al_id)),
            .groups = "drop")

# assign NA as other
al_class_sub_hemi_grp$al_subclass <- al_class_sub_hemi_grp$al_subclass %>%
                                     replace_na("undefined")

# Save to csv
if (write_csv){
  csv_path <- file.path(PATH_output, fig_folder_name, paste0(filename_al_csv,
                                                             ".csv"))
  write.csv(al_class_sub_hemi_grp, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

al_subgrp <- al_class_sub_hemi_grp %>%
  group_by(al_subclass) %>%
  arrange(al_subclass, desc(n_count)) %>%
  ungroup() %>%
  mutate(al_hemibrain_type = factor(al_hemibrain_type, levels = unique(al_hemibrain_type)))

# Bar chart - count of # of neurons
g <- ggplot(al_subgrp, aes(x = al_hemibrain_type, y = n_count)) +
  geom_col() +
  labs(y = "# of neurons") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5,
                                   color="black"),
        axis.text.y = element_text(size = 12, color="black"),
        axis.title = element_text(size = 14, color="black"),
        axis.title.x = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        plot.title = element_text(size = 16, hjust = 0.5, color="black"),
        legend.text = element_text(size = 10, color="black"),   
        legend.title = element_blank(),
        text = element_text(size = 12, color="black")) + 
  scale_y_continuous(breaks = seq(0, 2, by = 1), limits = c(0, 2), expand = c(0,0))

print(g)

if (write_plots){
  pdf_path <- file.path(PATH_output, fig_folder_name, paste0(filename_al, ".pdf"))
  ggsave(pdf_path, plot = g, width = 10, height = 8, units = "cm")
  cat("Saved PDF to:", normalizePath(pdf_path), "\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 5G - Top ALPN type
# ------------------------------------------------------------------------------
panel_name = "G_ALPN_top10"
filename_alpn = paste0(fig_folder_name, panel_name)

al_class_grp <- all_data %>%
  group_by(al_class, al_hemibrain_type) %>%
  summarise(sum_synapses = sum(al_n_synapses),
            .groups = "drop")

top10_tmp = al_class_grp[order(-al_class_grp$sum_synapses),]$al_hemibrain_type[1:10]
top_input = all_data[all_data$al_hemibrain_type %in% top10_tmp,]
# Create a data frame with id and their respective order
order_df = data.frame(al_hemibrain_type = top10_tmp, order = 1:10)
# Merge the order information
top_input = merge(top_input, order_df, by = "al_hemibrain_type")

# Save to csv
if (write_csv){
  csv_path <- file.path(PATH_output, fig_folder_name, paste0(filename_alpn,
                                                             ".csv"))
  write.csv(top_input, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}

# Top 4 only, for console print
top4_tmp <- top10_tmp[1:4]
top_input_print <- top_input[top_input$al_hemibrain_type %in% top4_tmp, ]

# Print unique al_id per orn_Behav_Signif group (top 4 only)
for (g in unique(top_input_print$orn_Behav_Signif)) {
  group_data <- top_input_print %>% filter(orn_Behav_Signif == g)
  unique_ids <- unique(group_data$al_id)
  
  cat("\n", g, "- count:", length(unique_ids), "\n")
  cat(paste(unique_ids, collapse = ", "), "\n")
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 5H - Top NSC type
# ------------------------------------------------------------------------------
nsc_top_grp <- all_data %>%
  group_by(name_nsc, nsc_id) %>%
  summarise(n_count = n(),
            .groups = "drop")

for (g in unique(all_data$name_nsc)) {
  group_data <- all_data %>% filter(name_nsc == g)
  unique_ids <- unique(group_data$nsc_id)
  
  cat("\n", g, "- count:", length(unique_ids), "\n")
  cat(paste(unique_ids, collapse = ", "), "\n")
}
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Fig 5I - DILP
# ------------------------------------------------------------------------------
# DILP
dilp <- all_data[all_data$name_nsc == "m_NSC_DILP", ]
nsc_uniq <- unique(dilp$nsc_id)#2
int_uniq <- unique(dilp$int_id)#1
al_uniq <- unique(dilp$al_id) #2
orn_uniq <- unique(dilp$orn_id) #13

# DILP only dataset
df_unique <- dilp %>%
  distinct(orn_id, orn_n_synapses, al_id, .keep_all = TRUE)

# Synapse count based on behavior
synp <- df_unique %>%
  group_by(orn_Behav_Signif) %>%
  summarise(synapse_sum = sum(orn_n_synapses))

# Behav count
result <- df_unique %>%
  group_by(orn_Behav_Signif) %>%
  summarise(n_count = n())

al_synp <- dilp %>%
  group_by(al_hemibrain_type, al_cell_type, al_subclass, al_side, int_id,
           al_nt_type_stringent_corr) %>%
  summarise(synapse_sum = sum(al_n_synapses),
            n_count = n(),
            synapse_unique = synapse_sum / n_count,
            .groups = "drop")

int_synp <- dilp %>%
  group_by(int_hemibrain_type, int_id, nsc_id, name_nsc, int_nt_type_stringent_corr) %>%
  summarise(synapse_sum = sum(int_n_synapses),
            n_count = n(),
            synapse_unique = synapse_sum / n_count,
            .groups = "drop")

# ------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Fig 5J - DH31
# ------------------------------------------------------------------------------
dh31 <- all_data[all_data$name_nsc == "l_NSC_DH31", ]
nsc_uniq <- unique(dh31$nsc_id)#2
int_uniq <- unique(dh31$int_id)#2
al_uniq <- unique(dh31$al_id) #6
orn_uniq <- unique(dh31$orn_id) #201

# Only unique combinations across ids
df_unique <- dh31 %>%
  distinct(orn_id, orn_n_synapses, al_id, .keep_all = TRUE)

synp <- df_unique %>%
  group_by(al_id, orn_Behav_Signif, al_hemibrain_type, al_side) %>%
  summarise(synapse_sum = sum(orn_n_synapses),
            .groups = "drop")

# Behav count
df_behav <- df_unique %>%
  distinct(orn_id, .keep_all = TRUE)

result <- df_behav %>%
  group_by(orn_Behav_Signif) %>%
  summarise(n_count = n())

al_synp <- dh31 %>%
  group_by(al_hemibrain_type, al_cell_type, al_subclass, al_side, int_id,
           al_nt_type_stringent_corr) %>%
  summarise(synapse_sum = sum(al_n_synapses),
            n_count = n(),
            synapse_unique = synapse_sum / n_count,
            .groups = "drop")

int_synp <- dh31 %>%
  group_by(int_hemibrain_type, int_id, nsc_id, name_nsc, int_nt_type_stringent_corr) %>%
  summarise(synapse_sum = sum(int_n_synapses),
            n_count = n(),
            synapse_unique = synapse_sum / n_count,
            .groups = "drop")
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Fig 5K - CAPA
# ------------------------------------------------------------------------------
capa <- all_data[all_data$name_nsc == "SEZ_NSC_CAPA", ]
nsc_uniq <- unique(capa$nsc_id)#2
int_uniq <- unique(capa$int_id)#2
al_uniq <- unique(capa$al_id) #4
orn_uniq <- unique(capa$orn_id) #105

## adjust for only unique combinations across ids
df_unique <- capa %>%
  distinct(orn_id, orn_n_synapses, al_id, .keep_all = TRUE)

synp <- df_unique %>%
  group_by(al_id, orn_Behav_Signif) %>%
  summarise(synapse_sum = sum(orn_n_synapses),
            .groups = "drop")

# adjust and get only unique ids
df_behav <- df_unique %>%
  distinct(orn_id, .keep_all = TRUE)

result <- df_behav %>%
  group_by(orn_Behav_Signif) %>%
  summarise(n_count = n())


# AL* names
al_synp <- capa %>%
  group_by(al_hemibrain_type, al_cell_type, al_subclass, al_side, int_id,
           al_nt_type_stringent_corr) %>%
  summarise(synapse_sum = sum(al_n_synapses),
            n_count = n(),
            synapse_unique = synapse_sum / n_count,
            .groups = "drop")


# INT 
int_synp <- capa %>%
  group_by(int_hemibrain_type, int_id, nsc_id, name_nsc, int_nt_type_stringent_corr) %>%
  summarise(synapse_sum = sum(int_n_synapses),
            n_count = n(),
            synapse_unique = synapse_sum / n_count,
            .groups = "drop")
# ------------------------------------------------------------------------------
