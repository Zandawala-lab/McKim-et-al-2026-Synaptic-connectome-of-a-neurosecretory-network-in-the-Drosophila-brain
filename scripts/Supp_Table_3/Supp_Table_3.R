# --- Supp Table 3 -------------------------------------------------------------
# ------------------------------------------------------------------------------
# Code recreates supplemental table 3 with NSC ids across datasets
# ------------------------------------------------------------------------------

# load packages ----------------------------------------------------------------
library(tidyverse)

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
fig_folder_name = "Supp_Table_3"
if(!dir.exists(file.path(PATH_output, fig_folder_name))){
  dir.create(file.path(PATH_output, fig_folder_name))
  cat("Created figure folder at:", normalizePath(file.path(PATH_output,
                                                           fig_folder_name)), "\n")
}

# ------------------------------------------------------------------------------
# Options for saving plots & csv files
# ------------------------------------------------------------------------------
# When running script for the first time, set to TRUE:
write_csv   = TRUE           # TRUE - save processed data associated w/figures
                             # FALSE - data not saved outside of R

# load data: -------------------------------------------------------------------
NSC <- read_delim(file.path(PATH_input, paste0("NSC_v", v, ".csv")),
                  delim = ",",
                  escape_double = FALSE,
                  col_types = cols(NSC_id = col_character()),
                  trim_ws = TRUE)

# First generated in script for Figure_1C_S1.R or Figure_3_S3.R
maleCNS_brain_NSC <- read_delim(file.path(PATH_input, "maleCNS_brain_NSC_meta_v0.9.csv"),
                                delim = ",",
                                escape_double = FALSE,
                                col_types = cols(bodyid = col_character()),
                                trim_ws = TRUE)

# First generated in script for Figure_1C_S1.R or Figure_3_S3.R
BANC_brain_NSC <- read_delim(file.path(PATH_input, "BANC_NSC_brain_meta_v888.csv"),
                             delim = ",",
                             escape_double = FALSE,
                             col_types = cols(pt_root_id = col_character()),
                             trim_ws = TRUE)

# ------------------------------------------------------------------------------
# sort and organize table according to NSC order
maleCNS_NSC_sorted <- maleCNS_brain_NSC %>%
  arrange(factor(MckimType, levels = unique(NSC$NSC_name))) %>%
  select(bodyid, MckimType, somaSide)

maleCNS_NSC_sorted <- maleCNS_NSC_sorted %>%
  mutate(somaSide = recode(somaSide, "L" = "left", "R" = "right"))

# adjust naming in BANC for consistency
#unique(BANC_brain_NSC$cell_type)
BANC_brain_NSC <- BANC_brain_NSC %>%
    mutate(cell_type = case_when(
      cell_type == "m-NSC" ~ "m_NSC_unknown",
      cell_type == "l_NSC" ~ "l_NSC_unknown",
      TRUE ~ cell_type
    ))

# ------------------------------------------------------------------------------
# Rename FAFB columns
NSC <- NSC %>%
  rename(
    FAFB_NSC_name = NSC_name,
    FAFB_NSC_id = NSC_id,
    FAFB_NSC_hemisphere = hemisphere
  )

# Rename BANC columns
BANC_NSC <- BANC_brain_NSC %>%
  rename(
    BANC_root_id = pt_root_id,
    BANC_cell_type = cell_type,
    BANC_side = side
  )

# Rename maleCNS columns
maleCNS_NSC <- maleCNS_NSC_sorted %>%
  rename(
    maleCNS_bodyid = bodyid,
    maleCNS_MckimType = MckimType,
    maleCNS_somaSide = somaSide
  )

# ------------------------------------------------------------------------------
# match rows by cell type for datasets with less entries 
pad_df_to_match_nsc <- function(df, type_col, nsc_counts) {
  padded_groups <- list()
  
  for (cell_type in names(nsc_counts)) {
    group <- df[df[[type_col]] == cell_type, ]
    n_missing <- nsc_counts[[cell_type]] - nrow(group)
    
    if (n_missing > 0) {
      empty_rows <- as.data.frame(matrix(NA, nrow = n_missing, ncol = ncol(df)))
      colnames(empty_rows) <- colnames(df)
      group <- rbind(group, empty_rows)
    }
    padded_groups[[cell_type]] <- group
  }
  
  do.call(rbind, padded_groups)
}

# Get per-type counts from NSC (the reference)
nsc_counts <- table(NSC$FAFB_NSC_name)[unique(NSC$FAFB_NSC_name)]  # preserves order

# Pad the other two
BANC_NSC_padded <- pad_df_to_match_nsc(BANC_NSC, "BANC_cell_type", nsc_counts)
maleCNS_NSC_padded <- pad_df_to_match_nsc(maleCNS_NSC, "maleCNS_MckimType", nsc_counts)

# ------------------------------------------------------------------------------
# Combine by columns
NSC_combined <- bind_cols(NSC, BANC_NSC_padded, maleCNS_NSC_padded)

# ------------------------------------------------------------------------------
# Final clean up and org
rownames(NSC_combined) <- NULL

NSC_datasets <- NSC_combined %>%
  select(FAFB_NSC_name, FAFB_NSC_id, BANC_root_id, maleCNS_bodyid)

NSC_datasets <- NSC_datasets %>%
  rename(NSC_name = FAFB_NSC_name,
         FAFB_id = FAFB_NSC_id,
         BANC_id = BANC_root_id,
         maleCNS_id = maleCNS_bodyid)

if (write_csv){
  csv_path <- file.path(PATH_output, fig_folder_name, 
                                  "Supp_Table_3_NSC_Across_Datasets.csv")
  write.csv(NSC_datasets, csv_path)
  cat("Saved CSV to:", normalizePath(csv_path), "\n")
}
