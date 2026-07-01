# -- Figure 9 &  S1 ------------------------------------------------------------
# ------------------------------------------------------------------------------
# Code for paracrine connectivity analyses
# Does not contain Figure 9E (separate script w/other single cell analyses)
# ------------------------------------------------------------------------------

# load packages ----------------------------------------------------------------
library(Seurat)
# must be version v4..4.0
library(SeuratDisk)
# must be version 0.0.0.9021
library(SeuratObject)
# must be version 4.1.4
library(ggplot2)
library(readr)
library(tidyverse)
library(circlize)

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
fig_folder_name = "Figure_9"
if(!dir.exists(file.path(PATH_output, fig_folder_name))){
  dir.create(file.path(PATH_output, fig_folder_name))
  cat("Created figure folder at:", normalizePath(file.path(PATH_output, fig_folder_name)), "\n")
}

# set colors -------------------------------------------------------------------
hormone_color_scale <- c(
  "Ms"      = "#FE7E00",
  "Hug"     = "#BD00B0",
  "Crz"     = "#8100FF",
  "sNPF"    = "#800080",
  "Dh31"    = "#00B4FF",
  "Tk"      = "#32CD32",
  "Dh44"    = "#ffe200",
  "Ilp_Avg" = "#BD0023"
)

track_color_scale <- c(
  "m_NSC_DILP"          = "#BD0023",
  "m_NSC_DH44"          = "#ffe200",
  "m_NSC_DMS"           = "#FE7E00",
  "l_NSC_CRZ"           = "#8100FF",
  "l_NSC_ITP"           = "#838383",
  "l_NSC_DH31"          = "#00B4FF",
  "SEZ_NSC_CAPA"        = "#003BBD",
  "SEZ_NSC_Hugin"       = "#BD00B0",
  "SEZ_NSC_CAPA_Ms"     = "#FE7E00",
  "m_NSC_DMS_Ms"        = "#FE7E00",
  "SEZ_NSC_Hugin_Hug"   = "#BD00B0",
  "l_NSC_CRZ_Crz"       = "#8100FF",
  "l_NSC_CRZ_sNPF"      = "#800080",
  "l_NSC_ITP_sNPF"      = "#800080",
  "l_NSC_ITP_Tk"        = "#32CD32",
  "l_NSC_DH31_Dh31"     = "#00B4FF",
  "m_NSC_DH44_Dh44"     = "#ffe200",
  "m_NSC_DILP_Ilp_Avg"  = "#BD0023"
)


# Define the color scale for Source_cluster
source_color_scale <- c(
  "m_NSC_DILP"    = "#BD0023",
  "m_NSC_DH44"    = "#ffe200ff",
  "m_NSC_DMS"     = "#FE7E00",
  "l_NSC_CRZ"     = "#8100FF",
  "l_NSC_ITP"     = "#838383",
  "l_NSC_DH31"    = "#00B4FF",
  "SEZ_NSC_CAPA"  = "#003BBD",
  "SEZ_NSC_Hugin" = "#BD00B0"
)


# Create a mapping from original Source_cluster to new labels
label_mapping <- c(
  "m_NSC_DILP"    = "DILP",
  "m_NSC_DH44"    = "DH44",
  "m_NSC_DMS"     = "DMS",
  "l_NSC_CRZ"     = "CRZ",
  "l_NSC_ITP"     = "ITP",
  "l_NSC_DH31"    = "DH31",
  "SEZ_NSC_CAPA"  = "CAPA",
  "SEZ_NSC_Hugin" = "Hugin"
)

# ------------------------------------------------------------------------------
# Options for loading data (bc raw files are large)
# ------------------------------------------------------------------------------
# When running script for the first time, use setting provided:
load_data = TRUE

# load processed data ----------------------------------------------------------
# File provided
load(file.path(PATH_input,"akhcc_stringent_seurat.rda"))
# ------------------------------------------------------------------------------


if (load_data) {
  load(file = file.path(PATH_input, "NSC_scRNA.rda"))
} else {
  # ----------------------------------------------------------------------------
  # Load raw data and process all
  # ----------------------------------------------------------------------------
  
  Brain_loom <- Connect(
    filename = file.path(PATH_input, "Aerts_Fly_AdultBrain_Filtered_57k.loom"),
    mode = 'r'
  )
  # download the loom file from scope (https://scope.aertslab.org/#/Davie_et_al_Cell_2018/Davie_et_al_Cell_2018%2FAerts_Fly_AdultBrain_Filtered_57k.loom/gene)
  Brain_loom
  Brain_loom[["/matrix"]]
  
  # gene * cell matrix
  Brain_mat <- Brain_loom[["/matrix"]][, ]
  Brain_mat <- Matrix::Matrix(Brain_mat, sparse = T)
  Brain_mat <- Matrix::t(Brain_mat)
  
  # extract the cell id and gene id and set them as the column and row names of the matrix
  Brain_cellid <- Brain_loom[["/col_attrs/CellID"]][]
  Brain_geneid <- Brain_loom[["/row_attrs/Gene"]][]
  colnames(Brain_mat) <- Brain_cellid
  rownames(Brain_mat) <- Brain_geneid
  
  # create seurat object from loom data
  Brain_seurat <- CreateSeuratObject(counts = Brain_mat, project = "adultbrain")
  View(Brain_seurat@meta.data)
  
  # check data quality
  VlnPlot(
    Brain_seurat,
    features = c("nFeature_RNA", "nCount_RNA"),
    ncol = 2
  )
  FeatureScatter(Brain_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
                + geom_smooth(method = 'lm')
  
  # normalize data
  Brain_seurat <- NormalizeData(Brain_seurat,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000)
  
  # scale data
  all.genes <- rownames(Brain_seurat)
  Brain_seurat <- ScaleData(Brain_seurat, features = all.genes)
  
  # ----------------------------------------------------------------------------
  # save entire dataset - large (7 GB)
  # saveRDS(Brain_seurat, file = file.path(PATH_input,"Brain_seurat.rds"))
  
  # load from here
  # Brain_seurat <- readRDS(file = file.path(PATH_input,"Brain_seurat.rds"))
  # ----------------------------------------------------------------------------
  
  # extract m_NSC_DILP
  m_NSC_DILP <- subset(x = Brain_seurat, subset = Ilp2 > 3 &
                         Ilp3 > 3 & Ilp5 > 3 & ChAT == 0)
  Cells(m_NSC_DILP)
  m_NSC_DILP[["new_cluster"]] <- "m_NSC_DILP"
  
  # extract m_NSC_DH44
  m_NSC_DH44 <- subset(x = Brain_seurat,
                       subset = Dh44 > 2 & CG13248 > 0 & CG13743 > 0 & Lkr > 0 & Phm > 0)
  Cells(m_NSC_DH44)
  m_NSC_DH44[["new_cluster"]] <- "m_NSC_DH44"
  
  # extract l_NSC_CRZ
  l_NSC_CRZ <- subset(
    x = Brain_seurat,
    subset = Crz > 3 &
      sNPF > 3 & Dh44 == 0 & ITP == 0 & ChAT == 0 & Gr64a == 0 & Phm > 0
  )
  Cells(l_NSC_CRZ)
  l_NSC_CRZ[["new_cluster"]] <- "l_NSC_CRZ"
  
  # extract l_NSC_DH31
  l_NSC_DH31 <- subset(x = Brain_seurat, subset = ITP > 2 &
                         Dh31 > 4 & amon > 0 & Phm > 0)
  Cells(l_NSC_DH31)
  l_NSC_DH31[["new_cluster"]] <- "l_NSC_DH31"
  
  # extract l_NSC_ITP
  l_NSC_ITP <- subset(x = Brain_seurat,
                      subset = Tk > 1 & sNPF > 1 & ITP > 1 & ImpL2 > 1 & Crz == 0)
  Cells(l_NSC_ITP)
  l_NSC_ITP[["new_cluster"]] <- "l_NSC_ITP"
  
  # extract m_NSC_DMS
  m_NSC_DMS <- subset(x = Brain_seurat, subset = Ms > 2 &
                        EcR > 0 & rk > 0 & amon > 0 & Phm > 0)
  Cells(m_NSC_DMS)
  m_NSC_DMS[["new_cluster"]] <- "m_NSC_DMS"
  
  # extract SEZ_NSC_CAPA
  SEZ_NSC_CAPA <- subset(x = Brain_seurat,
                         subset = Capa > 0 & CrzR > 0 & trp > 0 & amon > 0 & Phm > 0)
  Cells(SEZ_NSC_CAPA)
  SEZ_NSC_CAPA[["new_cluster"]] <- "SEZ_NSC_CAPA"
  
  # extract SEZ_NSC_Hugin
  SEZ_NSC_Hugin <- subset(x = Brain_seurat,
                          subset = Hug > 3 & `Dh44-R2` > 0 & amon > 0 & Phm > 0)
  Cells(SEZ_NSC_Hugin)
  SEZ_NSC_Hugin[["new_cluster"]] <- "SEZ_NSC_Hugin"
  
  # merge all NSC subtypes
  NSC <- merge(
    m_NSC_DILP,
    y = c(
      m_NSC_DH44,
      m_NSC_DMS,
      l_NSC_CRZ,
      l_NSC_DH31,
      l_NSC_ITP,
      SEZ_NSC_CAPA,
      SEZ_NSC_Hugin
    ),
    add.cell.ids = c("DILP", "DH44", "DMS", "CRZ", "DH31", "ITP", "CAPA", "Hugin"),
    merge.data = TRUE ,
    project = "NSC"
  )
  
  # scale data
  all.genes <- rownames(NSC)
  NSC <- ScaleData(NSC, features = all.genes)
  
  save(NSC, file = file.path(PATH_input,"NSC_scRNA.rda"))
}

# ------------------------------------------------------------------------------
# Figure 9 A - Dotplot of NSC markers
# ------------------------------------------------------------------------------
# create a named vector of custom y-labels with superscripts
custom_y_labels <- c("l-NSC^CRZ", "l-NSC^DH31", "l-NSC^ITP", "m-NSC^DH44", 
                     "m-NSC^DILP", "m-NSC^DMS", "SEZ-NSC^CAPA", "SEZ-NSC^Hugin")

NSC_markers <- c("amon","svr","Pal2","Phm","Cadps", "Crz", "sNPF", "Dh31", "ITP", 
                 "Tk", "Lk", "ImpL2", "Lkr", "Dh44", "Ilp2", "Ilp3", "Ilp5", "Ms",
                 "EcR", "rk", "Capa", "trp", "CrzR", "Hug", "Dh44-R2")

exp_mat <- as.matrix(NSC[["RNA"]]@data[NSC_markers,])

# cell meta data
meta <- NSC@meta.data %>% 
  select(new_cluster)
# merge in the expression data
meta <- bind_cols(meta, as.data.frame(t(exp_mat)))
# convert to tidy format
meta <- pivot_longer(meta, -new_cluster, names_to="Gene", values_to="Expression")
meta_summary <- meta %>%
  group_by(new_cluster, Gene) %>%
  summarise(Avg = mean(Expression),
            Pct = sum(Expression > 0) / length(Expression) * 100)
meta_summary$Gene <- factor(meta_summary$Gene, levels=NSC_markers)

# filter data points with less than 5% average expression
meta_summary_filtered <- meta_summary[meta_summary$Pct >= 5, ]

# scaled dotplot
r <- ggplot(meta_summary_filtered, aes(x=Gene, y=new_cluster)) +
  geom_point(aes(size = Pct, fill = Avg), color="black", shape=21) +
  scale_size("% expressed", range = c(1,8)) +
  scale_fill_gradient2(
    name = "Scaled expression",
    low = "grey",
    mid = "orange",
    high = "red",
    midpoint = 3.8, 
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill") +
  ylab(NULL) + xlab(NULL) +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "grey80", linetype = "dashed"),
        panel.grid.minor = element_line(color = "grey90", linetype = "dashed"),
        axis.text.x = element_text(face = "italic", size=12, angle=45, hjust=1, 
                                   color="black"),
        axis.text.y = element_text(size=12, color="black")) +
  scale_y_discrete(labels = parse(text = custom_y_labels))
print(r)

outpath <- file.path(PATH_output, fig_folder_name, paste0(fig_folder_name, "A.pdf"))
ggsave(outpath, plot = r, width = 27, height = 11, units = "cm")
cat("Saved PDF to:", normalizePath(outpath), "\n")

# ------------------------------------------------------------------------------
# Figure 9 B - Dotplot of hormone receptors in NSCs
# ------------------------------------------------------------------------------
hormone_receptors <- c("CapaR", "CCKLR-17D1", "CCKLR-17D3","CrzR","Dh31-R",
                       "Dh44-R1", "Dh44-R2","Gyc76C", "InR", "Lkr", "MsR1", 
                       "MsR2", "PK1-R", "PK2-R1", "PK2-R2", "sNPF-R", "TkR86C", 
                       "TkR99D")

exp_mat <- as.matrix(NSC[["RNA"]]@data[hormone_receptors,])

# cell meta data
meta <- NSC@meta.data %>% 
  select(new_cluster)
# merge in the expression data
meta <- bind_cols(meta, as.data.frame(t(exp_mat)))
# convert to tidy format
meta <- pivot_longer(meta, -new_cluster, names_to="Gene", values_to="Expression")
meta_summary <- meta %>%
  group_by(new_cluster, Gene) %>%
  summarise(Avg = mean(Expression),
            Pct = sum(Expression > 0) / length(Expression) * 100,
            .groups = "drop")
meta_summary$Gene <- factor(meta_summary$Gene, levels=hormone_receptors)

# filter data points with less than 5% average expression
meta_summary_filtered <- meta_summary[meta_summary$Pct >= 5, ]

outpath <- file.path(PATH_output, fig_folder_name, "receptor_expression_scaled.csv")
write.csv(meta_summary_filtered, file = outpath)
cat("Saved CSV to:", normalizePath(outpath), "\n")

# scaled dotplot
s <- ggplot(meta_summary_filtered, aes(x=Gene, y=new_cluster)) +
  geom_point(aes(size = Pct, fill = Avg), color="black", shape=21) +
  scale_size("% expressed", range = c(1,8)) +
  scale_fill_gradient2(
    name = "Scaled expression",
    low = "grey",
    mid = "orange",
    high = "red",
    midpoint = 0.65,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill") +
  ylab(NULL) + xlab(NULL) +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "grey80", linetype = "dashed"),
        panel.grid.minor = element_line(color = "grey90", linetype = "dashed"),
        axis.text.x = element_text(face = "italic", size=12, angle=45, hjust=1,
                                   color="black"),
        axis.text.y = element_text(size=12, color="black")) +
  scale_y_discrete(labels = parse(text = custom_y_labels))
print(s)

outpath <- file.path(PATH_output, fig_folder_name, paste0(fig_folder_name, "B.pdf"))
ggsave(outpath, plot = s, width = 27, height = 11, units = "cm")
cat("Saved PDF to:", normalizePath(outpath), "\n")

# ------------------------------------------------------------------------------
# Figure 9 S1 I - Dotplot of filtered neuropeptide expression
# ------------------------------------------------------------------------------
# Filter neuropeptides based on thresholding
# Prepare data for hormones (scaled data) for plotting
hormones <- c("Crz", "sNPF", "Dh31", "ITP", "Lk", "Tk","Dh44", 
              "Ilp2", "Ilp3", "Ilp5","Dsk","Ms", "Capa", "Hug")
exp_mat <- as.matrix(NSC[["RNA"]]@data[hormones,])

# cell meta data
meta <- NSC@meta.data %>% 
  select(new_cluster)
# merge in the expression data
meta <- bind_cols(meta, as.data.frame(t(exp_mat)))
# convert to tidy format
meta <- pivot_longer(meta, -new_cluster, names_to="Gene", values_to="Expression")
meta_summary <- meta %>%
  group_by(new_cluster, Gene) %>%
  summarise(Avg = mean(Expression),
            Pct = sum(Expression > 0) / length(Expression) * 100,
            .groups = "drop")
meta_summary$Gene <- factor(meta_summary$Gene, levels=hormones)

# filter data points with less than 50% average expression
meta_summary_filtered <- meta_summary[meta_summary$Pct >= 50, ]
meta_summary_filtered$expression = meta_summary_filtered$Avg * meta_summary_filtered$Pct
meta_summary_filtered$score = meta_summary_filtered$expression / 
                              max(meta_summary_filtered$expression)
# additional filtering
meta_summary_filtered <- meta_summary_filtered[meta_summary_filtered$score >= 0.25, ]
meta_summary_scaled <- subset(meta_summary_filtered, select = c(new_cluster, Gene, score))

t <- ggplot(meta_summary_scaled, aes(x=Gene, y=new_cluster)) +
  geom_point(aes(size = 4, fill = score), color="black", shape=21) +
  scale_size(NULL, range = c(1,8), labels = NULL) +
  scale_fill_gradient2(
    name = "Scaled expression",
    low = "grey",
    mid = "orange",
    high = "red",
    midpoint = 0.5,
    space = "Lab",
    guide = "colourbar") +
  ylab(NULL) + xlab(NULL) +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "grey80", linetype = "dashed"),
        panel.grid.minor = element_line(color = "grey90", linetype = "dashed"),
        axis.text.x = element_text(face = "italic", size=12, angle=45, hjust=1, 
                                   color="black"),
        axis.text.y = element_text(size=12, color="black")) +
  scale_y_discrete(labels = parse(text = custom_y_labels))
print(t)

outpath <- file.path(PATH_output, fig_folder_name, paste0(fig_folder_name, "_S1I.pdf"))
ggsave(outpath, plot = t, width = 27, height = 11, units = "cm")
cat("Saved PDF to:", normalizePath(outpath), "\n")
# ------------------------------------------------------------------------------

# Average DILP expression and combine into one row
meta_summary_scaled$source <- paste0(meta_summary_scaled$new_cluster,"_",meta_summary_scaled$Gene)
# Define Ilp sources
specific_sources <- c("m_NSC_DILP_Ilp2", "m_NSC_DILP_Ilp3", "m_NSC_DILP_Ilp5")
# Filter the data to get the scores for the specific sources
specific_scores <- meta_summary_scaled %>% filter(source %in% specific_sources) %>% select(score)
# Calculate the average score
average_score <- mean(specific_scores$score)
# Create a new data frame excluding Ilp sources
hormone_scaled <- meta_summary_scaled %>% filter(!source %in% specific_sources)
# Add a new row with the averaged Ilp score
hormone_scaled <- rbind(hormone_scaled, data.frame(new_cluster = "m_NSC_DILP", 
                                                   Gene = "Ilp_Avg", 
                                                   source = "m_NSC_DILP_Ilp_Avg", 
                                                   score = average_score))
# Rename columns for clarity
hormone_scaled <- hormone_scaled %>%
  rename(
    Source_cluster = new_cluster,
    Hormone = Gene,
    Hormone_score = score
  )

# ------------------------------------------------------------------------------
# Prepare data for hormone receptors (scaled data) for plotting
hormone_receptors <- c("CrzR", "sNPF-R", "Dh31-R", "Gyc76C", "TkR86C", "TkR99D", 
                       "Dh44-R1", "Dh44-R2", "InR", "MsR1", 
                       "MsR2", "CapaR", "PK1-R", "PK2-R1", "PK2-R2")

# receptors not included as no peptide present
# c("Lkr","CCKLR-17D1","CCKLR-17D3")

exp_mat <- as.matrix(NSC[["RNA"]]@data[hormone_receptors,])

# Cell meta data
meta <- NSC@meta.data %>% 
  select(new_cluster)
# Merge in the expression data
meta <- bind_cols(meta, as.data.frame(t(exp_mat)))
# Convert to tidy format
meta <- pivot_longer(meta, -new_cluster, names_to="Gene", values_to="Expression")
meta_summary <- meta %>%
  group_by(new_cluster, Gene) %>%
  summarise(Avg = mean(Expression),
            Pct = sum(Expression > 0) / length(Expression) * 100,
            .groups = "drop")
meta_summary$Gene <- factor(meta_summary$Gene, levels=hormone_receptors)

# Filter data points with less than 5% average expression
meta_summary_filtered <- meta_summary[meta_summary$Pct >= 5, ]
meta_summary_filtered$expression = meta_summary_filtered$Avg * meta_summary_filtered$Pct
meta_summary_filtered$score = meta_summary_filtered$expression / max(meta_summary_filtered$expression)

# Add a column with hormones for corresponding receptor
meta_summary_filtered <- meta_summary_filtered %>%
  mutate(Hormone = case_when(
    grepl("sNPF-R", Gene) ~ "sNPF",
    grepl("CrzR", Gene) ~ "Crz",
    grepl("Dh31-R", Gene) ~ "Dh31",
    grepl("Dh44-R1", Gene) ~ "Dh44",
    grepl("Dh44-R2", Gene) ~ "Dh44",
    grepl("MsR1", Gene) ~ "Ms",
    grepl("MsR2", Gene) ~ "Ms",
    grepl("TkR86C", Gene) ~ "Tk",
    grepl("TkR99D", Gene) ~ "Tk",
    grepl("PK2-R1", Gene) ~ "Hug",
    grepl("PK2-R2", Gene) ~ "Hug",
    grepl("InR", Gene) ~ "Ilp_Avg",
    TRUE ~ "Unknown"
  ))

meta_summary_filtered$connection <- paste0(meta_summary_filtered$new_cluster,"_",meta_summary_filtered$Hormone)
columns_to_select <- c("new_cluster","Gene","score","Hormone","connection")
receptor_scaled <- meta_summary_filtered[, columns_to_select]

# Retain only a single unique connection with the higher score value. 
receptor_unique <- receptor_scaled %>%
  group_by(connection) %>%
  filter(score == max(score)) %>%
  ungroup()  # Remove grouping to return to a normal data frame

# Rename columns for clarity
receptor_unique <- receptor_unique %>%
  rename(
    Receptor_cluster = new_cluster,
    Receptor = Gene,
    Receptor_score = score
  )

# Create a combined score matrix for plotting using receptor and hormone scores
# Perform a join operation between receptor_unique and hormone_scaled based on the Hormone column
joined_data <- receptor_unique %>%
  inner_join(hormone_scaled, by = "Hormone", relationship = "many-to-many")

result <- joined_data %>%
  mutate(connection_score = Receptor_score * Hormone_score) %>%
  select(source, Receptor_cluster, connection_score, Hormone, Source_cluster, connection)

# ------------------------------------------------------------------------------
# Figure 9 C - Plot chord diagram showing the connections
# ------------------------------------------------------------------------------

# Sort the result data frame by source
result <- result %>% arrange(source)

# Map the colors to the Hormone values
result$Hormone_color <- hormone_color_scale[result$Hormone]

# Map Source_cluster to colors in the result data frame
result$source_color <- source_color_scale[result$Source_cluster]

outpath <- file.path(PATH_output, fig_folder_name, paste0(fig_folder_name, "C.pdf"))
pdf(file = outpath, width = 13, height = 10)

# Create the chord diagram with specified colors
chordDiagram(result, 
             transparency = 0.3,
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.1),
             grid.col = track_color_scale,
             col = result$Hormone_color,
             symmetric = TRUE,
             directional = 1,
             direction.type = c("arrows"),
             link.arr.type = "big.arrow")

# Add the cluster annotation track with custom labels
circos.track(ylim = c(0, 1.5), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  cluster = result$Source_cluster[which(result$source == sector.index)][1]
  custom_label <- label_mapping[cluster]
  circos.rect(get.cell.meta.data("xlim")[1], 0, get.cell.meta.data("xlim")[2], 1,
              col = source_color_scale[cluster], border = NA, track.index = 1)
  circos.text(mean(get.cell.meta.data("xlim")), 8.3, custom_label,
              facing = "clockwise", niceFacing = TRUE, cex = 0.8, col = "white",
              font = 2)
}, track.height = 0.05, bg.border = NA)

# Add legend for sectors
par(xpd = TRUE) # Allow drawing outside the plot area
legend(x = 0.8, y = 0.8, 
       legend = names(hormone_color_scale), 
       fill = unname(hormone_color_scale), 
       bty = "n", 
       horiz = FALSE,
       cex = 0.9,
       title = "Hormone")

dev.off()
cat("Saved PDF to:", normalizePath(outpath), "\n")

# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Figure 9 S1 - A thru H - chord diagrams 
# ------------------------------------------------------------------------------

# Panel letter, hormone name (must match result$Hormone), 
# and filename label
chord_specs <- list(
  list(panel = "A", hormone = "Ms",      label = "Ms"),
  list(panel = "B", hormone = "Hug",     label = "Hugin"),
  list(panel = "C", hormone = "Crz",     label = "Crz"),
  list(panel = "D", hormone = "sNPF",    label = "sNPF"),
  list(panel = "E", hormone = "Dh31",    label = "Dh31"),
  list(panel = "F", hormone = "Tk",      label = "Tk"),
  list(panel = "G", hormone = "Dh44",    label = "Dh44"),
  list(panel = "H", hormone = "Ilp_Avg", label = "Ilp")
)

for (spec in chord_specs) {
  
  outpath <- file.path(PATH_output, fig_folder_name, 
                       paste0(fig_folder_name, "_S1", spec$panel, 
                              "_", spec$label, "_connectivity.pdf"))
  pdf(file = outpath, width = 13, height = 10)
  
  
  # Define a new color scale where all Hormone_colors except the current hormone are transparent
  result$Hormone_color_modified <- ifelse(result$Hormone == spec$hormone, 
                                          hormone_color_scale[result$Hormone], 
                                          "transparent")
  # Create the chord diagram with specified colors
  chordDiagram(result, 
               transparency = 0.3,
               annotationTrack = "grid",
               preAllocateTracks = list(track.height = 0.1),
               grid.col = track_color_scale,
               col = result$Hormone_color_modified,
               symmetric = TRUE,
               directional = 1,
               direction.type = c("arrows"),
               link.arr.type = "big.arrow")
  # Add the cluster annotation track with custom labels
  circos.track(ylim = c(0, 1.5), panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    cluster = result$Source_cluster[which(result$source == sector.index)][1]
    custom_label <- label_mapping[cluster]
    circos.rect(get.cell.meta.data("xlim")[1], 0, get.cell.meta.data("xlim")[2], 1,
                col = source_color_scale[cluster], border = NA, track.index = 1)
    circos.text(mean(get.cell.meta.data("xlim")), 8.3, custom_label,
                facing = "clockwise", niceFacing = TRUE, cex = 0.8, col = "white",
                font = 2)
  }, track.height = 0.05, bg.border = NA)
  
  dev.off()
  cat("Saved PDF to:", normalizePath(outpath), "\n")
}
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Figure 9 D - Dotplot of hormone receptors expressed in akh cells 
# ------------------------------------------------------------------------------

exp_mat <- as.matrix(akhcc_stringent[["RNA"]]@data[hormone_receptors,])

# Cell meta data
meta <- akhcc_stringent@meta.data %>% 
  select(sex)
# Merge in the expression data
meta <- bind_cols(meta, as.data.frame(t(exp_mat)))
# Convert to tidy format
meta <- pivot_longer(meta, -sex, names_to="Gene", values_to="Expression")
meta_summary <- meta %>%
  group_by(sex, Gene) %>%
  summarise(Avg = mean(Expression),
            Pct = sum(Expression > 0) / length(Expression) * 100,
            .groups = "drop")
meta_summary$Gene <- factor(meta_summary$Gene, levels=hormone_receptors)

# Filter data points with less than 5% average expression
meta_summary_filtered <- meta_summary[meta_summary$Pct >= 5, ]

# scaled dotplot
z <- ggplot(meta_summary_filtered, aes(x=Gene, y=sex)) +
  geom_point(aes(size = Pct, fill = Avg), color="black", shape=21) +
  scale_size("% expressed", range = c(1,8)) +
  scale_fill_gradient2(
    name = "Scaled expression",
    low = "grey",
    mid = "orange",
    high = "red",
    midpoint = 0.65,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill") +
  ylab(NULL) + xlab(NULL) +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "grey80", linetype = "dashed"),
        panel.grid.minor = element_line(color = "grey90", linetype = "dashed"),
        axis.text.y = element_text(face = "italic", size=12, hjust=1, color="black"),
        axis.text.x = element_text(size=12, color="black", angle=90)) +
  coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
print(z)

outpath <- file.path(PATH_output, fig_folder_name, paste0(fig_folder_name, "D.pdf"))
ggsave(outpath, plot = z, width = 8, height = 10, units = "cm")
cat("Saved PDF to:", normalizePath(outpath), "\n")
# ------------------------------------------------------------------------------