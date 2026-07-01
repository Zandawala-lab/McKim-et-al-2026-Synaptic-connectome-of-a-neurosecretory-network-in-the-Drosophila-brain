# -- Figure 9E & S2, S3, S4 -----------------------------------------------------
# ------------------------------------------------------------------------------
# Code for dot plots
# 9E: Neuropeptide receptor expression in adult tissues
# S2: Expression of hormone receptors released from brain NSCs
# S3: Expression of hormone receptors in peripheral tissues
# S4: Expression of hormone receptors in gut and reproductive tract
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Files are large, can start from beginning with download and follow preprocessing
# Email contact author (Meet Zandawala) for link to share tissueseurat.rda
# file needed to replicate analyses
# ------------------------------------------------------------------------------

# load packages ----------------------------------------------------------------
# need Seurat v5.4.0 and SeuratObject v5.3.0 for these analyses

library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(tidyr)
library(ggplot2)
library(readxl)
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
fig_folder_name = "Figure_9"
if(!dir.exists(file.path(PATH_output, fig_folder_name))){
  dir.create(file.path(PATH_output, fig_folder_name))
  cat("Created figure folder at:", normalizePath(file.path(PATH_output, 
                                                           fig_folder_name)), "\n")
}


# ------------------------------------------------------------------------------
# Options for loading data (bc raw files are large)
# ------------------------------------------------------------------------------
# When running script for the first time, use setting provided:
load_data = TRUE
# ------------------------------------------------------------------------------

if (load_data) {
  load(file = file.path(PATH_input, "tissueseurat.rda"))
  
} else {
  # ----------------------------------------------------------------------------
  # Load raw data and create object by processing the .loom and .h5ad files from flycellatlas
  # ----------------------------------------------------------------------------
  all_loom <- Connect(filename = file.path(PATH_input, "s_fca_biohub_all_wo_blood_10X.loom",
                                           mode = 'r'))
  # download stringent data from https://flycellatlas.org/#data
  
  # gene * cell matrix
  all_mat <- all_loom[["/matrix"]][, ]
  all_mat <- Matrix::Matrix(all_mat, sparse = T)
  all_mat <- Matrix::t(all_mat)
  
  # extract the cell id and gene id and set them as the column and row names of the matrix
  all_cellid <- all_loom[["/col_attrs/CellID"]][]
  all_geneid <- all_loom[["/row_attrs/Gene"]][]
  
  colnames(all_mat) <- all_cellid
  rownames(all_mat) <- all_geneid
  
  # figure out metadata stored in the file
  all_loom[["col_attrs"]]
  
  # pull metadata from loom
  attrs <- c(
    'Barcode',
    'CellID',
    'ClusterID',
    'n_counts',
    'n_genes',
    'percent_mito',
    'R_annotation',
    'R_annotation_broad',
    'age',
    'tissue',
    'annotation',
    'annotation_broad',
    'annotation_broad_extrapolated',
    'sex'
  )
  
  # extract attributes from loom
  attrs_all_loom <- map_dfc(attrs, ~ all_loom[[paste0("/col_attrs/", .)]][]) %>% 
                    as.data.frame()
  colnames(attrs_all_loom) <- attrs
  rownames(attrs_all_loom) <- all_cellid
  
  # create seurat object from loom data
  all_seurat <- CreateSeuratObject(counts = all_mat, meta.data = attrs_all_loom)
  View(all_seurat@meta.data)
  
  # download stringent data from https://flycellatlas.org/#data
  # convert h5ad file to h5seurat and load it
  Convert(
    file.path(PATH_input, "s_fca_biohub_all_wo_blood_10x.h5ad"),
    file.path(PATH_input, "all.h5seurat")
  )
  all_H5seurat <- LoadH5Seurat(file.path(PATH_input, "all.h5seurat"))
  
  # use the dimensionality reduction information of h5ad here
  # cell ID of h5ad is different from cell ID of loom - convert it to cell ID of loom
  all_seurat@reductions <- all_H5seurat@reductions
  rownames(all_seurat@reductions$pca@cell.embeddings) <- all_cellid
  rownames(all_seurat@reductions$tsne@cell.embeddings) <- all_cellid
  rownames(all_seurat@reductions$umap@cell.embeddings) <- all_cellid
  
  # check data quality
  VlnPlot(all_seurat,
          features = c("nFeature_RNA", "nCount_RNA"),
          ncol = 2)
  FeatureScatter(all_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  + geom_smooth(method = 'lm')
  
  # Normalize data
  all_seurat <- NormalizeData(all_seurat,
                              normalization.method = "LogNormalize",
                              scale.factor = 10000)
  
  # scale data since loom file did not have scaled data
  all.genes <- rownames(all_seurat)
  all_seurat <- ScaleData(all_seurat, features = all.genes)
  
  # ----------------------------------------------------------------------------
  save(all_seurat, file = file.path(PATH_input, "all_seurat.rda"))
  
  # load(file.path(PATH_input, "all_seurat.rda"))
  # ----------------------------------------------------------------------------
  
  # Subset all_seurat to exclude "head" and "body tissues.
  tissueseurat <- subset(
    all_seurat,
    tissue %in% c(
      "oenocyte",
      "fat_body",
      "haltere",
      "proboscis_and_maxpalp",
      "antenna",
      "trachea",
      "testis",
      "ovary",
      "gut",
      "malpighian_tubule",
      "body_wall",
      "heart",
      "male_reproductive_glands",
      "leg",
      "wing"
    )
  )
  
  # Check the number of unique cell annotations in the new tissueseurat dataset #167
  unique(tissueseurat@meta.data$annotation)
  
  # Load corrected tissue classifications from the manually edited excel file
  # Default tissue to annotation classification in seurat object is incorrect
  
  annotation_new <- read_excel(file.path(PATH_input, "uniqueannotation.xlsx"))
  
  # Change the names of the columns of the imported excel file/dataframe
  colnames(annotation_new)[1] <- "annotation_corr"
  colnames(annotation_new)[2] <- "tissue_corr"
  
  tissueseurat@meta.data$annotation_corr <- tissueseurat@meta.data$annotation
  seurat_metadata <- as.data.frame(tissueseurat@meta.data)
  
  # Convert row names to a column
  seurat_metadata <- seurat_metadata %>% rownames_to_column(var = "rowname")
  
  # Left join
  result <- seurat_metadata %>%
    left_join(annotation_new, by = "annotation_corr")
  
  # Convert the column back to row names
  result <- result %>% column_to_rownames(var = "rowname")
  
  # Update the Seurat object's metadata
  tissueseurat@meta.data <- result
  tissueseurat@meta.data
  
  # set corrected tissue as new idents
  Idents(tissueseurat) <- tissueseurat@meta.data$tissue_corr
  Idents(tissueseurat)
  
  # ----------------------------------------------------------------------------
  # Save the tissueseurat object - used for the featureplots and dotplots
  save(tissueseurat, file = file.path(PATH_input, "tissueseurat.rda"))
  # ----------------------------------------------------------------------------
}

# ------------------------------------------------------------------------------
# Figure 9E - Dotplot by tissue
# ------------------------------------------------------------------------------

hormone_receptors <- c("CrzR", "sNPF-R", "Gyc76C", "TkR86C", "TkR99D", 
                       "Dh31-R", "Dh44-R1", "Dh44-R2", "InR", "CCKLR-17D1", 
                       "CCKLR-17D3", "MsR1", "MsR2", "CapaR", "PK1-R", "PK2-R1",
                       "PK2-R2")

# sorts hormone receptors alphabetically
alpha <- sort(hormone_receptors)

c <- DotPlot(object = tissueseurat, features = alpha, group.by = 'tissue_corr', 
             idents = c('fat_body', 'nervous_system', 'haltere','proboscis_maxpalp',
                        'antenna','trachea','male_reproductive','female_reproductive',
                        'gut','Malpighian_tubule','body_wall','heart','general','leg',
                        'wing','salivary_gland'), dot.min = 0.05, scale = FALSE) +
  xlab(NULL) +
  ylab(NULL) +
  scale_size(range = c(3,8)) + RotatedAxis() + 
  scale_color_gradient(low = "black", high = "red") +
  theme(panel.grid.major = element_line(color = "grey80", linetype = "dashed"),
        panel.grid.minor = element_line(color = "grey90", linetype = "dashed"),
        axis.text.x = element_text(face = "italic")) 
print(c)

outpath <- file.path(PATH_output, fig_folder_name, paste0(fig_folder_name, "E.pdf"))
ggsave(outpath, plot = c, width = 22, height = 15, units = "cm")
cat("Saved PDF to:", normalizePath(outpath), "\n")
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Figure 9 S2 - Featureplots
# ------------------------------------------------------------------------------
# Featureplot showing expression of hormone receptors in tSNE space

# Map of variable/file label -> actual feature name in tissueseurat
receptor_features <- c(
  "CrzR"      = "CrzR",
  "sNPFR"     = "sNPF-R",
  "Gyc76C"    = "Gyc76C",
  "TkR86C"    = "TkR86C",
  "TkR99D"    = "TkR99D",
  "Dh31R"     = "Dh31-R",
  "Dh44R1"    = "Dh44-R1",
  "Dh44R2"    = "Dh44-R2",
  "InR"       = "InR",
  "CCKLR17D1" = "CCKLR-17D1",
  "CCKLR17D3" = "CCKLR-17D3",
  "MsR1"      = "MsR1",
  "MsR2"      = "MsR2",
  "CapaR"     = "CapaR",
  "PK1R"      = "PK1-R",
  "PK2R1"     = "PK2-R1",
  "PK2R2"     = "PK2-R2"
)

featureplots <- list()

for (label in names(receptor_features)) {
  feature_name <- receptor_features[[label]]
  
  p <- FeaturePlot(tissueseurat, features = feature_name, slot = "scale.data", 
                   reduction = "tsne", min.cutoff = "q10", max.cutoff = "q90", 
                   pt.size = 0.1, cols = c("#CCCCCC", "#FF9900", "#FF0000"),
                   raster = FALSE, order = TRUE) & 
    theme(plot.title = element_text(size = 8), legend.text = element_text(size = 7),
          legend.key.size = unit(3, "mm"), 
          panel.background = element_rect(colour = "black", size = 1, fill = NA),
          legend.box.spacing = unit(0.5, "mm")) & 
    NoAxes()
  
  featureplots[[label]] <- p
  
  outpath <- file.path(PATH_output, fig_folder_name, paste0(fig_folder_name,
                                                            "_S2_featureplot_",
                                                            label, ".pdf"))
  ggsave(outpath, plot = p, width = 10, height = 10, units = "cm")
  cat("Saved PDF to:", normalizePath(outpath), "\n")
}
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Dimplot showing tissue clusters in tSNE space
# Legend for the featureplot figure
# Alter the color of the dimplot to make each tissue more visible in tSNE space

dimplot <- DimPlot(tissueseurat, 
                   reduction = 'tsne', 
                   group.by = 'tissue_corr', 
                   pt.size = 0.1,
                   label = F,
                   cols = c(
                     "antenna"             = "#79baf2",
                     "body_wall"           = "#0dfa9e",
                     "exclude"             = "#fc9a99",
                     "fat_body"            = "#8c4613",
                     "female_reproductive" = "#b0c4dd",
                     "gut"                 = "#80ff00",
                     "haltere"             = "#0000fe",
                     "heart"               = "#2ff18b",
                     "leg"                 = "#028002",
                     "male_reproductive"   = "#fea500",
                     "Malpighian_tubule"   = "#bdb76b",
                     "nervous_system"      = "#ee82ef",
                     "proboscis_maxpalp"   = "#fe0000",
                     "salivary_gland"      = "#0bf2fc",
                     "trachea"             = "#ac8f14",
                     "wing"                = "#ffff3b"
                   ),
                   repel = TRUE,
                   raster = FALSE) + 
  labs(title = NULL)
print(dimplot)

outpath <- file.path(PATH_output, fig_folder_name, paste0(fig_folder_name,
                                                          "_S2_dimplot_tissueseurat.pdf"))
ggsave(outpath, plot = dimplot, width = 20, height = 15, units = "cm")
cat("Saved PDF to:", normalizePath(outpath), "\n")
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Figure 9 S3 & S4
# Dotplots showing hormone receptor expression in different tissues 
# ------------------------------------------------------------------------------
dotplot_specs <- list(
  list(label = "fatbodybodywall_dotplot",          sfig = "S3", idents = c('fat_body', 'body_wall'),               width = 22, height = 8),
  list(label = "legwinghaltere_dotplot",           sfig = "S3", idents = c('leg', 'wing', 'haltere'),              width = 22, height = 8),
  list(label = "salivaryglandprobmaxpalp_dotplot", sfig = "S3", idents = c('proboscis_maxpalp', 'salivary_gland'), width = 22, height = 10),
  list(label = "tracheageneral_dotplot",           sfig = "S3", idents = c('trachea', 'general'),                  width = 22, height = 10),
  list(label = "antenna_dotplot",                  sfig = "S3", idents = 'antenna',                                width = 22, height = 20),
  list(label = "MT_dotplot",                       sfig = "S3", idents = 'Malpighian_tubule',                      width = 22, height = 10),
  list(label = "heart_dotplot",                    sfig = "S3", idents = 'heart',                                  width = 22, height = 10),
  list(label = "nervoussystem_dotplot",            sfig = "S3", idents = 'nervous_system',                         width = 22, height = 7),
  list(label = "malereproductive_dotplot",         sfig = "S4", idents = 'male_reproductive',                      width = 22, height = 32),
  list(label = "femalereproductive_dotplot",       sfig = "S4", idents = 'female_reproductive',                    width = 22, height = 25),
  list(label = "gut_dotplot",                      sfig = "S4", idents = 'gut',                                    width = 22, height = 24)
)

dotplots <- list()

for (spec in dotplot_specs) {
  p <- DotPlot(object = tissueseurat, features = alpha, group.by = 'annotation', 
               idents = spec$idents, dot.min = 0.05, scale = FALSE) +
    xlab(NULL) +
    ylab(NULL) +
    scale_size(range = c(3, 8)) + RotatedAxis() + 
    scale_color_gradient(low = "black", high = "red") +
    theme(panel.grid.major = element_line(color = "grey80", linetype = "dashed"),
          panel.grid.minor = element_line(color = "grey90", linetype = "dashed"),
          axis.text.x = element_text(face = "italic"))
  
  dotplots[[spec$label]] <- p
  
  outpath <- file.path(PATH_output, fig_folder_name, paste0(fig_folder_name, "_",
                                                            spec$sfig, "_",
                                                            spec$label,
                                                            ".pdf"))
  ggsave(outpath, plot = p, width = spec$width, height = spec$height, units = "cm")
  cat("Saved PDF to:", normalizePath(outpath), "\n")
}
# ------------------------------------------------------------------------------