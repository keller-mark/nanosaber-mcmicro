library(Giotto)
library(Matrix)

sample_id <- "lung_2_1"

sample_dir <- file.path("~", "research", "dbmi", "nanosaber-mcmicro", "data", sample_id)
quant_data <- read.csv(file.path(sample_dir, "quantification", paste0("unmicst-", sample_id, ".csv")), header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)

quant_cols <- colnames(quant_data)
cellmask_cols <- quant_cols[grep("_cellMask$", quant_cols)]

raw_exprs <- t(as.matrix(quant_data[,cellmask_cols]))

spatial_locs <- quant_data[, c("X_centroid", "Y_centroid")]

gobject <- createGiottoObject(
  raw_exprs = raw_exprs,
  spatial_locs = spatial_locs,
  norm_expr = raw_exprs,
)

gobject <- createSpatialNetwork(
  gobject,
  name = "spatial_network",
  dimensions = "all",
  method = "kNN",
  maximum_distance_delaunay = "auto",
  options = "Pp",
  Y = TRUE,
  j = TRUE,
  S = 0,
  minimum_k = 0,
  knn_method = "dbscan",
  k = 4,
  maximum_distance_knn = NULL,
  verbose = F,
  return_gobject = TRUE,
)

spatCorObject <- detectSpatialCorGenes(
  gobject,
  method = "network",
  expression_values = "normalized",
  subset_genes = NULL,
  spatial_network_name = "spatial_network",
  network_smoothing = NULL,
  min_cells_per_grid = 4,
  cor_method = "pearson"
)

heatmSpatialCorGenes(
  gobject,
  spatCorObject,
  use_clus_name = NULL,
  show_cluster_annot = TRUE,
  show_row_dend = T,
  show_column_dend = F,
  show_row_names = T,
  show_column_names = T,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "heatmSpatialCorGenes",
)



showNetworks(gobject, verbose = TRUE)
spatPlot(gobject)
