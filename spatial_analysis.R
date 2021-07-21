library(Giotto)
library(Matrix)

sample_id <- "lung_2_1"
sample_name <- "Lung 2.1"

sample_dir <- file.path("~", "research", "dbmi", "nanosaber-mcmicro", "data", sample_id)
quant_data <- read.csv(file.path(sample_dir, "quantification", paste0("unmicst-", sample_id, ".csv")), header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)

cell_class_data <- read.csv(file.path(sample_dir, "flowcore", paste0(sample_id, ".giotto_cell_classes.csv")), header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)

quant_cols <- colnames(quant_data)
cellmask_cols <- quant_cols[grep("_cellMask$", quant_cols)]

raw_exprs <- t(as.matrix(quant_data[,cellmask_cols]))
rownames(raw_exprs) <- gsub("_cellMask", "", rownames(raw_exprs))

spatial_locs <- quant_data[, c("X_centroid", "Y_centroid")]

gobject <- createGiottoObject(
  raw_exprs = raw_exprs,
  spatial_locs = spatial_locs,
  norm_expr = raw_exprs,
)

gobject <- addCellMetadata(
  gobject,
  cell_class_data[, c("name")],
  vector_name = "fc_class"
)


# detectSpatialCorGenes not working with the "grid" method
#######

# gobject <- createSpatialGrid(
#   gobject,
#   name = "spatial_grid",
#   method = "default",
#   sdimx_stepsize = 1,
#   sdimy_stepsize = 1,
#   sdimz_stepsize = NULL,
#   minimum_padding = 1,
#   return_gobject = TRUE
# )
# 
# annotateSpatialGrid(
#   gobject,
#   spatial_grid_name = "spatial_grid",
#   cluster_columns = "fc_class"
# )
# 
# spatPlot(gobject, cell_color = "fc_class", point_size = 2)
# 
# gridSpatCorObj <- detectSpatialCorGenes(
#   gobject,
#   method = "grid",
#   expression_values = "normalized",
#   spatial_grid_name = "spatial_grid",
#   min_cells_per_grid = 4,
#   cor_method = "pearson"
# )




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

rankSpatialCorGroups(
  gobject,
  spatCorObject,
  use_clus_name = "name",
  show_plot = NA,
  return_plot = FALSE,
  save_plot = NA,
  save_param = list(),
  default_save_name = "rankSpatialCorGroups"
)



spatCorObject <- clusterSpatialCorGenes(
  spatCorObject,
  name = "spat_clus",
  hclust_method = "ward.D",
  k = 10,
  return_obj = TRUE
)


rankSpatialCorGroups(
  gobject,
  spatCorObject,
  use_clus_name = "spat_clus",
  show_plot = NA,
  return_plot = FALSE,
  save_plot = NA,
  save_param = list(),
  default_save_name = "rankSpatialCorGroups"
)



gobject <- createSpatialDelaunayNetwork(
  gobject,
  method = "deldir",
  dimensions = "all",
  name = "Delaunay_network",
  maximum_distance = "auto",
  minimum_k = 0,
  options = "Pp",
  Y = TRUE,
  j = TRUE,
  S = 0,
  verbose = T,
  return_gobject = TRUE,
)
plotStatDelaunayNetwork(
  gobject,
  method = "deldir",
  dimensions = "all",
  maximum_distance = "auto",
  minimum_k = 0,
  options = "Pp",
  Y = TRUE,
  j = TRUE,
  S = 0,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "plotStatDelaunayNetwork",
)

CPscore <- cellProximityEnrichment(
  gobject,
  spatial_network_name = "spatial_network",
  cluster_column = "fc_class",
  number_of_simulations = 1000,
  adjust_method = "none",
  set_seed = TRUE,
  seed_number = 1234
)
cellProximityBarplot(
  gobject,
  CPscore,
  min_orig_ints = 5,
  min_sim_ints = 5,
  p_val = 0.05,
  show_plot = NA,
  return_plot = NA,
  save_plot = NA,
  save_param = list(),
  default_save_name = "cellProximityBarplot"
)







