# https://www.researchgate.net/post/Does_anyone_know_how_to_convert_a_csv_file_to_a_fcs_file

library(Biobase)
library(flowCore)

mcmicro_quantification_csv_to_fcs <- function(sample_id) {
  sample_dir <- file.path("~", "research", "dbmi", "nanosaber-mcmicro", "data", sample_id)
  
  # Read the CSV from mcmicro.
  exprs_data <- read.csv(file.path(sample_dir, "quantification", paste0("unmicst-", sample_id, ".csv")), header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)
  
  # Prepare the column name metadata.
  column_names = dimnames(exprs_data)[[2]]
  
  meta <- data.frame(name = column_names, desc = column_names)
  meta$range <- apply(apply(exprs_data[,column_names], 2, range), 2, diff)
  meta$minRange <- apply(exprs_data[,column_names], 2, min)
  meta$maxRange <- apply(exprs_data[,column_names], 2, max)
  head(meta)
  
  # flowFrame is the internal representation of a FCS file.
  matrix_exprs_data <- data.matrix(exprs_data[,column_names])
  ff <- flowCore::flowFrame(exprs = matrix_exprs_data, parameters = AnnotatedDataFrame(meta))
  # Save the flowFrame to an FCS file.
  dir.create(file.path(sample_dir, "flowcore"), showWarnings = FALSE)
  flowCore::write.FCS(ff, file.path(sample_dir, "flowcore", paste(sample_id, 'fcs', sep='.')))
}

mcmicro_quantification_csv_to_fcs("lung_1_1")
mcmicro_quantification_csv_to_fcs("lung_2_1")
mcmicro_quantification_csv_to_fcs("lung_2_2")


