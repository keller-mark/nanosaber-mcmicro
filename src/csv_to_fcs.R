# https://www.researchgate.net/post/Does_anyone_know_how_to_convert_a_csv_file_to_a_fcs_file
# https://bioconductor.org/packages/release/bioc/html/flowCore.html
# https://www.rdocumentation.org/packages/flowCore/versions/1.38.2/topics/write.FCS

library(Biobase)
library(flowCore)

# Read the CSV from mcmicro.
exprs_data <- read.csv(snakemake@input[[1]], header = TRUE, sep = ",", stringsAsFactors = FALSE, row.names = 1)


# Prepare the column name metadata.
column_names = dimnames(exprs_data)[[2]]

meta <- data.frame(name = column_names, desc = column_names)
meta$range <- apply(apply(exprs_data[,column_names], 2, range), 2, diff)
meta$minRange <- apply(exprs_data[,column_names], 2, min)
meta$maxRange <- apply(exprs_data[,column_names], 2, max)

# flowFrame is the internal representation of a FCS file.
matrix_exprs_data <- data.matrix(exprs_data[,column_names])
ff <- flowCore::flowFrame(exprs = matrix_exprs_data, parameters = AnnotatedDataFrame(meta))
# Save the flowFrame to an FCS file.
flowCore::write.FCS(ff, snakemake@output[[1]])

