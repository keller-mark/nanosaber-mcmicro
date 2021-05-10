# https://www.researchgate.net/post/Does_anyone_know_how_to_convert_a_csv_file_to_a_fcs_file

library(flowCore)

# Read the csv
exprs_data <- read.csv(file.path("~", "research", "dbmi", "nanosaber-mcmicro", "data", "lung_2_1", "quantification", "unmicst-lung.csv"), header=TRUE, sep=",", stringsAsFactors = F, row.names = 1)

# you need to prepare some metadata

column_names = dimnames(exprs_data)[[2]]

meta <- data.frame(name = column_names, desc = column_names)
meta$range <- apply(apply(exprs_data[,column_names], 2, range), 2, diff)
meta$minRange <- apply(exprs_data[,column_names], 2, min)
meta$maxRange <- apply(exprs_data[,column_names], 2, max)
head(meta)

# all these are required for the following steps to work
# a flowFrame is the internal representation of a FCS file
matrix_exprs_data <- data.matrix(exprs_data[,column_names])
ff <- flowCore::flowFrame(exprs = matrix_exprs_data, parameters = AnnotatedDataFrame(meta))
# now you can save it back to the filesystem
flowCore::write.FCS(ff, paste('lung_2_1', 'fcs', sep='.'))