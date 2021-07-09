import os
from os.path import join
import pandas as pd
import numpy as np
from scipy.stats import zscore


def filter_quantification(input_file, marker_file, output_file, output_out_file):
  markers_df = pd.read_csv(marker_file)
  quant_df = pd.read_csv(input_file, index_col=0)

  filter_cols = [
    "Area",
  ]

  new_cols = []

  zscore_threshold = 3.0

  quant_out_df = pd.DataFrame()

  for filter_col in filter_cols:
    vals = quant_df[filter_col].values
    quant_df[f"{filter_col}_zscore"] = np.abs(zscore(vals))
    quant_df[f"{filter_col}_outliers"] = (quant_df[f"{filter_col}_zscore"] > zscore_threshold)

    quant_out_df = quant_out_df.append(quant_df.loc[quant_df[f"{filter_col}_outliers"]], ignore_index=False)
    quant_df = quant_df.loc[~quant_df[f"{filter_col}_outliers"]]

    new_cols.append(f"{filter_col}_zscore")
    new_cols.append(f"{filter_col}_outliers")
  
  print(quant_df.head())


  quant_df = quant_df[list(set(quant_df.columns.values.tolist()) - set(new_cols))]
  
  quant_df.to_csv(output_file, index=True)
  quant_out_df.to_csv(output_out_file, index=True)

if __name__ == "__main__":
  filter_quantification(
    snakemake.input["quantification"],
    snakemake.input["markers"],
    snakemake.output["filtered_in"],
    snakemake.output["filtered_out"],
  )
  
  