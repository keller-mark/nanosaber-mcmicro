import os
from os.path import join
import pandas as pd
import numpy as np
from scipy.stats import zscore


def filter_quantification(input_file, marker_file, output_file, output_out_file, sample_id):
  markers_df = pd.read_csv(marker_file)
  quant_df = pd.read_csv(input_file, index_col=0)

  if sample_id == "lung_1_1":
    # Fix the incorrect numbering of the DAPI columns for the lung 1.1
    # image which lacks a DAPI channel for the first cycle.
    quant_df["DAPI_3_cellMask"] = quant_df["DAPI_2_cellMask"]
    quant_df["DAPI_2_cellMask"] = quant_df["DAPI_1_cellMask"]
    quant_df = quant_df.drop(columns=["DAPI_1_cellMask"])

  filter_cols = [
    "Area",
  ]

  quant_df["to_remove"] = False

  new_cols = [ "to_remove" ]

  zscore_threshold = 3.0

  for filter_col in filter_cols:
    vals = quant_df[filter_col].values
    quant_df[f"{filter_col}_zscore"] = np.abs(zscore(vals))
    quant_df[f"{filter_col}_outliers"] = (quant_df[f"{filter_col}_zscore"] > zscore_threshold)

    quant_df["to_remove"] = (quant_df["to_remove"]) | (quant_df[f"{filter_col}_outliers"])

    new_cols.append(f"{filter_col}_zscore")
    new_cols.append(f"{filter_col}_outliers")
  
  # TODO: Filter out cells which had quantification of zero for any DAPI channel

  # TODO: Filter out cells which did not register well (DAPI quantification does not align across cycles)
  
  quant_df = quant_df.loc[~quant_df["to_remove"]]
  quant_out_df = quant_df.loc[quant_df["to_remove"]]

  quant_df = quant_df[list(set(quant_df.columns.values.tolist()) - set(new_cols))]
  
  quant_df.to_csv(output_file, index=True)
  quant_out_df.to_csv(output_out_file, index=True)

if __name__ == "__main__":
  filter_quantification(
    snakemake.input["quantification"],
    snakemake.input["markers"],
    snakemake.output["filtered_in"],
    snakemake.output["filtered_out"],
    snakemake.wildcards["sample_id"]
  )
  
  