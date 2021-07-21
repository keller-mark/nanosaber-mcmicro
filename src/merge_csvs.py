import os
from os.path import join
import pandas as pd
import numpy as np

def merge_csvs(sample_ids, input_files, output_file):
  
  df = pd.DataFrame(index=[], columns=["sample_id"])
  for sample_id, input_file in zip(sample_ids, input_files):
    sample_df = pd.read_csv(input_file, index_col=0)
    sample_df["sample_id"] = sample_id
    df = df.append(sample_df, ignore_index=True)
  
  df.to_csv(output_file, index=True)
  

if __name__ == "__main__":
  merge_csvs(
    snakemake.params["sample_ids"],
    snakemake.input,
    snakemake.output[0],
  )
  
  