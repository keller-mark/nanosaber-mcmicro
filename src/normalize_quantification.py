import os
from os.path import join
import pandas as pd
import numpy as np
from scipy.stats import iqr


def normalize_quantification(sample_ids, input_files, marker_files, output_files, to_iqr=False):
  
  markers_df = pd.DataFrame(index=[], columns=["sample_id"])
  for sample_id, marker_file in zip(sample_ids, marker_files):
    df = pd.read_csv(marker_file)
    df["sample_id"] = sample_id
    markers_df = markers_df.append(df, ignore_index=True)
    
  sample_id_to_input_file = dict(zip(sample_ids, input_files))
  sample_id_to_output_file = dict(zip(sample_ids, output_files))
  
  marker_values = {}
  marker_extents = {}
  
  sample_id_to_df = {}
  for sample_id in sample_ids:
    sample_id_to_df[sample_id] = pd.read_csv(sample_id_to_input_file[sample_id], index_col=0)
  
  # Compute the extent for each marker
  for marker_id, marker_df in markers_df.groupby("marker_name"):
    for sample_id, sample_marker_df in marker_df.groupby("sample_id"):
      cycle_nums = sample_marker_df["cycle_number"].astype(int).tolist()
            
      for cycle_num in cycle_nums:
        if sample_id == "lung_1_1":
          cycle_num -= 1
        quant_col = f"{marker_id}_{cycle_num}_cellMask" if marker_id == "DAPI" else f"{marker_id}_cellMask"
        
        quant_df = sample_id_to_df[sample_id]
        #print(quant_df.columns)
        quant_vals = quant_df[quant_col].values
        
        if marker_id not in marker_values:
          marker_values[marker_id] = quant_vals
        else:
          marker_values[marker_id] = np.concatenate((marker_values[marker_id], quant_vals), axis=None)
    
    if to_iqr:
      marker_extents[marker_id] = (np.percentile(marker_values[marker_id], 25), np.percentile(marker_values[marker_id], 75))
    else:
      marker_extents[marker_id] = (np.amin(marker_values[marker_id]), np.amax(marker_values[marker_id]))
    
  
  # Normalize the sample quantification dataframe for each marker column
  for marker_id, marker_df in markers_df.groupby("marker_name"):
    for sample_id, sample_marker_df in marker_df.groupby("sample_id"):
      cycle_nums = sample_marker_df["cycle_number"].astype(int).tolist()
      
      print(sample_id, marker_id, cycle_nums)
      
      for cycle_num in cycle_nums:
        if sample_id == "lung_1_1":
          cycle_num -= 1
        quant_col = f"{marker_id}_{cycle_num}_cellMask" if marker_id == "DAPI" else f"{marker_id}_cellMask"
        
        quant_df = sample_id_to_df[sample_id].copy()
        quant_vals = quant_df[quant_col].values

        if to_iqr:
          sample_25_pct = np.percentile(quant_vals, 25)
          sample_75_pct = np.percentile(quant_vals, 75)
          sample_range = (sample_75_pct - sample_25_pct)

          full_25_pct = marker_extents[marker_id][0]
          full_75_pct = marker_extents[marker_id][1]
          full_range = full_75_pct - full_25_pct
          
          # Normalize the column
          if sample_range > 0:
            quant_df[quant_col] = (quant_df[quant_col] * (full_range/sample_range)) - (sample_25_pct - full_25_pct)

        else:
          sample_min = np.amin(quant_vals)
          sample_max = np.amax(quant_vals)
          sample_range = (sample_max - sample_min)
          
          full_min = marker_extents[marker_id][0]
          full_max = marker_extents[marker_id][1]
          full_range = full_max - full_min


          
          # Normalize the column
          if sample_range > 0:
            quant_df[quant_col] = (quant_df[quant_col] * (full_range/sample_range)) - (sample_min - full_min)
          
        # Store the updated dataframe
        sample_id_to_df[sample_id] = quant_df
  
  # Save each normalized sample df to file
  for sample_id, output_file in zip(sample_ids, output_files):
    #sample_id_to_df[sample_id]["CellID"] = sample_id_to_df[sample_id]["CellID"].astype(int)
    sample_id_to_df[sample_id].to_csv(output_file, index=True)

if __name__ == "__main__":
  normalize_quantification(
    snakemake.params["sample_ids"],
    snakemake.input["quantification"],
    snakemake.input["markers"],
    snakemake.output,
  )
  
  