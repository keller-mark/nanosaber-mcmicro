import os
import platform
from os.path import join
import pandas as pd

O2_USER = "mk596"
# Check if this is running on O2
IS_O2 = (platform.system() == "Linux" and O2_USER != None)

# Directory / file constants
DATA_DIR = ("data" if not IS_O2 else join(os.sep, "n", "scratch3", "users", O2_USER[0], O2_USER, "lr", "data"))

FIJI_CMD = "java -jar -Xmx4096m /Applications/Fiji.app/jars/ij-1.53c.jar -ijpath /Applications/Fiji.app/ -batch"

SAMPLE_IDS = [
  "lung_1_1",
  "lung_2_1",
  "lung_2_2"
]

def get_markers_df(sample_id):
  return pd.read_csv(join(DATA_DIR, sample_id, "markers.csv"))
  
print(get_markers_df("lung_2_1"))

rule all:
  input:
    expand(
      join(DATA_DIR, "{sample_id}", "pre-raw", "{sample_id}-cycle-01-channel-01.tif"),
      sample_id=SAMPLE_IDS
    ),
    expand(
      join(DATA_DIR, "{sample_id}", "bunwarpj", "{sample_id}-cycle-01-channel-01_direct_transf.txt"),
      sample_id=SAMPLE_IDS
    ),
    expand(
      join(DATA_DIR, "{sample_id}", "raw", "{sample_id}-cycle-01.tif"),
      sample_id=SAMPLE_IDS
    ),
    expand(
      join(DATA_DIR, "{sample_id}", "quantification", "unmicst-{sample_id}.csv"),
      sample_id=SAMPLE_IDS
    ),
    expand(
      join(DATA_DIR, "{sample_id}", "quantification", "unmicst-{sample_id}.normalized.csv"),
      sample_id=SAMPLE_IDS
    ),
    expand(
      join(DATA_DIR, "{sample_id}", "quantification", "unmicst-{sample_id}.normalized.fcs"),
      sample_id=SAMPLE_IDS
    )


rule fcs_gating:
  input:
    expand(
      join(DATA_DIR, "{sample_id}", "quantification", "unmicst-{sample_id}.normalized.fcs"),
      sample_id=SAMPLE_IDS
    )
  output:
    join(DATA_DIR, "{sample_id}", "quantification", "unmicst-{sample_id}.gating.png")
  script:
    join("src", "fcs_gating.py")

rule csv_to_fcs:
  input:
    join(DATA_DIR, "{sample_id}", "quantification", "unmicst-{sample_id}.normalized.csv")
  output:
    join(DATA_DIR, "{sample_id}", "quantification", "unmicst-{sample_id}.normalized.fcs")
  script:
    join("src", "csv_to_fcs.R")


rule normalize_quantification:
  input:
    quantification=expand(
      join(DATA_DIR, "{sample_id}", "quantification", "unmicst-{sample_id}.csv"),
      sample_id=SAMPLE_IDS
    ),
    markers=expand(
      join(DATA_DIR, "{sample_id}", "markers.csv"),
      sample_id=SAMPLE_IDS
    ),
  params:
    sample_ids=SAMPLE_IDS
  output:
    expand(
      join(DATA_DIR, "{sample_id}", "quantification", "unmicst-{sample_id}.normalized.csv"),
      sample_id=SAMPLE_IDS
    )
  script:
    join("src", "normalize_quantification.py")
    
    
rule mcmicro:
  input:
    join(DATA_DIR, "{sample_id}", "registration", "{sample_id}.ome.tif")
  output:
    join(DATA_DIR, "{sample_id}", "quantification", "unmicst-{sample_id}.csv")
  params:
    seg_channel_index=(lambda w: 3 if w.sample_id == "lung_1_1" else 4)
  shell:
    """
    nextflow run labsyspharm/mcmicro \
      --in ./data/{wildcards.sample_id} \
      --start-at probability-maps \
      --unmicst-opts '--scalingFactor 1.7 --channel {params.seg_channel_index}' \
      --s3seg-opts '--logSigma 0 15' \
      -w ./work
    """

rule concatenate_stacked_registered:
  input:
    join(DATA_DIR, "{sample_id}", "raw", "{sample_id}-cycle-01.tif"),
    join(DATA_DIR, "{sample_id}", "raw", "{sample_id}-cycle-02.tif"),
    join(DATA_DIR, "{sample_id}", "raw", "{sample_id}-cycle-03.tif")
  output:
    join(DATA_DIR, "{sample_id}", "registration", "{sample_id}.ome.tif")
  shell:
    """
    {FIJI_CMD} ./src/concatenate_registered.js {wildcards.sample_id}
    """

rule stack_registered:
  input:
    join(DATA_DIR, "{sample_id}", "pre-raw", "{sample_id}-cycle-01-channel-01-registered.tif"),
    join(DATA_DIR, "{sample_id}", "pre-raw", "{sample_id}-cycle-02-channel-01.tif"),
    join(DATA_DIR, "{sample_id}", "pre-raw", "{sample_id}-cycle-03-channel-01-registered.tif")
  output:
    join(DATA_DIR, "{sample_id}", "raw", "{sample_id}-cycle-01.tif"),
    join(DATA_DIR, "{sample_id}", "raw", "{sample_id}-cycle-02.tif"),
    join(DATA_DIR, "{sample_id}", "raw", "{sample_id}-cycle-03.tif")
  shell:
    """
    {FIJI_CMD} ./src/stack_registered.js {wildcards.sample_id}__cycle-01 && \
    {FIJI_CMD} ./src/stack_registered.js {wildcards.sample_id}__cycle-02 && \
    {FIJI_CMD} ./src/stack_registered.js {wildcards.sample_id}__cycle-03
    """

rule register:
  input:
    join(DATA_DIR, "{sample_id}", "bunwarpj", "{sample_id}-cycle-01-channel-01_direct_transf.txt"),
    join(DATA_DIR, "{sample_id}", "bunwarpj", "{sample_id}-cycle-03-channel-01_direct_transf.txt"),
    join(DATA_DIR, "{sample_id}", "pre-raw", "{sample_id}-cycle-01-channel-01.tif"),
    join(DATA_DIR, "{sample_id}", "pre-raw", "{sample_id}-cycle-02-channel-01.tif"),
    join(DATA_DIR, "{sample_id}", "pre-raw", "{sample_id}-cycle-03-channel-01.tif")
  output:
    join(DATA_DIR, "{sample_id}", "pre-raw", "{sample_id}-cycle-01-channel-01-registered.tif"),
    join(DATA_DIR, "{sample_id}", "pre-raw", "{sample_id}-cycle-03-channel-01-registered.tif")
  shell:
    """
    bash ./src/register.sh {wildcards.sample_id}
    """

rule compute_registration:
  input:
    join(DATA_DIR, "{sample_id}", "pre-raw", "{sample_id}-cycle-01-channel-01.tif"),
    join(DATA_DIR, "{sample_id}", "pre-raw", "{sample_id}-cycle-02-channel-01.tif"),
    join(DATA_DIR, "{sample_id}", "pre-raw", "{sample_id}-cycle-03-channel-01.tif")
  output:
    join(DATA_DIR, "{sample_id}", "bunwarpj", "{sample_id}-cycle-01-channel-01_direct_transf.txt"),
    join(DATA_DIR, "{sample_id}", "bunwarpj", "{sample_id}-cycle-03-channel-01_direct_transf.txt")
  shell:
    """
    {FIJI_CMD} ./src/compute_registration.ijm {wildcards.sample_id}
    """

rule rename_raw_files:
  output:
    join(DATA_DIR, "{sample_id}", "pre-raw", "{sample_id}-cycle-01-channel-01.tif"),
    join(DATA_DIR, "{sample_id}", "pre-raw", "{sample_id}-cycle-02-channel-01.tif"),
    join(DATA_DIR, "{sample_id}", "pre-raw", "{sample_id}-cycle-03-channel-01.tif")
  params:
    raw_sample_dir=(lambda w: "Tumor tissue_" + w.sample_id.split("_")[1] + "." + w.sample_id.split("_")[2])
  shell:
    """
    bash ./src/rename.sh {wildcards.sample_id} "./data/raw/{params.raw_sample_dir}" ./data/{wildcards.sample_id}/pre-raw
    """