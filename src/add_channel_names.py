import argparse
from os.path import join
from bs4 import BeautifulSoup
import pandas as pd

def add_channel_names(in_xml_file, out_xml_file, markers_df):
  
  with open(in_xml_file) as f:
    in_xml = f.read()
  
  soup = BeautifulSoup(in_xml, 'xml')
  
  channel_els = soup.find_all('Channel')
  
  markers_df["channel_name"] = markers_df.apply(lambda row: f"{row['marker_name']} (Cycle {row['cycle_number']})", axis='columns')
  channel_names = markers_df["channel_name"].values.tolist()
  
  for name, el in zip(channel_names, channel_els):
    el['Name'] = name
    
  with open(out_xml_file, "w") as f:
    f.write(str(soup))

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument(
    '-s',
    '--sample_id',
    type=str,
    required=True,
    help='The sample ID, e.g. lung_2_1'
  )
  args = parser.parse_args()
  
  in_xml_file = join("data", args.sample_id, "registration", f"{args.sample_id}.in.ome.xml")
  out_xml_file = join("data", args.sample_id, "registration", f"{args.sample_id}.out.ome.xml")
  
  markers_df = pd.read_csv(join("data", args.sample_id, "markers.csv"))
    
  add_channel_names(
    in_xml_file,
    out_xml_file,
    markers_df
  )