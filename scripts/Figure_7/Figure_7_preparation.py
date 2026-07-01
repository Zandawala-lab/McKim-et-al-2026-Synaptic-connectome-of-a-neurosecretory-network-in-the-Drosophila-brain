# -- Figure 7  preparation -----------------------------------------------------
# ------------------------------------------------------------------------------
# Code to be run in preparation for Figure 7 - get synaptic output of NSCs
# ------------------------------------------------------------------------------
import caveclient
import numpy as np
from scipy import spatial
import pandas as pd
import os

datastack_name = "flywire_fafb_production"
client = caveclient.CAVEclient(datastack_name)

v = str(pd.read_csv("./input/version.csv", sep=",")["version"].values[0])

NSC = pd.read_csv(("./input/NSC_v"+v+".csv"), sep=",")

NSC_id = NSC["NSC_id"]
syn_df = client.materialize.query_table("synapses_nt_v1", 
                                  filter_in_dict={"pre_pt_root_id" : NSC_id},
                                  materialization_version=v)
syn_id = syn_df["id"]
valid_syn_df = client.materialize.query_table("valid_synapses_nt_v2", 
                                  filter_in_dict={"target_id": syn_id},
                                  merge_reference=False)
valid_id = valid_syn_df["target_id"]
result = syn_df.query("id in @valid_id")

outpath = "./input/tmp/NSC_output_filtered_v"+v+".csv"
result.to_csv(outpath)

print(f"Saved filtered NSC output data to: {os.path.abspath(outpath)}")
