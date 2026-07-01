# -- Figure 3 preparation ---
# ------------------------------------------------------------------------------
# Code to be run in preparation for Figure 3 - get input to NSCs
# ------------------------------------------------------------------------------
import caveclient
import numpy as np
from scipy import spatial
import pandas as pd

datastack_name = "flywire_fafb_production"
client = caveclient.CAVEclient(datastack_name)

v = str(pd.read_csv("./input/version.csv", sep=",")["version"].values[0])

NSC = pd.read_csv(("./input/NSC_v"+v+".csv"), sep=",")

NSC_id = NSC["NSC_id"]
syn_df = client.materialize.query_table("synapses_nt_v1", 
                                  filter_in_dict={"post_pt_root_id" : NSC_id},
                                  materialization_version=v)
syn_id = syn_df["id"]

# Chunk syn_id to avoid size limits
chunk_size = 5000
syn_id_list = syn_id.tolist()
valid_chunks = []
n_chunks = -(-len(syn_id_list) // chunk_size)

for i in range(0, len(syn_id_list), chunk_size):
    chunk = syn_id_list[i:i + chunk_size]
    print(f"Chunk {i // chunk_size + 1} / {n_chunks}")
    valid_chunk_df = client.materialize.query_table(
        "valid_synapses_nt_v2",
        filter_in_dict={"target_id": chunk},
        merge_reference=False
    )
    valid_chunks.append(valid_chunk_df)

valid_syn_df = pd.concat(valid_chunks, ignore_index=True)
valid_id = valid_syn_df["target_id"]
result = syn_df.query("id in @valid_id")

result.to_csv("./input/tmp/NSC_input_filtered_v"+v+".csv")



