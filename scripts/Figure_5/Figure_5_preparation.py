# -- Figure 5 preparation -----------------------------------------------------
#-------------------------------------------------------------------------------
# Code which has to be run in preparation for Figure 5
#-------------------------------------------------------------------------------
import caveclient
import numpy as np
from scipy import spatial
import pandas as pd
import time
import os

datastack_name = "flywire_fafb_production"
client = caveclient.CAVEclient(datastack_name)

v = str(pd.read_csv("./input/version.csv", sep=",")["version"].values[0])

interneuron = pd.read_csv(("./input/al_interneurons_input_filtered_ids_v"+v+".csv"), sep=",")

interneuron_id = interneuron["pre_pt_root_id"]
syn_df = client.materialize.query_table("synapses_nt_v1", 
                                  filter_in_dict={"post_pt_root_id" : interneuron_id},
                                  materialization_version=v)
syn_id = syn_df["id"]
# Original code - if too much to query all at once, use the batched version below
# valid_syn_df = client.materialize.query_table("valid_synapses_nt_v2", 
#                                   filter_in_dict={"target_id": syn_id},
#                                   merge_reference=False)
#                                   

syn_id_list = syn_id.tolist()

batch_size = 50000
results = []
failed_batches = []

for i in range(0, len(syn_id_list), batch_size):
    batch = syn_id_list[i:i+batch_size]
    
    for attempt in range(3):
        try:
            start = time.time()
            df_batch = client.materialize.query_table(
                "valid_synapses_nt_v2",
                filter_in_dict={"target_id": batch},
                merge_reference=False
            )
            elapsed = time.time() - start
            print(f"Batch {i}-{i+len(batch)}: {len(df_batch)} rows in {elapsed:.1f}s")
            results.append(df_batch)
            break
        except Exception as e:
            print(f"Batch {i}-{i+len(batch)} failed (attempt {attempt+1}): {e}")
            time.sleep(3)
    else:
        failed_batches.append((i, batch))
    
    time.sleep(0.5)

valid_syn_df = pd.concat(results, ignore_index=True)
print(f"\nTotal rows retrieved: {len(valid_syn_df)}")
print(f"Failed batches: {len(failed_batches)}")

valid_id = valid_syn_df["target_id"]
result = syn_df.query("id in @valid_id")

outpath = "./input/tmp/input_to_interneurons_filtered_v"+v+".csv"
result.to_csv(outpath)

print(f"Saved filtered input data to: {os.path.abspath(outpath)}")



