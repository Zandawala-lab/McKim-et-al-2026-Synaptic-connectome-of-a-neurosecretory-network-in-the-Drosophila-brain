# -- Figure_2_S4_preparation ---
# ------------------------------------------------------------------------------
# Code to run in preparation for Figure 2_S4 - gets nuclei volume data
# ------------------------------------------------------------------------------
import caveclient
import pandas as pd

datastack_name = "flywire_fafb_public"
client = caveclient.CAVEclient(datastack_name)

#check public versions available
client.materialize.get_versions()

#check which annotation tables are available
client.materialize.get_tables()

v = str(pd.read_csv("./input/version.csv", sep=",")["version"].values[0])

nuclei_df = client.materialize.query_table("nuclei_v1", materialization_version=v)
nuclei_df

nuclei_df.to_csv("./input/tmp/nuclei_volume_v"+v+".csv")

#https://github.com/seung-lab/FlyConnectome/blob/main/CAVE%20tutorial.ipynb
#check here on how to filter the data



