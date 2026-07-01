# -- Figure 3 S1 ---------------------------------------------------------------
# ------------------------------------------------------------------------------
# Code to run for Figure 3 S1 to obtain the brain mesh 
# ------------------------------------------------------------------------------

import navis
import cloudvolume
import pandas as pd
import os

# read in input version file
v = str(pd.read_csv("./input/version.csv", sep=",")["version"].values[0])

id = 1
navis.patch_cloudvolume()

vol = cloudvolume.CloudVolume('precomputed://gs://flywire_neuropil_meshes/whole_neuropil/brain_mesh_v3', use_https=True, progress=False)
 
m = vol.mesh.get(id, as_navis=True)

# check dir exists, if not, make them
outdir = "./output/neurons_v"+v+"/brainmesh"
os.makedirs(outdir, exist_ok=True)

# save mesh obj to: output/neurons_v783/brainmesh/ dir
filepath = outdir+"/brainmesh.obj"
navis.write_mesh(m, filepath)
print(f"Saved brainmesh file to: {os.path.abspath(filepath)}")
