# -- Figure 3 S1 NSC meshes ----------------------------------------------------
# ------------------------------------------------------------------------------
# Code to obtain the NSC meshes for Figure 3 S1
# ------------------------------------------------------------------------------

import navis
import cloudvolume
import pandas as pd
import os

# read in input version file
v = str(pd.read_csv("./input/version.csv", sep=",")["version"].values[0])

# flywire.set_default_dataset("public")

# read in NSC data
NSC = pd.read_csv("./input/NSC_v"+v+".csv", sep=",", dtype=str)
NSC_id = NSC["NSC_id"]

navis.patch_cloudvolume()

vol = cloudvolume.CloudVolume('graphene://https://prodv1.flywire-daf.com/segmentation/1.0/fly_v31',
use_https=True, progress=False)
 
m = vol.mesh.get(NSC_id, as_navis=True)

# check dir exists, if not, make it
outdir = "./output/neurons_v"+v+"/NSC_mesh"
os.makedirs(outdir, exist_ok=True)

navis.write_mesh(m, outdir, filetype = "obj")

print(f"Saved {len(m)} NSC mesh files to: {os.path.abspath(outdir)}")

# clean up the Python environment, as otherwise variables may interfere with the R environment
del NSC
del NSC_id
del vol
del m

