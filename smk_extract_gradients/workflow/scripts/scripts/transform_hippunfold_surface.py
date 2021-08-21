import nibabel as nib
import numpy as np
import os

# Read in hippunfold surface
surf_in = nib.load(snakemake.input.surf)

# Extract coordinates
for darray in [0,1]:
    if isinstance(surf_in.darrays[darray].data[0,0],np.float32):
        surf_in_coords = surf_in.darrays[darray].data
        idx = darray

# Find eventual nan vertices and change to 0
mask = np.all(np.isnan(surf_in_coords), axis=1)

surf_out = surf_in
surf_out.darrays[idx].data[mask] = 0

# Save to file
nib.save(surf_out, snakemake.output.surf)

# Apply transform
xfm_command = "module load connectome-workbench ; wb_command -surface-apply-warpfield {} {} {} -fnirt {}".format(
    snakemake.output.surf, snakemake.input.xfm_rev, snakemake.output.surf, snakemake.input.xfm_fwd)
 
res = os.system(xfm_command)

# Read in transformed surface
surf_out = nib.load(snakemake.output.surf)
surf_out.darrays[idx].data[mask] = np.nan

# Save to file
nib.save(surf_out, snakemake.output.surf)