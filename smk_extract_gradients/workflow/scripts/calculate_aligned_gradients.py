import nibabel as nib
import pandas as pd
import numpy as np
import math

from brainspace.gradient import GradientMaps

# Load affinity matrix
norm_angle_matrix = np.load(snakemake.input.affinity_matrix)

# Load average affinity matrix
avg_norm_angle_matrix = np.load(snakemake.input.avg_affinity_matrix)

# Calculate gradients based on normalized angle matrix
# Kernel = none as input matrix is already affinity matrix
gp = GradientMaps(n_components=snakemake.params.n_gradients, kernel=None, alignment='procrustes', random_state=0)
gp.fit([avg_norm_angle_matrix, norm_angle_matrix], diffusion_time=0)

# Save aligned gradients to gifti file
gii = nib.gifti.GiftiImage()
gii_unaligned = nib.gifti.GiftiImage()

for g in range(0,snakemake.params.n_gradients):
    gii.add_gifti_data_array(
        nib.gifti.GiftiDataArray(
            data=gp.aligned_[1][:,g].astype(np.float32),
            meta={
                'AnatomicalStructurePrimary':'CortexLeft' if snakemake.wildcards.hemi is 'L' else 'CortexRight',
                'Name':'Gradient {}'.format(g+1)
                }
            )
    ) 

    gii_unaligned.add_gifti_data_array(
        nib.gifti.GiftiDataArray(
            data=gp.gradients_[1][:,g].astype(np.float32),
            meta={
                'AnatomicalStructurePrimary':'CortexLeft' if snakemake.wildcards.hemi is 'L' else 'CortexRight',
                'Name':'Gradient {}'.format(g+1)
                }
            )
    ) 

nib.save(gii, snakemake.output.gradient_maps)
nib.save(gii_unaligned, snakemake.output.gradient_unaligned)