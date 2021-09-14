import nibabel as nib
import pandas as pd
import numpy as np
import math

from nilearn.connectome import ConnectivityMeasure
from brainspace.gradient import GradientMaps
from scipy.spatial.distance import pdist, squareform


# Load affinity matrix
norm_angle_matrix = np.load(snakemake.input.affinity_matrix)

# Calculate gradients based on normalized angle matrix
# Kernel = none as input matrix is already affinity matrix
gm = GradientMaps(n_components=snakemake.params.n_gradients, kernel=None, random_state=0)
gm.fit(norm_angle_matrix, diffusion_time=0)

# Save gradients to gifti file
gii = nib.gifti.GiftiImage()

for g in range(0,len(gm.gradients_.T)):
    gii.add_gifti_data_array(
        nib.gifti.GiftiDataArray(
            data=gm.gradients_.T[g].astype(np.float32),
            meta={
                'AnatomicalStructurePrimary':'CortexLeft' if snakemake.wildcards.hemi is 'L' else 'CortexRight',
                'Name':'Gradient {}'.format(g+1)
                }
            )
    ) 

nib.save(gii, snakemake.output.gradient_maps)