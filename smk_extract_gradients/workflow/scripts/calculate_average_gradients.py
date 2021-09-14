import nibabel as nib
import pandas as pd
import numpy as np

from brainspace.gradient import GradientMaps

# Get number of subjects
n_subjects = len(snakemake.input.affinity_matrix)

# Load subject info
subject_info = pd.read_table(snakemake.input.subject_info, header=None)

# Load first subject's affinity matrix for shape
affinity_matrix_data = np.load(snakemake.input.affinity_matrix[0])

# Create empty matrix to fill
affinity_matrix_concat = np.zeros((
    affinity_matrix_data.shape[0],
    affinity_matrix_data.shape[1],
    n_subjects)
)

# Fill empty matrix using subjects' affinity matrices
for s in range(0,n_subjects):
    affinity_matrix_concat[:,:,s] = np.load(snakemake.input.affinity_matrix[s])

# Average affinity matrices
affinity_matrix_avg = np.mean(affinity_matrix_concat,2)

# Save matrix
np.save(snakemake.output.affinity_matrix, affinity_matrix_avg)

# Calculate gradients based on average normalized angle matrix
# Kernel = none as input matrix is already affinity matrix
gm = GradientMaps(n_components=snakemake.params.n_gradients, kernel=None, random_state=0)
gm.fit(affinity_matrix_avg, diffusion_time=0)

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