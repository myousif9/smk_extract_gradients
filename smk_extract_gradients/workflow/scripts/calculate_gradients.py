import nibabel as nib
import pandas as pd
import numpy as np

from brainspace.gradient import GradientMaps

# Get number of subjects
n_subjects = len(snakemake.input.correlation_matrix)

# Load subject info
# subject_info = pd.read_table(snakemake.input.subject_info, header=0)

# Load first subject's correlation matrix for shape
correlation_matrix_data = np.load(snakemake.input.correlation_matrix[0])

# Create empty matrix to fill
correlation_matrix_concat = np.zeros((
    correlation_matrix_data.shape[0],
    correlation_matrix_data.shape[1],
    n_subjects)
)

# Fill empty matrix using subjects' correlation matrices
for s in range(0,n_subjects):
    correlation_matrix_concat[:,:,s] = np.load(snakemake.input.correlation_matrix[s])

# Average correlation matrices
correlation_matrix_avg = np.mean(correlation_matrix_concat,2)

# Save matrix
np.save(snakemake.output.avg_correlation_matrix, correlation_matrix_avg)

# Calculate average gradients
gm = GradientMaps(n_components=snakemake.params.n_gradients, kernel='cosine', random_state=0)
gm.fit(correlation_matrix_avg, diffusion_time=0)

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

nib.save(gii, snakemake.output.avg_gradient_maps)

# Calculate each subject gradient
for s in range(0,n_subjects):
    gp = GradientMaps(n_components=snakemake.params.n_gradients, kernel='cosine', alignment='procrustes', random_state=0)
    gp.fit(correlation_matrix[:,:,s}, reference=gm.gradients_, diffusion_time=0)

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


