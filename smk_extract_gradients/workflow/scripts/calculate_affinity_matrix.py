import nibabel as nib
import pandas as pd
import numpy as np
import math

from nilearn.connectome import ConnectivityMeasure
from brainspace.gradient import GradientMaps
from scipy.spatial.distance import pdist, squareform

def pull_data(cifti,structure):
    cifti_file = nib.load(cifti)
    header_data = cifti_file.header.get_axis(1)
    for brain_structure in header_data.iter_structures():
        if brain_structure[0] == 'CIFTI_STRUCTURE_{}_LEFT'.format(structure.upper()):
            start = brain_structure[1].start
        if brain_structure[0] == 'CIFTI_STRUCTURE_{}_RIGHT'.format(structure.upper()):
            stop = brain_structure[1].stop
    data = cifti_file.get_fdata()[:,start:stop]
    return np.transpose(data)

def generate_correlation_map(x, y):
    """Correlate each n with each m.

    Parameters
    ----------
    x : np.array
      Shape N X T.

    y : np.array
      Shape M X T.

    Returns
    -------
    np.array
      N X M array in which each element is a correlation coefficient.

    """
    mu_x = x.mean(1)
    mu_y = y.mean(1)
    n = x.shape[1]
    if n != y.shape[1]:
        raise ValueError('x and y must ' +
                         'have the same number of timepoints.')
    s_x = x.std(1, ddof=n - 1)
    s_y = y.std(1, ddof=n - 1)
    cov = np.dot(x,
                 y.T) - n * np.dot(mu_x[:, np.newaxis],
                                  mu_y[np.newaxis, :])
    return cov / np.dot(s_x[:, np.newaxis], s_y[np.newaxis, :])

# Define some variables for testing
subject = snakemake.wildcards.subject
hemi    = snakemake.wildcards.hemi
density = snakemake.wildcards.density

# Load hippocampal rfMRI data
rfmri_hipp_file = snakemake.input.rfmri_hipp
rfmri_hipp_gii  = nib.load(rfmri_hipp_file)

rfmri_hipp_data_rest = np.zeros((len(rfmri_hipp_gii.darrays[0].data),len(rfmri_hipp_gii.darrays)))

for i in range(0,len(rfmri_hipp_gii.darrays)):
    rfmri_hipp_data_rest[:,i] = rfmri_hipp_gii.darrays[i].data

# for r, run in enumerate(snakemake.input.rfmri_ctx):
rfmri_ctx_file  = snakemake.input.rfmri_ctx
rfmri_ctx_data_rest = pull_data(rfmri_ctx_file,'cortex')
# if r == 0:
#     rfmri_ctx_data_rest = pull_data(rfmri_ctx_file,'cortex')
# else:
#     rfmri_ctx_data_rest = np.hstack((rfmri_ctx_data_rest,pull_data(rfmri_ctx_file,'cortex')))
        
# Compute hipp vertex-wise correlation matrix first
print(rfmri_hipp_data_rest.shape)
print(rfmri_ctx_data_rest.shape)

correlation_matrix = generate_correlation_map(rfmri_hipp_data_rest,rfmri_ctx_data_rest)

# Save to npy file
np.save(snakemake.output.correlation_matrix, correlation_matrix)

# Transform correlation matrix to cosine similarity and then normalized angle matrix
dist_upper        = pdist(correlation_matrix,'cosine')
cosine_similarity = 1-squareform(dist_upper)
cosine_similarity = np.nan_to_num(cosine_similarity)
norm_angle_matrix = 1-(np.arccos(cosine_similarity)/math.pi)

# Save to npy file
np.save(snakemake.output.affinity_matrix, norm_angle_matrix)

# # Calculate gradients based on normalized angle matrix
# # Kernel = none as input matrix is already affinity matrix
# gm = GradientMaps(n_components=snakemake.params.n_gradients, kernel=None, random_state=0)
# gm.fit(norm_angle_matrix, diffusion_time=0)

# # Save gradients to gifti file
# for g in range(0,len(gm.gradients_.T)):
#     gii = nib.gifti.GiftiImage()
#     gii.add_gifti_data_array(
#         nib.gifti.GiftiDataArray(data=gm.gradients_.T[g].astype(np.float32))
#     ) 

#     nib.save(gii, snakemake.output.gradient_maps[g])