import nibabel as nib
import numpy as np
from brainspace.gradient import GradientMaps
from scipy.spatial.distance import pdist, squareform

def compute_affinity_matrix(correlation_matrix):
  
  correlation_matrix = np.nan_to_num(correlation_matrix)

  # Transform correlation matrix to cosine similarity and then normalized angle matrix
  dist_upper        = pdist(correlation_matrix,'cosine')
  cosine_similarity = 1-squareform(dist_upper)
  cosine_similarity = np.nan_to_num(cosine_similarity)
  
  return cosine_similarity
#   norm_angle_matrix = 1-(np.arccos(cosine_similarity)/math.pi)

if __name__ == '__main__':
    affinity_mat = snakemake.input.affinity_matrix
    avg_corr_out = snakemake.output.avg_correlation_matrix
    n_gradients = snakemake.params.n_gradients
    hemi = snakemake.wildcards.hemi

    # Get number of subjects
    n_subjects = len(affinity_mat)

    # Load first subject's correlation matrix for shape
    affinity_matrix_data = np.load(affinity_mat[0])

    # Create empty matrix to fill
    affinity_matrix_concat = np.zeros((
        affinity_matrix_data.shape[0],
        affinity_matrix_data.shape[1],
        n_subjects)
    )

    # Fill empty matrix using subjects' correlation matrices
    for s in range(0,n_subjects):
        affinity_matrix_concat[:,:,s] = np.nan_to_num(np.load(affinity_mat[s]))

    # Average correlation matrices
    affinity_matrix_avg = np.mean(affinity_matrix_concat,2)

    # Save matrix
    np.save(avg_corr_out, affinity_matrix_avg)
    
    affinity_reference =  compute_affinity_matrix(affinity_matrix_avg)

    # Calculate average gradients
    gm = GradientMaps(n_components=n_gradients, kernel=None, random_state=0)
    gm.fit(affinity_reference, diffusion_time=0)

    # Save gradients to gifti file
    gii = nib.gifti.GiftiImage()

    for g in range(0,len(gm.gradients_.T)):
        gii.add_gifti_data_array(
            nib.gifti.GiftiDataArray(
                data=gm.gradients_.T[g].astype(np.float32),
                meta={
                    'AnatomicalStructurePrimary':'CortexLeft' if hemi == 'L' else 'CortexRight',
                    'Name':'Gradient {}'.format(g+1)
                    }
                )
        ) 

    nib.save(gii, snakemake.output.avg_gradient_maps)
    
    for s in range(0,n_subjects):
        sub_affinity_mat = affinity_matrix_concat[:,:,s]

        gp = GradientMaps(n_components=n_gradients, kernel=None, alignment='procrustes', random_state=0)
        gp.fit([affinity_reference,sub_affinity_mat], diffusion_time=0)

        # Save aligned gradients to gifti file
        gii = nib.gifti.GiftiImage()
        gii_unaligned = nib.gifti.GiftiImage()

        for g in range(0,n_gradients):
            gii.add_gifti_data_array(
                nib.gifti.GiftiDataArray(
                    data=gp.aligned_[1][:,g].astype(np.float32),
                    meta={
                        'AnatomicalStructurePrimary':'CortexLeft' if hemi == 'L' else 'CortexRight',
                        'Name':'Gradient {}'.format(g+1)
                        }
                    )
            ) 

            gii_unaligned.add_gifti_data_array(
                nib.gifti.GiftiDataArray(
                    data=gp.gradients_[1][:,g].astype(np.float32),
                    meta={
                        'AnatomicalStructurePrimary':'CortexLeft' if hemi == 'L' else 'CortexRight',
                        'Name':'Gradient {}'.format(g+1)
                        }
                    )
            ) 

        nib.save(gii, snakemake.output.gradient_maps[s])
        nib.save(gii_unaligned, snakemake.output.gradient_unaligned[s])


