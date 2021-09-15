import nibabel as nb
import numpy as np

def csv2gifti(csv_file, gii_file, out_file, itk_lps = True):
    gii = nb.load(gii_file)
    vertices = gii.get_arrays_from_intent('NIFTI_INTENT_POINTSET')[0].data
    faces = gii.get_arrays_from_intent('NIFTI_INTENT_TRIANGLE')[0].data 
    
    data = np.loadtxt(
        csv_file, delimiter=",", skiprows=1, usecols=(0, 1, 2)
    )

    if itk_lps:  # ITK: flip X and Y around 0
        data[:, :2] *= -1
    
    new_gii = nb.gifti.GiftiImage(header=gii.header,meta = gii.meta)

    new_gii.add_gifti_data_array(
        nb.gifti.GiftiDataArray(
            data=data[:, :3].astype(vertices.dtype),
            intent='NIFTI_INTENT_POINTSET'
        )
    )
    
    new_gii.add_gifti_data_array(
        nb.gifti.GiftiDataArray(
                data = faces.astype(faces.dtype),
                intent = 'NIFTI_INTENT_TRIANGLE'
            )
    )
    
    new_gii.to_filename(out_file)

csv2gifti(snakemake.input.surf_csv,snakemake.input.surf_gii,snakemake.ouptut.surf)