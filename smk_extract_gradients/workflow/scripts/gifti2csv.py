import nibabel as nb
import numpy as np

def gifti2csv(gii_file,out_file, itk_lps = True):
        gii = nb.load(gii_file)
        data = gii.get_arrays_from_intent('NIFTI_INTENT_POINTSET')[0].data

        if itk_lps:  # ITK: flip X and Y around 0
            data[:, :2] *= -1

        # antsApplyTransformsToPoints requires 5 cols with headers
        csvdata = np.hstack((data, np.zeros((data.shape[0], 3))))
        
        np.savetxt(
            out_file,
            csvdata,
            delimiter=",",
            header="x,y,z,t,label,comment",
            fmt=["%.5f"] * 4 + ["%d"] * 2,
        )

gifti2csv(snakemake.input.surf,snakemake.output.surf)