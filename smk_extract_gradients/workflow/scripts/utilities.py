import nibabel as nb
import numpy as np
import pandas as pd
import os

def gifti2csv(gii_file, out_file, itk_lps = True):
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

def get_fmriprep_dir(fmri_path):
    fmri_path_list = fmri_path.split('/')
    
    for idx, path_item in enumerate(fmri_path_list):
        if 'fmriprep' == path_item.strip():
            return '/'.join(fmri_path_list[0:idx])

def gen_cohort(fmri_path,fmriprep_dir,output_file,wildcards):
    
    cohort_dict = {}
    
    for idx, wildcard in enumerate(wildcards.values()):
        if 'subject' in wildcard:
            cohort_col = len(cohort_dict.values())
            cohort_dict['id'+str(cohort_col)] = [wildcards['subject']]
        elif 'ses' in wildcard:
            cohort_col = len(cohort_col) + 1
            cohort_dict['id'+str(cohort_col)] = [wildcards['ses']]
        elif 'run' in wildcard:
            cohort_col = len(cohort_col) + 1
            cohort_dict['id'+str(cohort_col)] = [wildcards['run']]
    
    cohort_dict['img'] = fmri_path.replace(fmriprep_dir,'').strip('/')

    cohort_df = pd.Dataframe.from_dict(cohort_dict)

    if os.path.exists(cohort_df) == False:
        cohort_df.to_csv(output_file, index=0)
    else:
        old_cohort_df =  pd.read_csv(output_file)
        pd.concat([old_cohort_df,cohort_df]).to_csv(output_file, index = 0)

