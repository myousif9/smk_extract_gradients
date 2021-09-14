import os
import subprocess
import shutil
from niworkflows.interfaces.surf import CSVToGifti, GiftiToCSV

def surftransform_gii(gii_surf,output_surf, warp, affine, new_space):
    """Takes gifti surface and applies ants transforms
    """
    # convert gii to csv
    result_GiftiToCSV = GiftiToCSV(in_file=gii_surf,
                                   itk_lps=True).run()
    csv_surf = result_GiftiToCSV.outputs.out_file
    
    shutil.move(csv_surf,'/'.join(output_surf.split('/')[0:-1]))
    csv_surf = '/'.join(output_surf.split('/')[0:-1]+[csv_surf])
    
    csv_surf_transformed = csv_surf.replace('space-T1w','space-'+new_space)
    
    # defining command for transfomation
    ants_cmmd = ["antsApplyTransformsToPoints", 
    "-d", "3", 
    "-i", csv_surf, 
    "-o", csv_surf_transformed, 
    "-t", warp,
    "-t", affine, 
    "-n", "NearestNeighbor"]
    
    # loading ants, will be replaced with ants container
    subprocess.run('module load StdEnv/2020'.split(' ') )
    subprocess.run('module load gcc/9.3.0 
    ants/2.3.5

    # running ants command
    os.system(' '.join(ants_cmmd))

    # convert csv to gii
    result_CSVToGifti = CSVToGifti(in_file=csv_surf_transformed,
                                   gii_file=gii_surf,
                                   itk_lps=True).run()
    gii_surf_transformed = result_CSVToGifti.outputs.out_file
    os.rename(gii_surf_transformed,output_surf)


if __name__ == '__main__':
    surftransform_gii(snakemake.input.surf, snakemake.output.surf, snakemake.input.fwd_warp, snakemake.input.fwd_affine, snakemake.params.new_space)
