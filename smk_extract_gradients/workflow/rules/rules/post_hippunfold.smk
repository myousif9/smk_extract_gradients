# from niworkflows.interfaces.surf import GiftiToCSV, CSVToGifti
import os

rule correct_nan_vertices:
    input: 
        surf = join(config['hippunfold'],'sub-{subject}/surf_T1w/sub-{subject}_hemi-{hemi}_space-T1w_den-{density}_midthickness.surf.gii')
    output:
        surf = 'deriv/post_hippunfold/sub-{subject}/surf_T1w/sub-{subject}_hemi-{hemi}_space-T1w_den-{density}_desc-nancorrect_midthickness.surf.gii'
    group: 'subj'
    script: '../scripts/fix_nan_vertices.py'

rule decompose_transform:
    input:
        transform = expand(join(config['bids_dir'],'\derivatives\fmriprep\sub-{{subject}}\sub={{subject}}_from-{spaces2space}_mode-image_xfm.h5'),spaces2space=['MNI152NLin2009cAsym_to-T1w','T1w_to-MNI152NLin2009cAsym'])
    params:
        transform_name_prefix = str(input.transform).replace('.h5',''),
        warp = '00_' +  params.transform_name_prefix + '_DisplacementFieldTransform.nii.gz',
        affine = '01_' + params.transform_name_prefix + '_AffineTransform.mat',
        transform_dir = 'deriv/post_fmriprep/sub-{subject}/MNINonLinear/',
    output:
        warp = join(params.transform_dir, params.warp)
        affine = join(params.transform_dir, params.affine)
    envmodules:
        "ants/2.3.5"
    group: 'subj'
    run:
        current_dir = os.getcwd()
        os.chdir(params.transform_dir)
        shell('CompositeTransformUtil --disassemble {input.transform} {params.tranform_name_prefix}')
        os.chdir(current_dir)

# rule gifti_to_csv:
#     input:
#         surf = 'deriv/post_hippunfold/{subject}_V1_MR/surf_T1w/sub-{subject}_hemi-{hemi}_space-T1w_den-{density}_desc-nancorrect_midthickness.surf.gii'
#     output:
#         surf_csv =  'deriv/post_hippunfold/{subject}_V1_MR/surf_T1w/sub-{subject}_hemi-{hemi}_space-T1w_den-{density}_desc-nancorrect_midthickness.surfpoints.csv'
#     group: 'subj'
#     run:
#         gifti2csv = GiftiToCSV()
#         gifti2csv.inputs.in_file = input.surf
#         gifti2csv_run = gifti2csv.run()
#         os.rename(str(gifti2csv.outputs).replace('out_file = ','').strip(),input.surf_csv)


# rule transform_hippunfold_Surface:
#     input:
#         transform = config['input_path']['transforms'],
#         surf_csv =  'deriv/post_hippunfold/{subject}_V1_MR/surf_T1w/sub-{subject}_hemi-{hemi}_space-T1w_den-{density}_desc-nancorrect_midthickness.surfpoints.csv'
#     output:
#         surf_csv = 'deriv/post_hippunfold/{subject}_V1_MR/surf_T1w/sub-{subject}_hemi-{hemi}_space-MNI_den-{density}_desc-nancorrect_midthickness.surfpoints.csv'
#     group: 'subj'
#     envmodules:
#         "ants/2.3.5"
#     shell: 
#         """
#         antsApplyTransfromsToPoints -d 3 -i {input.surf_csv} -o {output.surf_csv} -t {input.transfrom}
#         """

# rule csv_to_gifti:
#     input:
#         surf = 'deriv/post_hippunfold/{subject}_V1_MR/surf_T1w/sub-{subject}_hemi-{hemi}_space-T1w_den-{density}_desc-nancorrect_midthickness.surf.gii',
#         surf_csv = 'deriv/post_hippunfold/{subject}_V1_MR/surf_T1w/sub-{subject}_hemi-{hemi}_space-MNI_den-{density}_desc-nancorrect_midthickness.surfpoints.csv',
#     output:
#         surf = 'deriv/post_hippunfold/{subject}_V1_MR/surf_T1w/sub-{subject}_hemi-{hemi}_space-MNI_den-{density}_desc-nancorrect_midthickness.surf.gii'
#     group: 'subj'
#     run:
#         csv2gifti = CSVToGifti()
#         csv2gifti.inputs.in_file = input.surf_csv
#         csv2gifti.inputs.gii_files = input.surf
#         csv2gifti_run = csv2gifti.run()
#         os.rename(str(csv2gifti.outputs).replace('out_file = ','').strip(),input.surf_csv)

rule transform_hippunfold_Surface:
    input:
        surf = 'deriv/post_hippunfold/sub-{subject}/surf_T1w/sub-{subject}_hemi-{hemi}_space-T1w_den-{density}_desc-nancorrect_midthickness.surf.gii',
        xfm_fwd = 'deriv/post_fmriprep/sub-{subject}/MNINonLinear/00_sub-{subject}_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm_DisplacementFieldTransform.nii.gz'),
        xfm_rev = 'deriv/post_fmriprep/sub-{subject}/MNINonLinear/00_sub-{subject}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm_DisplacementFieldTransform.nii.gz')
    output:
        surf = 'deriv/post_hippunfold/{subject}_V1_MR/surf_MNI/sub-{subject}_hemi-{hemi}_space-MNI_den-{density}_desc-nancorrect_midthickness.surf.gii'
    envmodules: 'connectome-workbench/1.4.1_20191117.simg'
    group: 'subj'  
    shell:
        """
        wb_command -surface-apply-warpfield {input.surf} {input.xfm_rev} {output.surf} -fnirt {input.xfm_fwd}
        """

rule create_average_surface:
    input:
        expand('deriv/post_hippunfold/{subject}_V1_MR/surf_MNI/sub-{subject}_hemi-{{hemi}}_space-MNI_den-{{density}}_desc-nancorrect_midthickness.surf.gii',subject=subjects)
    output:
        surf = 'deriv/post_hippunfold/sub-avg/surf_MNI/sub-avg_hemi-{hemi}_space-MNI_den-{density}_desc-nancorrect_midthickness.surf.gii'
    conda:
        "envs/environment.yml"
    group: 'subj'
    script: '../scripts/create_average_surface.py'
