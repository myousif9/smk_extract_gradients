rule map_rfmri_hippunfold_surface:
    input:
        surf = 'deriv/post_hippunfold/sub-{subject}/surf_MNI/sub-{subject}_hemi-{hemi}_space-MNI_den-{density}_desc-nancorrect_midthickness.surf.gii',
        rfmri = join(config['fmri_clean_dir'],'sub-{subject}/run-1/regress/sub-{subject}_run-1_residualised.nii.gz')
    output:
        rfmri = 'deriv/post_fmriprep/sub-{subject}/surface_maps/sub-{subject}_hemi-{hemi}_space-MNI_den-{density}_fMRI_CONCAT_ALL.func.gii'
    group: 'subj'
    threads: 8
    resources:
        mem_mb = 16000,
        time = 60    
    shell:
        """
        wb_command -volume-to-surface-mapping {input.rfmri} {input.surf} {output.rfmri} -trilinear
        """

rule calculate_gradients:
    input:
        rfmri_hipp = 'deriv/post_fmriprep/sub-{subject}/surface_maps/sub-{subject}_hemi-{hemi}_space-MNI_den-{density}_fMRI_CONCAT_ALL.func.gii',
        rfmri_ctx = join(config['hcp']['func'],'sub-{subject}/MNINonLinear/Results/rfMRI_REST/rfMRI_REST_Atlas_MSMAll_hp0_clean.dtseries.nii'),
        rfmri_info = join(config['hcp']['func_ext'],'sub-{subject}/MNINonLinear/Results/fMRI_CONCAT_ALL/fMRI_CONCAT_ALL_Runs.csv') 
    params:
        n_gradients = config['n_gradients']
    output:
        gradient_maps = expand('deriv/post_fmriprep/sub-{{subject}}/surface_maps/sub-{{subject}}_hemi-{{hemi}}_space-MNI_den-{{density}}_gradient-{g}.func.gii',g=range(1,config['n_gradients']+1)),
        correlation_matrix = 'deriv/post_fmriprep/sub-{subject}/correlation_matrix/sub-{subject}_hemi-{hemi}_space-MNI_den-{density}_correlationmatrix.npy',
        affinity_matrix = 'deriv/post_fmriprep/sub-{subject}/correlation_matrix/sub-{subject}_hemi-{hemi}_space-MNI_den-{density}_affinitymatrix.npy'
    threads: 8
    resources:
        mem_mb = 16000,
        time = 30
    group: 'subj'
    script: '../scripts/calculate_connectivity_matrix.py'

rule calculate_average_gradients:
    input:
        affinity_matrix = expand('deriv/post_hcp/sub-{subject}/correlation_matrix/sub-{subject}_hemi-{{hemi}}_space-MNI_den-{{density}}_affinitymatrix.npy',subject=subjects)
        subject_info = 'resources/subject_info.csv'
    output:
        gradient_maps = 'deriv/post_hcp/sub-avg/surface_maps/sub-avg_hemi-{hemi}_space-MNI_den-{density}_gradients.func.gii',
        affinity_matrix = 'deriv/post_hcp/sub-avg/correlation_matrix/sub-avg_hemi-{hemi}_space-MNI_den-{density}_affinitymatrix.npy'
    params:
        n_gradients = 10
    group: 'group'
    script: '../scripts/calculate_average_gradients.py'

# rule map_myelin_hippunfold_surface:
#     input:
#         surf = 'deriv/hippunfold/{subject}_V1_MR/surf_T1w/sub-{subject}_hemi-{hemi}_space-T1w_den-{density}_midthickness.surf.gii',
#         rfmri = 'deriv/HCP/{subject}_V1_MR/T1w/T1wDividedByT2w.nii.gz'
#     output:
#         rfmri = 'deriv/myelin/{subject}_V1_MR/gifti_map/sub-{subject}_hemi-{hemi}_space-MNI_den-{density}_myelin.func.gii'
#     shell:
#         """
#         wb_command -volume-to-surface-mapping {input.rfmri} {input.surf} {output.rfmri} -trilinear
#         """
