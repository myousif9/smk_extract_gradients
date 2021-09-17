from os.path import join

# This script transforms hippunfold midthickness surfaces into micapipe's space-rsfmri. 

# TODO:
# I think ideally the rsfmri and midthickness surfaces should BOTH be transformed into space-MNI152 (but the rsfmri should not be ubpsampled if possible)
# temporary file names should maybe be updated to something more systematic
# session (eg. ses-pre-01) should be made a wildcard
# a flag could be added (eg. `--is_micapipe`) to use this file/folder struture instead of the normal bids grabber

rule register_t1w_to_nativepro:
    input:
        t1w = join(config['hippunfold'],'sub-{subject}/anat/sub-{subject}_desc-preproc_T1w.nii.gz'),
        nativepro = join(config['micapipe'], 'sub-{subject}/ses-pre/anat/sub-{subject}_ses-pre_space-nativepro_t1w.nii.gz'),
    output:
        reg = join(config['deriv'],'sub-{subject}/t1w_to_nativepro.nii.gz'),
        xfm = join(config['deriv'],'sub-{subject}/t1w_to_nativepro_xfm.txt'),
    shell: 'singularity exec /data/mica3/jordand/singularity/hippunfold_v0.6.0.sif reg_aladin -flo {input.t1w} -ref {input.nativepro} -res {output.reg} -aff {output.xfm} -rigOnly'

rule xfm_itk2ras:
    input:
        xfm_itk = join(config['deriv'],'sub-{subject}/t1w_to_nativepro_xfm.txt'),
    output:
        xfm_ras = join(config['deriv'],'sub-{subject}/t1w_to_nativepro_xfm_ras.txt'),
    group: 'subj'
    shell: 'c3d_affine_tool {input.xfm_itk} -inv -o {output.xfm_ras}'

rule transform_t1w_to_nativepro:
    input:
        surf = join(config['hippunfold'],'sub-{subject}/surf_T1w/sub-{subject}_hemi-{hemi}_space-T1w_den-{density}_midthickness.surf.gii'),
        xfm_ras = join(config['deriv'],'sub-{subject}/t1w_to_nativepro_xfm_ras.txt'),
    output:
        surf = join(config['deriv'],'sub-{subject}/tmp_hemi-{hemi}_den-{density}_space-nativepro_midthickness.surf.gii')
    group: 'subj'  
    shell: 'wb_command -surface-apply-affine {input.surf} {input.xfm_ras} {output.surf}'



rule template_xfm_itk2ras1:
    input:
        xfm_itk = join(config['micapipe'], 'sub-{subject}/ses-pre/xfm/sub-{subject}_ses-pre_rsfmri_from-rsfmri_to-nativepro_mode-image_desc-affine_0GenericAffine.mat'),
    output:
        xfm_ras = join(config['deriv'], 'sub-{subject}/tmp1.txt')
    group: 'subj'
    shell: 'c3d_affine_tool -itk {input.xfm_itk} -o {output.xfm_ras}'

rule transform_hippunfold_surface1:
    input:
        surf = join(config['deriv'],'sub-{subject}/tmp_hemi-{hemi}_den-{density}_space-nativepro_midthickness.surf.gii'),
        xfm_ras = join(config['deriv'], 'sub-{subject}/tmp1.txt')
    output:
        surf = join(config['deriv'],'sub-{subject}/tmp1_hemi-{hemi}_den-{density}_space-nativepro_midthickness.surf.gii')
    group: 'subj'  
    shell: 'wb_command -surface-apply-affine {input.surf} {input.xfm_ras} {output.surf}'



rule template_xfm_itk2ras2:
    input:
        xfm_itk = join(config['micapipe'], 'sub-{subject}/ses-pre/xfm/sub-{subject}_ses-pre_rsfmri_from-nativepro_rsfmri_to-rsfmri_mode-image_desc-SyN_0GenericAffine.mat'),
    output:
        xfm_ras = join(config['deriv'], 'sub-{subject}/tmp2.txt')
    group: 'subj'
    shell: 'c3d_affine_tool -itk {input.xfm_itk} -inv -o {output.xfm_ras}'

rule transform_hippunfold_surface2:
    input:
        surf = join(config['deriv'],'sub-{subject}/tmp1_hemi-{hemi}_den-{density}_space-nativepro_midthickness.surf.gii'),
        xfm_ras = join(config['deriv'], 'sub-{subject}/tmp2.txt')
    output:
        surf = join(config['deriv'],'sub-{subject}/tmp2_hemi-{hemi}_den-{density}_midthickness.surf.gii')
    group: 'subj'  
    shell: 'wb_command -surface-apply-affine {input.surf} {input.xfm_ras} {output.surf}'



rule transform_hippunfold_surface_deform:
    input:
        surf = join(config['deriv'],'sub-{subject}/tmp2_hemi-{hemi}_den-{density}_midthickness.surf.gii'),
        warp = join(config['micapipe'], 'sub-{subject}/ses-pre/xfm/sub-{subject}_ses-pre_rsfmri_from-nativepro_rsfmri_to-rsfmri_mode-image_desc-SyN_1InverseWarp.nii.gz'),
    output:
        surf = join(config['deriv'],'sub-{subject}/surf_T1w/sub-{subject}_hemi-{hemi}_space-T1w_den-{density}_midthickness_nonNaNcorrect.surf.gii')
    group: 'subj'  
    shell: 'wb_command -surface-apply-warpfield {input.surf} {input.warp} {output.surf}'


