import os
import pandas as pd
from pathlib import Path

rule map_rfmri_hippunfold_surface:
    input:
        check = bids(
            root = "work",
            datatype = "anat",
            hemi = "{hemi}",
            den = "{density}",
            suffix = "structure.done",
            **config["subj_wildcards"]
        ),
        # surf = 'deriv/post_hippunfold/sub-{subject}/surf_MNI/sub-{subject}_hemi-{hemi}_space-MNI152NLin2009cAsym_den-{density}_desc-nancorrect_midthickness.surf.gii',
        surf = bids(
            root = "results",
            datatype = "anat",
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            desc = "nancorrect",
            suffix = "midthickness.surf.gii",
            **config["subj_wildcards"]
        ),
        rfmri = join(config['fmri_clean_dir'],'sub-{subject}/run-1/regress/sub-{subject}_run-1_residualised.nii.gz')
    output:
        # rfmri = 'deriv/post_fmriprep/sub-{subject}/surface_maps/sub-{subject}_hemi-{hemi}_space-MNI152NLin2009cAsym_den-{density}_bold.func.gii'
        rfmri = bids(
            root = "results",
            datatype = "func",
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            suffix = "bold.func.gii",
            **config["subj_wildcards"]
        )
    group: 'subj'
    threads: 8
    resources:
        mem_mb = 16000,
        time = 60    
    shell:
        """
        MODULEPATH=/project/6050199/software/transparentsingularity/modules:$MODULEPATH
        export MODULEPATH

        module load connectome-workbench

        wb_command -volume-to-surface-mapping {input.rfmri} {input.surf} {output.rfmri} -trilinear
        """

rule calculate_affinity_matrix:
    input:
        # rfmri_hipp = 'deriv/post_fmriprep/sub-{subject}/surface_maps/sub-{subject}_hemi-{hemi}_space-MNI152NLin2009cAsym_den-{density}_bold.func.gii',
        rfmri_hipp = bids(
            root = "results",
            datatype = "func",
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            suffix = "bold.func.gii",
            **config["subj_wildcards"]
        ),
        rfmri_ctx = join(config['fmri_clean_dir'],'sub-{subject}/run-1/regress/sub-{subject}_run-1_residualised_space-fsLR_den-91k_bold.dtseries.nii'),
    params:
        n_gradients = config['n_gradients']
    output:
        # correlation_matrix = 'deriv/post_fmriprep/sub-{subject}/correlation_matrix/sub-{subject}_hemi-{hemi}_space-MNI152NLin2009cAsym_den-{density}_correlationmatrix.npy',
        correlation_matrix = bids(
            root = "results",
            datatype = "func",
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            suffix = "correlationmatrix.npy",
            **config["subj_wildcards"]
        ),
        # affinity_matrix = 'deriv/post_fmriprep/sub-{subject}/correlation_matrix/sub-{subject}_hemi-{hemi}_space-MNI152NLin2009cAsym_den-{density}_affinitymatrix.npy',
        affinity_matrix = bids(
            root = "results",
            datatype = "func",
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            suffix = "affinitymatrix.npy",
            **config["subj_wildcards"]
        )
    threads: 8
    resources:
        mem_mb = 16000,
        time = 30
    group: 'subj'
    script: '../scripts/calculate_affinity_matrix.py'

# rule concat_subj_matrix_tsv:
#     input:
#         affinity_matrix = bids(
#             root = "results",
#             datatype = "func",
#             hemi = "L",
#             space = "MNI152NLin2009cAsym",
#             den = "2k",
#             suffix = "affinitymatrix.npy",
#             **config["subj_wildcards"]
#         )
#     params:
#         tsv = bids(
#             root = "work",
#             prefix = "group",
#             suffix = "complete.tsv"
#         ),
#     output:
#         tsv_done = bids(
#             root = "work",
#             suffix = "aggregate.done",
#             **config["subj_wildcards"]
#         )
#     group: "subj"
#     run:
#         print()
#         if os.path.exists(params.tsv) == False:
#             pd.DataFrame.from_dict({'participant_id':[wildcards.subject]}).to_csv(params.tsv,sep='\t',index=False)
#         else:
#            df = pd.read_csv(params.tsv,sep='\t')
#            df_newsub = pd.DataFrame.from_dict({'participant_id':[wildcards.subject]})
#            df.append(df_newsub).to_csv(params.tsv,sep='\t',index=0)
#            Path(output.tsv_done).touch()

rule calculate_average_gradients:
    input:
        # affinity_matrix = expand('deriv/post_hcp/sub-{subject}/correlation_matrix/sub-{subject}_hemi-{{hemi}}_space-MNI_den-{{density}}_affinitymatrix.npy',subject=config['participant_label'])
        affinity_matrix = expand(bids(
            root = "results",
            datatype = "func",
            hemi = "{{hemi}}",
            space = "MNI152NLin2009cAsym",
            den = "{{density}}",
            suffix = "affinitymatrix.npy",
            **config["subj_wildcards"]),
            subject = config['input_lists']['reverse_transform']['subject'],
        )
        # subject_info = 'resources/subject_info.csv'
    output:
        # gradient_maps = 'deriv/post_hcp/sub-avg/surface_maps/sub-avg_hemi-{hemi}_space-MNI_den-{density}_gradients.func.gii',
        # affinity_matrix = 'deriv/post_hcp/sub-avg/correlation_matrix/sub-avg_hemi-{hemi}_space-MNI_den-{density}_affinitymatrix.npy'
        gradient_maps = bids(
            root = "results",
            datatype = "group",
            prefix = 'sub_avg',
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            suffix = "gradients.func.gii"),
        
        affinity_matrix = bids(
            root = "results",
            datatype = "group",
            prefix = "sub_avg",
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            suffix = "affinitymatrix.npy")

    params:
        n_gradients = config['n_gradients']
    group: 'calc_gradients'
    script: '../scripts/calculate_average_gradients.py'

rule calculate_aligned_gradients:
    input:
        # affinity_matrix = 'deriv/post_fmriprep/sub-{subject}/correlation_matrix/sub-{subject}_hemi-{hemi}_space-MNI_den-{density}_affinitymatrix.npy',
        affinity_matrix = bids(
            root = "results",
            datatype = "func",
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            suffix = "affinitymatrix.npy",
            **config["subj_wildcards"]),
        # avg_affinity_matrix = 'deriv/post_fmriprep/sub-avg/correlation_matrix/sub-avg_hemi-{hemi}_space-MNI_den-{density}_affinitymatrix.npy'
        
        avg_affinity_matrix = bids(
            root = "results",
            datatype = "group",
            prefix = "sub_avg",
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            suffix = "affinitymatrix.npy",
        ),
    output:
        # gradient_maps = 'deriv/post_fmriprep/sub-{subject}/surface_maps/sub-{subject}_hemi-{hemi}_space-MNI_den-{density}_alignedgradients.func.gii'
        gradient_maps = bids(
            root = "results",
            datatype = "func",
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            density = "{density}",
            desc = "aligned",
            suffix = "gradients.func.gii",
            **config["subj_wildcards"]
        )
    params:
        n_gradients = config['n_gradients']
    group: 'calc_gradients'
    script: '../scripts/calculate_aligned_gradients.py'

