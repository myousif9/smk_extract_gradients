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
        rfmri = bids(
            root = "results",
            datatype = "func",
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            suffix = "bold.func.gii",
            **config["subj_wildcards"]
        )
    container: config['singularity']['autotop']
    group: 'subj'
    threads: 8
    resources:
        mem_mb = 16000,
        time = 60    
    shell:
        """
        # replace with container
        # MODULEPATH=/project/6050199/software/transparentsingularity/modules:$MODULEPATH
        # export MODULEPATH

        # module load connectome-workbench

        wb_command -volume-to-surface-mapping {input.rfmri} {input.surf} {output.rfmri} -trilinear
        """

rule calculate_affinity_matrix:
    input:
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
        correlation_matrix = bids(
            root = "results",
            datatype = "func",
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            suffix = "correlationmatrix.npy",
            **config["subj_wildcards"]
        ),
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

rule calculate_average_gradients:
    input:
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
    output:
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
        affinity_matrix = bids(
            root = "results",
            datatype = "func",
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            suffix = "affinitymatrix.npy",
            **config["subj_wildcards"]),

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
        gradient_maps = bids(
            root = "results",
            datatype = "func",
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            density = "{density}",
            desc = "aligned",
            suffix = "gradients.func.gii",
            **config["subj_wildcards"]
        ),
        gradient_unaligned = bids(
            root = "results",
            datatype = "func",
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            density = "{density}",
            desc = "unaligned",
            suffix = "gradients.func.gii",
            **config["subj_wildcards"]
        )
    params:
        n_gradients = config['n_gradients']
    group: 'calc_gradients'
    script: '../scripts/calculate_aligned_gradients.py'

rule set_func_structure:
    input:
        # surf = bids(
        #     root = "results",
        #     datatype = "anat",
        #     hemi = "{hemi}",
        #     space = "MNI152NLin2009cAsym",
        #     den = "{density}",
        #     desc = "nancorrect",
        #     suffix = "midthickness.surf.gii",
        #     **config["subj_wildcards"]
        # ),
        gradient_maps = expand(bids(
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
        structure = "CORTEX_LEFT" if "{hemi}" == "L" else "CORTEX_RIGHT"
    output:
        check = bids(
            root = "work",
            datatype = "func",
            hemi = "{hemi}",
            den = "{density}",
            suffix = "func.done",
            **config['subj_wildcards']
        ),
    container: config['singularity']['autotop']
    group: 'subj'
    shell: 
        """
        # replace with container
        # MODULEPATH=/project/6050199/software/transparentsingularity/modules:$MODULEPATH
        # export MODULEPATH

        # module load connectome-workbench
        wb_command -set-structure {input.gradient_maps} {params.structure} -surface-type ANATOMICAL
        touch {output.check}
        """

# rule create_average_surface:
#     input:
#         expand('deriv/post_hippunfold/sub-{subject}/surf_MNI/sub-{subject}_hemi-{{hemi}}_space-MNI_den-{{density}}_desc-nancorrect_midthickness.surf.gii',subject=subject_list)
#     output:
#         surf = 'deriv/post_hippunfold/sub-avg/surf_MNI/sub-avg_hemi-{hemi}_space-MNI_den-{density}_desc-nancorrect_midthickness.surf.gii'
#     group: 'average'
#     script: '../scripts/create_average_surface.py'