rule map_rfmri_hippunfold_surface:
    input:
        check_struct = rules.set_surf_structure.output.check,
        surf = rules.csv2gifti.output.surf,
        fmri = rules.clean_fmri_reorganize.output.fmri_volume
    output:
        rfmri = bids(
            root = "results",
            datatype = "func",
            task =  '{task}',
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            suffix = "bold.func.gii",
            **subj_wildcards
            )
    container: config['singularity']['autotop']
    group: 'subj'
    threads: 8
    resources:
        mem_mb = 16000,
        time = 60    
    shell:
        """
        wb_command -volume-to-surface-mapping {input.fmri} {input.surf} {output.rfmri} -trilinear
        """

rule calculate_affinity_matrix:
    input:
        rfmri_hipp = rules.map_rfmri_hippunfold_surface.output.rfmri,
        rfmri_ctx = rules.clean_fmri_reorganize.output.fmri_surf
    params:
        n_gradients = config['n_gradients'],
        # rfmri_ctx = lambda wildcards: join('results/xcpengine/', fmri_path_cohort( input.fmri_cohort_path )[1])
    output:
        correlation_matrix = bids(
            root = "results",
            datatype = "func",
            task = '{task}',
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            suffix = "correlationmatrix.npy",
            **subj_wildcards
            ),
        affinity_matrix = bids(
            root = "results",
            datatype = "func",
            task = '{task}',
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            suffix = "affinitymatrix.npy",
            **subj_wildcards
            )
    threads: 8
    resources:
        mem_mb = 16000,
        time = 30
    group: 'subj'
    script: '../scripts/calculate_affinity_matrix.py'

            
rule calculate_gradients:
    input:
        correlation_matrix = expand(bids(
            root = "results",
            datatype = "func",
            task = '{{task}}',
            hemi = "{{hemi}}",
            space = "MNI152NLin2009cAsym",
            den = "{{density}}",
            suffix = "correlationmatrix.npy",
            **subj_wildcards),
            subject = subjects
            ),
    output:
        avg_gradient_maps = bids(
            root = "results",
            datatype = "group",
            prefix = 'sub-avg',
            task = '{task}',
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            suffix = "gradients.func.gii"
            ),
        avg_correlation_matrix = bids(
            root = "results",
            datatype = "group",
            prefix = "sub-avg",
            task = '{task}',
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            suffix = "correlationmatrix.npy"
            ),
        gradient_maps = expand(bids(
            root = "results",
            datatype = "func",
            task = "{{task}}",
            hemi = "{{hemi}}",
            space = "MNI152NLin2009cAsym",
            den = "{{density}}",
            desc = "aligned",
            suffix = "gradients.func.gii",
            **subj_wildcards),
            subject = subjects
            ),
        gradient_unaligned = expand(bids(
            root = "results",
            datatype = "func",
            task =  "{{task}}",
            hemi = "{{hemi}}",
            space = "MNI152NLin2009cAsym",
            den = "{{density}}",
            desc = "unaligned",
            suffix = "gradients.func.gii",
            **subj_wildcards),
            subject = subjects
            )
    params:
        n_gradients = config['n_gradients']
    group: 'calc_gradients'
    script: '../scripts/calculate_gradients.py'


rule set_func_structure:
    input:
        gradient_maps = rules.calculate_gradients.output.gradient_maps,
    params:
        structure = "CORTEX_LEFT" if "{hemi}" == "L" else "CORTEX_RIGHT"
    output:
        check = bids(
            root = "work",
            datatype = "func",
            task =  '{task}',
            hemi = "{hemi}",
            den = "{density}",
            suffix = "func.done",
            **subj_wildcards
            ),
    container: config['singularity']['autotop']
    group: 'calc_gradients'
    shell: 
        """
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
