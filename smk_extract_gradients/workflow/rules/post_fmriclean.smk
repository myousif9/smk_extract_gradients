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
    log: bids(root = 'logs',**subj_wildcards, task = '{task}', hemi = '{hemi}', den = '{density}', suffix = 'map-rfmri-hippunfold-surface.txt')
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
    log: bids(root = 'logs',**subj_wildcards, task = '{task}', hemi = '{hemi}', den = '{density}', suffix = 'calculate-affinity-matrix.txt')
    script: '../scripts/calculate_affinity_matrix.py'

# compile affinity matrices that exist and 

# affinity_path =bids(
#     root = "results",
#     datatype = "func",
#     task = '{task}',
#     hemi = "{hemi}",
#     space = "MNI152NLin2009cAsym",
#     den = "{density}",
#     suffix = "affinitymatrix.npy",
#     **subj_wildcards)

# affinity_mat_subjects = {}

# for idx, subj in enumerate(fmri_input_list['subject']):
#     for hemi in config['hemi']:
#         if os.path.exists(affinity_path.format(subject=subj, task=config['task'],hemi=hemi, density=config['density'])):
#             affinity_exist_subjects.add(subj)
#         else:
#             affinity_exist_subjects.discard(subj)

# affinity_mat_subjects = list(affinity_mat_subjects)
            
rule calculate_average_gradients:
    input:
        affinity_matrix = expand(bids(
            root = "results",
            datatype = "func",
            task = '{{task}}',
            hemi = "{{hemi}}",
            space = "MNI152NLin2009cAsym",
            den = "{{density}}",
            suffix = "affinitymatrix.npy",
            **subj_wildcards),
            subject = fmri_input_list['subject']
            ),
    output:
        gradient_maps = bids(
            root = "results",
            datatype = "group",
            prefix = 'sub_avg',
            task = '{task}',
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            suffix = "gradients.func.gii"
            ),
        affinity_matrix = bids(
            root = "results",
            datatype = "group",
            prefix = "sub_avg",
            task = '{task}',
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            suffix = "affinitymatrix.npy"
            ),
    params:
        n_gradients = config['n_gradients']
    group: 'calc_gradients'
    log: bids(root = 'logs', sub = 'avg', task = '{task}', hemi = '{hemi}', den = '{density}', suffix = 'calculate-average-gradients.txt')
    script: '../scripts/calculate_average_gradients.py'

rule calculate_aligned_gradients:
    input:
        affinity_matrix = rules.calculate_affinity_matrix.output.affinity_matrix,
        avg_affinity_matrix =  rules.calculate_average_gradients.output.affinity_matrix if config['reference_gradient'] == None else config['reference_gradient'],
    output:
        gradient_maps = bids(
            root = "results",
            datatype = "func",
            task = '{task}',
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            desc = "aligned",
            suffix = "gradients.func.gii",
            **subj_wildcards
            ),
        gradient_unaligned = bids(
            root = "results",
            datatype = "func",
            task =  '{task}',
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            desc = "unaligned",
            suffix = "gradients.func.gii",
            **subj_wildcards
            ),
    params:
        n_gradients = config['n_gradients']
    group: 'calc_gradients'
    log: bids(root = 'logs',**subj_wildcards, task = '{task}', hemi = '{hemi}', den = '{density}', suffix = 'calculate-aligned-gradients.txt')
    script: '../scripts/calculate_aligned_gradients.py'

rule set_func_structure:
    input:
        gradient_maps = rules.calculate_aligned_gradients.output.gradient_maps,
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
    log: bids(root = 'logs',**subj_wildcards, task = '{task}', hemi = '{hemi}', den = '{density}', suffix = 'set-func-structure.txt')
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