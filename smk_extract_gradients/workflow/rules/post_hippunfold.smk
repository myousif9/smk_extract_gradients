rule correct_nan_vertices:
    input: 
        surf = lambda wildcards: hippunfold_surf_list[wildcards.subject]
        # surf = bids(
        #     root = join(config['hippunfold_dir'],'results'),
        #     datatype = "surf_T1w",
        #     hemi = "{hemi}",
        #     space = "T1w",
        #     den = "{density}",
        #     suffix = "midthickness.surf.gii",
        #     **subj_wildcards
        # ),
    output:
        surf = bids(
            root = "results",
            datatype = "anat",
            hemi = "{hemi}",
            space = "T1w",
            den = "{density}",
            desc = "nancorrect",
            suffix = "midthickness.surf.gii",
            **subj_wildcards
        ),
    group: 'struct_subj'
    log: bids(root = 'logs',**subj_wildcards,  hemi = '{hemi}', den = '{density}', suffix = 'correct_nan_verticies.txt')
    script: '../scripts/fix_nan_vertices.py'

rule gifti2csv:
    input:
        surf = bids(
            root = "results",
            datatype = "anat",
            hemi = "{hemi}",
            space = "T1w",
            den = "{density}",
            desc = "nancorrect",
            suffix = "midthickness.surf.gii",
            **subj_wildcards
        ),
    output:
        surf = bids(
            root = "work",
            datatype = "anat",
            hemi = "{hemi}",
            space = "T1w",
            den = "{density}",
            desc = "nancorrect",
            suffix = "midthickness.surfpoints.csv",
            **subj_wildcards
        ),
    group: 'struct_subj'
    log: bids(root = 'logs',**subj_wildcards, hemi = '{hemi}', den = '{density}', suffix = 'gifti2csv.txt')
    run:
        gifti2csv(input.surf,output.surf)

rule apply_transform:
    input:
        rvr_transform = config['input_path']['reverse_transform'],
        surf = bids(
            root = "work",
            datatype = "anat",
            hemi = "{hemi}",
            space = "T1w",
            den = "{density}",
            desc = "nancorrect",
            suffix = "midthickness.surfpoints.csv",
            **subj_wildcards
        ),
    output:
        surf = bids(
            root = "work",
            datatype = "anat",
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            desc = "nancorrect",
            suffix = "midthickness.surfpoints.csv",
            **subj_wildcards
        ),
    container: config['singularity']['ants']
    group: 'struct_subj'
    log: bids(root = 'logs',**subj_wildcards, hemi = '{hemi}', den = '{density}', suffix = 'apply_transform.txt')
    shell:
        """ 
        # replace with container
        # module load StdEnv/2020
        # module load gcc/9.3.0 
        # module load ants/2.3.5

        antsApplyTransformsToPoints -d 3 -i {input.surf} -o {output.surf} -t {input.rvr_transform}
        """

rule csv2gifti:
    input:
        surf_csv = bids(
            root = "work",
            datatype = "anat",
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            desc = "nancorrect",
            suffix = "midthickness.surfpoints.csv",
            **subj_wildcards
        ),
        surf_gii = bids(
            root = "results",
            datatype = "anat",
            hemi = "{hemi}",
            space = "T1w",
            den = "{density}",
            desc = "nancorrect",
            suffix = "midthickness.surf.gii",
            **subj_wildcards
        ),
    output:
        surf = bids(
            root = "results",
            datatype = "anat",
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            desc = "nancorrect",
            suffix = "midthickness.surf.gii",
            **subj_wildcards
        ),
    group: 'struct_subj'
    log: bids(root = 'logs',**subj_wildcards, hemi = '{hemi}', den = '{density}', suffix = 'csv2gifti.txt')
    run:
        csv2gifti(input.surf_csv, input.surf_gii, output.surf)

rule set_surf_structure:
    input:
        surf = bids(
            root = "results",
            datatype = "anat",
            hemi = "{hemi}",
            space = "MNI152NLin2009cAsym",
            den = "{density}",
            desc = "nancorrect",
            suffix = "midthickness.surf.gii",
            **subj_wildcards
        ),
    params:
        structure = "CORTEX_LEFT" if "{hemi}" == "L" else "CORTEX_RIGHT"
    output:
        check = bids(
            root = "work",
            datatype = "anat",
            hemi = "{hemi}",
            den = "{density}",
            suffix = "structure.done",
            **subj_wildcards
        ),
    container: config['singularity']['autotop']
    group: 'struct_subj'
    log: bids(root = 'logs',**subj_wildcards, hemi = '{hemi}', den = '{density}', suffix = 'set_surf_structure.txt')
    shell: 
        """
        # replace with container
        # MODULEPATH=/project/6050199/software/transparentsingularity/modules:$MODULEPATH
        # export MODULEPATH

        # module load connectome-workbench
        wb_command -set-structure {input.surf} {params.structure} -surface-type ANATOMICAL
        touch {output.check}
        """


