rule correct_nan_vertices:
    input: 
        surf = lambda wildcards: hippunfold_surf_list[wildcards.subject]
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
    group: 'subj'
    script: '../scripts/fix_nan_vertices.py'

rule gifti2csv:
    input:
        surf = rules.correct_nan_vertices.output.surf
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
    group: 'subj'
    run:
        gifti2csv(input.surf,output.surf)

rule apply_transform:
    input:
        rvr_transform = config['input_path']['reverse_transform'],
        surf = rules.gifti2csv.output.surf
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
    group: 'subj'
    shell:
        """ 
        antsApplyTransformsToPoints -d 3 -i {input.surf} -o {output.surf} -t {input.rvr_transform}
        """

rule csv2gifti:
    input:
        surf_csv = rules.apply_transform.output.surf,
        surf_gii = rules.correct_nan_vertices.output.surf
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
    group: 'subj'
    run:
        csv2gifti(input.surf_csv, input.surf_gii, output.surf)

rule set_surf_structure:
    input:
        surf = rules.csv2gifti.output.surf
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
    group: 'subj'
    shell: 
        """
        wb_command -set-structure {input.surf} {params.structure} -surface-type ANATOMICAL
        touch {output.check}
        """


