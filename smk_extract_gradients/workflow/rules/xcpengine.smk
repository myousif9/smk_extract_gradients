from scripts.utilities import get_friprep_dir, gen_cohort

rule gen_cohort:
    input:
        rfmri = config['input_path']['bold_volume'],
    params:
        fmriprep_dir = get_friprep_dir(config['input_path']['bold_volume']) 
    output:
        cohort = bids(
            root = "work",
            run = '{run}',
            suffix = 'cohort.csv',
            **config['subj_wildcards']
        )
    group: 'subj'
    run:
        gen_cohort(input.rfmri, params.fmriprep_dir, output.cohort, wildcards)
        
rule run_xcpengine:
    input:
        rfmri = config['input_path']['bold_volume'],

        cohort = bids(
            root = "work",
            run = '{run}',
            suffix = 'cohort.csv',
            **config['subj_wildcards']
        )
    params:
        fmriprep_dir = get_friprep_dir(config['input_path']['bold_volume']),
        pipeline_design = '../resources/fc-36p.dsn'
    output:
        work_dir = directory(
            bids(
                root = "work",
                dataype = 'xcpengine'
                **config['subj_wildcards']
            )
        ),
        output_dir = directory(
            bids(
                root = "result",
                dataype = 'xcpengine'
                **config['subj_wildcards']
            )
        )
    group: 'subj'
    shell:
        """
        xcpengine-singularity \
        --image \
        -d {param.pipeline_design} \
        -c {input.cohort} \
        -r {params.fmriprep_dir} \
        -i {output.work_dir} \
        -o {output.output_dir}
        """