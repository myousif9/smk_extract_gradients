get_fmriprep_dir = lambda fmri_path:  [ '/'.join(fmri_path.split('/')[0:idx+1]) for idx, item in enumerate(fmri_path.split('/')) if 'fmriprep' in item][0]

rule gen_cohort:
    input:
        rfmri = lambda wildcards: fmri_img_list[wildcards.subject]
    params:
        fmriprep_path = get_fmriprep_dir(config['input_path']['bold_volume']),
        # subject = lambda wildcards: wildcards.subject,
        # run = next(iter([ item.replace('run-','') for item in config['input_path']['bold_volume'].split('/')[-1].split('_') if 'run' in item ])),
        # run = dict(x.split('-') for x in config['input_path']['bold_volume'].split('/')[-1].split('_') if len(x.split('-')) == 2 )['run'],
        session = config['input_lists']['bold_volume']['session'] if 'session' in config['input_lists']['bold_volume'].keys() else False
    output:
        cohort = bids(
            root = "work",
            datatype = "func",
            task = '{task}',
            suffix = 'cohort.csv',
            **subj_wildcards)
    group: 'xcpengine_subj'
    log: bids(root = 'logs',**subj_wildcards, task = '{task}', suffix = 'gen_cohort.txt')
    script:
        '../scripts/gen_cohort.py'

rule run_xcpengine:
    input:
        rfmri = lambda wildcards: fmri_img_list[wildcards.subject],
        cohort = bids(
            root = "work",
            datatype = "func",
            task = '{task}',
            suffix = 'cohort.csv',
            **subj_wildcards),
        pipeline_design = os.path.join(config['snakemake_dir'], config['fmri_cleaning_design']['36p']) # need to make this more modular depending on cleaning design chosen
    params:
        fmriprep_dir = get_fmriprep_dir(config['input_path']['bold_volume']),
        work_dir = "work/xcpengine/",
        output_dir = "results/xcpengine/"
    output:
        xcp_done = bids(
            root = "work",
            datatype = "func",
            task = "{task}",
            suffix = "fmriclean.done",
            **subj_wildcards),
    container: config['singularity']['xcpengine']
    group: 'xcpengine_subj'
    resources:
        mem_mb = 32000,
        time = 120
    log: bids(root = 'logs',**subj_wildcards, task = '{task}', suffix = 'fmri_cleaning.txt')
    shell:
        """
        xcpEngine -d {input.pipeline_design} -c {input.cohort} -r {params.fmriprep_dir} -i {params.work_dir} -o {params.output_dir}
        touch {output.xcp_done}
        """

rule clean_fmri_reorganize:
    input:
        xcp_done = expand(bids(
            root = "work",
            datatype = "func",
            task = '{{task}}',
            suffix = "fmriclean.done",
            **subj_wildcards),
            subject = subjects,
            ),
        cohort = bids(
            root = "work",
            datatype = "func",
            task = '{task}',
            suffix = 'cohort.csv',
            **subj_wildcards)
    output:
        fmri =  bids(
            root = 'results',
            datatype = 'func',
            task =  '{task}',
            desc =  'cleaned',
            suffix =  'bold.nii.gz',
            **subj_wildcards),
    group: 'xcpengine_group'
    log: bids(root = 'logs',**subj_wildcards, task = '{task}', suffix = 'clean_fmri_reorganize.txt')
    run:
        input_fmri_path = join('results/xcpengine/', fmri_path_cohort({input.cohort})[0])
        shutil.copyfile(input_fmri_path, {output.fmri})

