get_fmriprep_dir = lambda fmri_path:  [ '/'.join(fmri_path.split('/')[0:idx+1]) for idx, item in enumerate(fmri_path.split('/')) if 'fmriprep' in item][0]

rule gen_cohort:
    input:
        rfmri = lambda wildcards: fmri_img_dict[wildcards.subject]
    params:
        fmriprep_path = get_fmriprep_dir(config['input_path']['bold_volume']) if config['fmriprep_dir'] == None else config['fmriprep_dir'],
        # multiple sessions will become point of error, need new solution for this
        session = config['input_lists']['bold_volume']['session'] if 'session' in config['input_lists']['bold_volume'].keys() else False
    output:
        cohort = bids(
            root = "work",
            datatype = "func",
            task = '{task}',
            suffix = 'cohort.csv',
            **subj_wildcards)
    group: 'subj'
    log: bids(root = 'logs',**subj_wildcards, task = '{task}', suffix = 'gen_cohort.txt')
    script:
        '../scripts/gen_cohort.py'

rule run_xcpengine:
    input:
        rfmri = lambda wildcards: fmri_img_dict[wildcards.subject],
        cohort = rules.gen_cohort.output.cohort,
        pipeline_design = os.path.join(config['snakemake_dir'], config['fmri_cleaning_design']['36p']) # need to make this more modular depending on cleaning design chosen
    params:
        fmriprep_dir = get_fmriprep_dir(config['input_path']['bold_volume']) if config['fmriprep_dir'] == None else config['fmriprep_dir'],
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
    group: 'subj'
    resources:
        mem_mb = 32000,
        time = 1440
    log: bids(root = 'logs',**subj_wildcards, task = '{task}', suffix = 'fmri_cleaning.txt')
    shell:
        """
        xcpEngine -d {input.pipeline_design} -c {input.cohort} -r {params.fmriprep_dir} -i {params.work_dir} -o {params.output_dir}
        touch {output.xcp_done}
        """

rule clean_fmri_reorganize:
    input:
        xcp_done = bids(
            root = "work",
            datatype = "func",
            task = '{task}',
            suffix = "fmriclean.done",
            **subj_wildcards),
        cohort = rules.gen_cohort.output.cohort
    output:
        fmri_volume =  bids(
            root = 'results',
            datatype = 'func',
            task =  '{task}',
            desc =  'cleaned',
            suffix =  'bold.nii.gz',
            **subj_wildcards
            ),
        fmri_surf =  bids(
            root = 'results',
            datatype = 'func',
            task =  '{task}',
            desc =  'cleaned',
            suffix =  'bold.surf.gii',
            **subj_wildcards
            ),
    group: 'subj'
    log: bids(root = 'logs',**subj_wildcards, task = '{task}', suffix = 'clean_fmri_reorganize.txt')
    run:
        cohort_path = fmri_path_cohort(input.cohort)
        fmri_volume_path = join('results/xcpengine/', cohort_path[0])
        fmri_surf_path = join('results/xcpengine/', cohort_path[1])
        shutil.copyfile(fmri_volume_path, output.fmri_volume)
        shutil.copyfile(fmri_surf_path, output.fmri_surf)